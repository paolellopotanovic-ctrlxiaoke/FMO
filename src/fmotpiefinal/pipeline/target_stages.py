#!/usr/bin/env python
"""Produce label-blind Galectin3 complex poses and frozen FMO2/PIEDA predictions."""
from __future__ import annotations

import argparse
import os
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
import hashlib
import json
from pathlib import Path
import re
import shutil
import sys
from dataclasses import replace


from rdkit import Chem
import numpy as np
from scipy.spatial import cKDTree

from fmotpiefinal.core.fmo_input import FMODeckOptions, render_fmo2_pieda
from fmotpiefinal.core.fragmentation import make_facio_like_fragments
from fmotpiefinal.core.parser import aggregate_tpie, find_selected_tie, parse_pie
from fmotpiefinal.core.preparation import prepare_protein_pdb2pqr
from fmotpiefinal.core.bound_geometry import load_bound_ligand, shared_pocket
from fmotpiefinal.core.runner import run_gamess, pie_roundoff_tolerance
from fmotpiefinal.core.structure import Protein, read_pdb, read_sdf, residue_key
from fmotpiefinal.core.water_hydrogen_transfer import preserve_source_waters
from fmotpiefinal.core.fmo_soscf import enable_fmo_soscf, verify_fmo_soscf_log
from fmotpiefinal.core.fmo_soscf_onset import soscf_onset_deck


def copy_prepared_receptor(source, destination):
    source, destination = Path(source), Path(destination)
    if destination.exists() and digest(source) != digest(destination):
        raise ValueError("prepared receptor does not match frozen reference")
    if not destination.exists():
        shutil.copyfile(source, destination)
    for suffix in (".components.json", ".components.sdf", ".model.json"):
        original, copied = source.with_suffix(suffix), destination.with_suffix(suffix)
        if original.exists():
            if copied.exists() and digest(original) != digest(copied):
                raise ValueError(f"registered receptor sidecar changed: {suffix}")
            if not copied.exists():
                shutil.copyfile(original, copied)
        elif copied.exists():
            raise ValueError(f"unregistered receptor sidecar at destination: {suffix}")


def select_fixed_pocket(protein: Protein, reference: dict, target: str,
                        names: list[str], cutoff: float) -> Protein:
    keys = reference.get("residue_keys")
    if (reference.get("schema") != "fmotpiefinal-fixed-pocket-membership-v1" or
            reference.get("label_free") is not True or reference.get("target") != target or
            reference.get("names") != names or reference.get("cutoff_angstrom") != cutoff or
            not isinstance(keys, list) or len(keys) != len(set(keys)) or
            not keys or any(not isinstance(key, str) for key in keys)):
        raise ValueError("frozen pocket membership does not match the full target")
    chosen = set(keys)
    present = {residue_key(atom) for atom in protein.atoms}
    if not chosen <= present:
        raise ValueError("frozen pocket residue missing from candidate receptor")
    pocket = Protein(tuple(atom for atom in protein.atoms if residue_key(atom) in chosen))
    if len(pocket.atoms) != reference["n_reference_atoms"]:
        raise ValueError("frozen pocket changed atom count or hydrogen representation")
    return pocket


def validate_retained_component_region(protein_pdb, protein, ligands, cutoff):
    metadata = Path(protein_pdb).with_suffix(".components.json")
    if metadata.exists():
        registered = json.loads(metadata.read_text())
        special_ids = {model["id"] for model in registered["models"]}
    else:
        metadata = Path(protein_pdb).with_suffix(".model.json")
        if not metadata.exists():
            return None
        registered = json.loads(metadata.read_text())
        if registered.get("schema") in {"retained-halide-model-v1", "retained-monatomic-ion-model-v1"}:
            special_ids = set(registered["matched_ion_templates"])
        elif registered.get("schema") == "retained-special-chemistry-model-v1":
            from fmotpiefinal.core.special_chemistry import registered_phosphoresidue_model
            registered_phosphoresidue_model(protein_pdb)
            special_ids = set(registered["matched_special_templates"])
        elif registered.get("schema") == "source-terminal-cap-retention-v1":
            special_ids = {row["prepared_id"] for row in registered["retained"]}
        else:
            raise ValueError("unregistered retained-component model schema")
    if registered.get("prepared_sha256") != digest(protein_pdb):
        raise ValueError("registered special-component fingerprint changed")
    def component_key(atom):
        return f"{atom.chain.strip() or '_'}:{atom.residue_name}:{atom.residue_number}"

    special_points = [atom.xyz for atom in protein.atoms if component_key(atom) in special_ids]
    present_ids = {component_key(atom) for atom in protein.atoms}
    if not special_ids or not special_ids <= present_ids:
        raise ValueError("registered special-component atom inventory is missing")
    ligand_points = [atom.xyz for ligand in ligands for atom in ligand.atoms]
    minimum = float(min(cKDTree(ligand_points).query(special_points)[0]))
    if minimum <= cutoff:
        raise ValueError("retained MM component intersects actual ligand-defined QM region; special QM model required")
    return {"status": "MM_ONLY_COMPONENTS_OUTSIDE_QM_REGION", "component_ids": sorted(special_ids),
            "minimum_contact_angstrom": minimum, "cutoff_angstrom": cutoff,
            "component_model_sha256": digest(metadata), "prepared_sha256": digest(protein_pdb)}


def verify_fmo_deck_heavy_geometry(deck: str, threshold_angstrom: float = 0.5) -> None:
    if deck.count(" $FMOXYZ\n") != 1:
        raise ValueError("FMO deck lacks a unique coordinate block")
    lines = deck.split(" $FMOXYZ\n", 1)[1].split(" $END", 1)[0].splitlines()
    coordinates = []
    indices = []
    for atom_index, line in enumerate(lines, start=1):
        tokens = line.split()
        if len(tokens) != 5:
            raise ValueError("FMO coordinate block is malformed")
        if tokens[0].upper() in ("H", "D"):
            continue
        coordinates.append([float(value) for value in tokens[-3:]])
        indices.append(atom_index)
    if not coordinates or not np.isfinite(coordinates).all():
        raise ValueError("FMO heavy-atom coordinates are invalid")
    points = np.asarray(coordinates)
    overlaps = cKDTree(points).query_pairs(threshold_angstrom)
    if overlaps:
        left, right = min(overlaps, key=lambda pair: float(np.linalg.norm(points[pair[0]] - points[pair[1]])))
        distance = float(np.linalg.norm(points[left] - points[right]))
        raise ValueError(f"FMO heavy-atom/cap overlap at deck atoms {indices[left]} and "
                         f"{indices[right]}: {distance:.3f} A; audit peptide numbering and HOP caps")


def receptor_at_minimized_geometry(prepared: Protein, complex_path: Path) -> Protein:
    """Use matched, validated receptor coordinates from an OpenMM complex."""
    receptor = [atom for atom in read_pdb(complex_path).atoms
                if not (atom.record == "HETATM" and atom.residue_name == "UNK"
                        and atom.chain == "X")]
    if len(receptor) != len(prepared.atoms):
        raise ValueError("minimized receptor atom count differs from PDB2PQR receptor")
    identity = lambda atom: (atom.atom_name, atom.residue_name, atom.residue_number,
                             atom.chain, atom.element)
    if any(identity(left) != identity(right) for left, right in
           zip(prepared.atoms, receptor, strict=True)):
        raise ValueError("minimized receptor atom identity/order differs from prepared receptor")
    return Protein(tuple(replace(atom, xyz=pose.xyz) for atom, pose in
                         zip(prepared.atoms, receptor, strict=True)))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def equivalent_frozen_protocol(existing: dict, requested: dict) -> bool:
    path_fields = ("poses_root", "prepared_protein_override", "ligand_sdf_override",
                   "reference_pocket_root", "reference_water_hydrogens")
    for field in path_fields:
        if field in existing and field in requested:
            if Path(existing[field]).resolve() != Path(requested[field]).resolve():
                return False
    return ({key: value for key, value in existing.items() if key not in path_fields} ==
            {key: value for key, value in requested.items() if key not in path_fields} and
            set(existing) == set(requested))


def directory_key(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", name)


def sanitize(source: Path, output: Path, expected: int = 26) -> list[str]:
    molecules = list(Chem.SDMolSupplier(str(source), removeHs=False))
    if len(molecules) != expected or any(molecule is None for molecule in molecules):
        raise ValueError(f"holdout must contain {expected} readable molecules")
    names = [molecule.GetProp("_Name").strip() for molecule in molecules]
    if (len(set(names)) != expected or
            len({directory_key(name) for name in names}) != expected or not all(names)):
        raise ValueError("holdout ligand names must be unique and nonempty")
    if any(molecule.GetNumConformers() != 1 or
           molecule.GetNumAtoms() == molecule.GetNumHeavyAtoms()
           for molecule in molecules):
        raise ValueError("holdout ligands require explicit-H bound conformers")
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        sanitized = list(Chem.SDMolSupplier(str(output), removeHs=False))
        if any(molecule is None or molecule.GetPropNames() for molecule in sanitized):
            raise ValueError("sanitized SDF contains properties")
        if ([molecule.GetProp("_Name").strip() for molecule in sanitized] != names or
                [Chem.GetFormalCharge(molecule) for molecule in sanitized] !=
                [Chem.GetFormalCharge(molecule) for molecule in molecules] or
                [Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=True)
                 for molecule in sanitized] !=
                [Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=True)
                 for molecule in molecules]):
            raise ValueError("existing sanitized SDF has different molecules")
        return names
    with Chem.SDWriter(str(output)) as writer:
        for molecule in molecules:
            molecule.ClearComputedProps()
            for prop in list(molecule.GetPropNames()):
                molecule.ClearProp(prop)
            writer.write(molecule)
    if any(molecule.GetPropNames() for molecule in Chem.SDMolSupplier(str(output), removeHs=False)):
        raise ValueError("holdout sanitization leaked SDF properties")
    return names


def score_one(name: str, protein, pocket, ligand, args, root: Path, signature):
    prepared, fragments, metadata = make_facio_like_fragments(
        pocket, ligand, include_waters=True, reference_protein=protein, use_hop=True,
        merge_disulfides=args.merge_disulfides,
        repair_peptide_numbering=args.repair_peptide_numbering)
    receptor_signature = (prepared, fragments[:-1], metadata["hop_pairs_gamess_indices"])
    if receptor_signature != signature:
        raise ValueError(f"shared receptor or caps changed: {name}")
    fragment = fragments[-1].fragment_id
    deck = render_fmo2_pieda(
        prepared, ligand, fragments, title=f"blind {args.target_name} {name}",
        dftb_parameter_dir="../params",
        fragment_bonds=[tuple(pair) for pair in signature[2]],
        deck_options=FMODeckOptions(use_pcm=args.use_pcm, dampex=args.dampex),
    )
    verify_fmo_deck_heavy_geometry(deck)
    frozen_scf = "MAXIT=200"
    if args.maxit != 200:
        if deck.count(frozen_scf) < 1:
            raise ValueError("cannot identify the frozen MAXIT control")
        control_line = " $CONTRL RUNTYP=ENERGY SCFTYP=RHF ISPHER=1 UNITS=ANGS MAXIT=200 $END"
        if deck.count(control_line) != 1:
            raise ValueError("cannot identify the $CONTRL MAXIT control")
        deck = deck.replace(control_line,
                            control_line.replace(frozen_scf, f"MAXIT={args.maxit}"), 1)
    if args.scf_solver == "diis":
        before = "DIIS=.FALSE. SOSCF=.TRUE."
        if deck.count(before) != 1:
            raise ValueError("cannot identify the frozen SCF solver line")
        deck = deck.replace(before, "DIIS=.TRUE. SOSCF=.FALSE.", 1)
    if args.nguess != 2:
        before = "NGUESS=2"
        if deck.count(before) != 1:
            raise ValueError("cannot identify the frozen FMO initial-guess control")
        deck = deck.replace(before, f"NGUESS={args.nguess}", 1)
    if args.mconv_dimer is not None and args.scf_solver != "fmo_soscf":
        before = "MCONV(2)=530 MCONV(4)=530"
        if deck.count(before) != 1:
            raise ValueError("cannot identify the frozen FMO dimer converger control")
        deck = deck.replace(before, f"MCONV(2)=530 MCONV(4)={args.mconv_dimer}", 1)
    if args.scc_mixer != "broyden":
        before = f"$DFTB SCC=.TRUE. DFTB3=.TRUE. DAMPXH=.TRUE. DAMPEX={args.dampex:.2f} $END"
        method = 2 if args.scc_mixer == "charge_diis" else -1
        after = ("$DFTB SCC=.TRUE. DFTB3=.TRUE. DAMPXH=.TRUE. "
                 f"DAMPEX={args.dampex:.2f} ITYPMX={method} $END")
        if deck.count(before) != 1 or len(after) > 80:
            raise ValueError("cannot identify the frozen SCC-DFTB charge mixer line")
        deck = deck.replace(before, after, 1)
    if args.scf_solver == "fmo_soscf":
        deck = enable_fmo_soscf(deck)
        deck = soscf_onset_deck(deck, args.soscf_onset)
        if args.mconv_dimer is not None:
            deck = deck.replace("MCONV(4)=576", f"MCONV(4)={args.mconv_dimer}", 1)
    if args.cnvdmp == 0:
        if deck.count("CNVDMP=50") != 1:
            raise ValueError("cannot identify the adaptive monomer precision control")
        deck = deck.replace("CNVDMP=50", "CNVDMP=0", 1)
    folder = root / "qm" / directory_key(name)
    folder.mkdir(parents=True, exist_ok=True)
    input_path = folder / "input.inp"
    run = folder / "run"
    if run.exists():
        previous = json.loads((run / "result.json").read_text())
        if previous.get("status") != "SUCCEEDED":
            old_hash = previous.get("input_sha256", "unknown")[:12]
            archive = folder / f"run.failed_{old_hash}"
            retry_index = 1
            while archive.exists():
                archive = folder / f"run.failed_{old_hash}_{retry_index}"
                retry_index += 1
            run.rename(archive)
            if input_path.exists():
                archived_input = folder / f"input.failed_{old_hash}.inp"
                retry_index = 1
                while archived_input.exists():
                    archived_input = folder / f"input.failed_{old_hash}_{retry_index}.inp"
                    retry_index += 1
                input_path.rename(archived_input)
    if input_path.exists() and input_path.read_text() != deck:
        raise ValueError(f"refusing changed successful QM deck for {name}")
    input_path.write_text(deck)
    parameter_link = folder / "params"
    if run.exists():
        record = json.loads((run / "result.json").read_text())
        if record.get("input_sha256") != digest(input_path) or record.get("status") != "SUCCEEDED":
            raise ValueError(f"existing calculation is stale or failed: {name}")
        if parameter_link.resolve(strict=True) != args.dftb_parameter_dir:
            raise ValueError(f"existing DFTB parameters changed: {name}")
    else:
        record = run_gamess(input_path, run, gamess_dir=args.gamess_dir,
                            library_dir=args.library_dir, ligand_fragment=fragment,
                            ncores=args.ncores,
                            timeout_seconds=args.timeout,
                            parameter_dir=args.dftb_parameter_dir)
    if record["status"] != "SUCCEEDED":
        raise RuntimeError(f"QM failure for {name}: {record['status']}")
    if args.scf_solver == "fmo_soscf":
        verify_fmo_soscf_log(run / "gamess.log")
    pairs = parse_pie(run / "gamess.log")
    selected = [pair for pair in pairs if fragment in (pair.fragment_i, pair.fragment_j)]
    if len(selected) != fragment - 1 or any(
            None in (pair.electrostatic_kcal_mol, pair.charge_transfer_kcal_mol,
                     pair.solvation_kcal_mol) for pair in selected):
        raise ValueError(f"incomplete PIEDA components for {name}")
    tpie = aggregate_tpie(pairs, fragment)
    official = find_selected_tie(run / "gamess.log", fragment)
    if official is not None and abs(official - tpie) > pie_roundoff_tolerance(len(selected)):
        raise ValueError(f"PIEDA total does not match GAMESS selected TIE for {name}")
    electrostatic = sum(pair.electrostatic_kcal_mol for pair in selected)
    transfer = sum(pair.charge_transfer_kcal_mol for pair in selected)
    solvation = sum(pair.solvation_kcal_mol for pair in selected)
    return {"ligand": name, "input_sha256": digest(input_path),
            "log_sha256": digest(run / "gamess.log"), "status": "SUCCEEDED",
            "n_receptor_fragments": fragment - 1, "tpie": tpie,
            "ES": electrostatic, "ES+SOLV": electrostatic + solvation,
            "ES+CT+SOLV": electrostatic + transfer + solvation}


def collect_cohort_scores(jobs, workers, expected, scoring, failure_gate=None):
    rows = []
    failures = []
    stopped = False

    def gate_triggered():
        return failure_gate is not None and failure_gate([name for name, error in failures])

    if workers == 1:
        for index, name, ligand in jobs:
            try:
                row = scoring(name, ligand)
            except Exception as error:
                failures.append((name, error))
                print(f"[{index}/{expected}] {name}: FAILED {error}", flush=True)
                if gate_triggered():
                    stopped = True
                    break
            else:
                rows.append(row)
                print(f"[{index}/{expected}] {name}: tPIE={row['tpie']:.3f}", flush=True)
        if failures:
            name, error = failures[0]
            if stopped:
                raise RuntimeError(f"QM failure for {name}: repeated fragment SCF failure; "
                                   "remaining cases deferred to whole-cohort rescue") from error
            raise RuntimeError(f"QM failure for {name}: {error}; "
                               f"{len(failures)} independent frozen failures; "
                               f"all {len(jobs)} requested cases attempted") from error
        return rows
    with ThreadPoolExecutor(max_workers=workers) as executor:
        remaining = iter(jobs)
        futures = {}

        def submit_next() -> bool:
            try:
                index, name, ligand = next(remaining)
            except StopIteration:
                return False
            future = executor.submit(scoring, name, ligand)
            futures[future] = (index, name)
            return True

        initial = min(workers, 2) if failure_gate is not None else workers
        probing = failure_gate is not None
        for _ in range(min(initial, len(jobs))):
            submit_next()
        while futures:
            complete, _ = wait(futures, return_when=FIRST_COMPLETED)
            for future in complete:
                index, name = futures.pop(future)
                try:
                    row = future.result()
                except Exception as error:
                    failures.append((name, error))
                    print(f"[{index}/{expected}] {name}: FAILED {error}", flush=True)
                else:
                    rows.append(row)
                    print(f"[{index}/{expected}] {name}: tPIE={row['tpie']:.3f}", flush=True)
            stopped = stopped or bool(gate_triggered())
            if probing and futures:
                continue
            probing = False
            if not stopped:
                for _ in range(workers - len(futures)):
                    if not submit_next():
                        break
        if failures:
            name, error = failures[0]
            if stopped:
                raise RuntimeError(f"QM failure for {name}: repeated fragment SCF failure; "
                                   "remaining cases deferred to whole-cohort rescue") from error
            raise RuntimeError(f"QM failure for {name}: {error}; "
                               f"{len(failures)} independent frozen failures; "
                               f"all {len(jobs)} requested cases attempted") from error
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("prepare", "minimize", "score"), required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--gamess-dir", type=Path)
    parser.add_argument("--library-dir", type=Path)
    parser.add_argument("--dftb-parameter-dir", type=Path)
    parser.add_argument("--experiment", default="E36 external target/source holdout")
    parser.add_argument("--set-name", default="galectin3_extra")
    parser.add_argument("--target-name", default="Galectin3")
    parser.add_argument("--expected", type=int, default=26)
    parser.add_argument("--cutoff", type=float, default=5.0)
    parser.add_argument("--merge-disulfides", action="store_true")
    parser.add_argument("--repair-peptide-numbering", action="store_true")
    parser.add_argument("--pocket-heavy-only", action="store_true",
                        help="select complete receptor residues using heavy-atom ligand distances")
    parser.add_argument("--exclude-pocket-residue", metavar="CHAIN:RESNAME:NUMBER",
                        help="registered QM-failing peripheral residue; retain the full prepared receptor")
    parser.add_argument("--scf-solver", choices=("soscf", "diis", "fmo_soscf"), default="soscf")
    parser.add_argument("--nguess", type=int, default=2)
    parser.add_argument("--mconv-dimer", type=int)
    parser.add_argument("--scc-mixer", choices=("broyden", "charge_diis", "gamess_soscf"),
                        default="broyden")
    parser.add_argument("--soscf-onset", type=float, choices=(0.25, 1.0), default=0.25,
                        help="bounded density-SOSCF startup threshold; does not relax final accuracy")
    parser.add_argument("--maxit", type=int, default=200)
    parser.add_argument("--cnvdmp", type=int, choices=(0, 50), default=50,
                        help="0 disables adaptive initial SCF precision; final thresholds remain unchanged")
    parser.add_argument("--dampex", type=float, default=4.0)
    parser.add_argument("--no-pcm", dest="use_pcm", action="store_false",
                        help="diagnostic gas-phase FMO arm; not equivalent to PCM production")
    parser.set_defaults(use_pcm=True)
    parser.add_argument("--only-ligand")
    parser.add_argument("--workers", type=int, default=1,
                        help="independent ligand jobs to run concurrently")
    parser.add_argument("--ncores", type=int, default=1,
                        help="GAMESS cores per ligand job")
    parser.add_argument("--restraint-k", type=float, default=10.0)
    parser.add_argument("--fix-ligand-heavy", action="store_true",
                        help="diagnostic: preserve exact source ligand heavy coordinates while relaxing hydrogens")
    parser.add_argument("--minimization-threads", type=int,
                        help="freeze OpenMM CPU thread count for reproducible rescue minimization")
    parser.add_argument("--repair-clashes", action="store_true",
                        help="allow restrained protein relaxation for diagnosed clash repair")
    parser.add_argument("--prepared-protein", type=Path)
    parser.add_argument("--reference-water-hydrogens", type=Path,
                        help="preserve validated water hydrogen orientations from treated source PDB")
    parser.add_argument("--ligand-sdf", type=Path,
                        help="optional label-free ligand SDF override after structure repair")
    parser.add_argument("--poses-root", type=Path)
    parser.add_argument("--reference-pocket-poses", type=Path)
    parser.add_argument("--fixed-pocket-membership", type=Path,
                        help="diagnostic: freeze receptor residue keys from a fingerprinted baseline pocket")
    parser.add_argument("--allow-reference-pocket-delta", action="store_true",
                        help="diagnostic only: use reference ligand poses for pocket selection despite receptor hash delta")
    parser.add_argument("--timeout", type=float, default=3600)
    args = parser.parse_args(argv)
    if args.scf_solver == "fmo_soscf" and (args.scc_mixer != "gamess_soscf" or
                                           args.mconv_dimer is not None):
        parser.error("fmo_soscf requires --scc-mixer gamess_soscf and no dimer override")
    if args.soscf_onset != 0.25 and args.scf_solver != "fmo_soscf":
        parser.error("nondefault SOSCF onset requires the explicit fmo_soscf solver")
    if not 1 <= args.maxit <= 200:
        parser.error("this GAMESS build requires MAXIT between 1 and 200")
    if (args.restraint_k <= 0 or args.timeout <= 0 or args.expected < 1 or args.cutoff <= 0
            or args.maxit < 1 or args.dampex <= 0 or args.workers < 1 or args.ncores < 1
            or args.nguess < 0 or args.minimization_threads is not None and
            args.minimization_threads < 1):
        parser.error("restraint-k, timeout, expected and cutoff must be positive")
    if args.fixed_pocket_membership is not None and (args.pocket_heavy_only or
                                                      args.reference_pocket_poses is not None or
                                                      args.exclude_pocket_residue is not None):
        parser.error("fixed pocket membership cannot mix other pocket selection overrides")
    if args.fix_ligand_heavy and (args.repair_clashes or args.poses_root is not None):
        parser.error("fixed ligand heavy atoms require local minimization and a fixed receptor")
    source_sdf = args.ligand_sdf or args.source / f"{args.set_name}_ligands.sdf"
    source_protein = args.source / f"{args.set_name}_protein.pdb"
    output = args.output
    clean_sdf = output / "ligands_no_properties.sdf"
    protein_pdb = output / f"{args.target_name.lower()}_prepared.pdb"
    names = sanitize(source_sdf, clean_sdf, expected=args.expected)
    cleaned_molecules = list(Chem.SDMolSupplier(str(clean_sdf), removeHs=False))
    ligand_charges = {molecule.GetProp("_Name").strip(): int(Chem.GetFormalCharge(molecule))
                      for molecule in cleaned_molecules}
    charge_counts = {str(charge): list(ligand_charges.values()).count(charge)
                     for charge in sorted(set(ligand_charges.values()))}
    manifest_path = output / "frozen_protocol.json"
    frozen = {"experiment": args.experiment, "n": args.expected,
              "source_sdf_sha256": digest(source_sdf), "source_pdb_sha256": digest(source_protein),
              "sanitized_sdf_sha256": digest(clean_sdf), "ligand_names": names,
              "protein_preparation": "PDB2PQR AMBER pH7.4; retain waters, exclude chloride",
              "minimization": (
                  f"free protein/waters OpenMM/OpenFF k={args.restraint_k:g} kJ/mol/A^2"
                  if args.repair_clashes else
                  f"fixed protein/waters OpenMM/OpenFF k={args.restraint_k:g} kJ/mol/A^2"),
              "qm": f"{args.cutoff:g}A union pocket, FMO2-DFTB3/PCM selected PIEDA; "
                    + ("per-case minimized receptor and ligand" if args.repair_clashes
                       else "fixed prepared receptor and minimized ligand"),
              "scores": ["tpie", "ES", "ES+SOLV", "ES+CT+SOLV"],
              "cohort_failure_policy": "two-case-preflight-and-repeated-monomer-gate-v1",
              "ligand_formal_charges": ligand_charges,
              "ligand_charge_counts": charge_counts,
              "ligand_state_policy": "preserve explicit-H molecular states and formal charges from the property-stripped source SDF; no automatic pKa or microstate changes",
              "mixed_charge_scores": "diagnostic only until solution-state protonation populations are modeled at the assay pH"}
    if args.workers != 1:
        frozen["execution_workers"] = args.workers
    if args.ncores != 1:
        frozen["gamess_ncores_per_job"] = args.ncores
    if args.cutoff != 5.0:
        frozen["cutoff_angstrom"] = args.cutoff
    if args.pocket_heavy_only:
        frozen["pocket_heavy_only"] = True
    if args.exclude_pocket_residue is not None:
        if args.pocket_heavy_only or not re.fullmatch(r"[^:]+:[A-Z]{3}:[0-9]+", args.exclude_pocket_residue):
            parser.error("selective residue exclusion needs a residue identity and the all-atom pocket")
        frozen["excluded_pocket_residue"] = args.exclude_pocket_residue
    if args.merge_disulfides:
        frozen["merge_disulfides"] = True
    if args.repair_peptide_numbering:
        frozen["repair_peptide_numbering"] = True
    if (args.target_name, args.set_name) != ("Galectin3", "galectin3_extra"):
        frozen.update(target_name=args.target_name, source_prefix=args.set_name)
    if args.scf_solver != "soscf":
        frozen["scf_solver"] = args.scf_solver
    if args.nguess != 2:
        frozen["nguess"] = args.nguess
    if args.mconv_dimer is not None:
        frozen["mconv_dimer"] = args.mconv_dimer
    if args.scc_mixer != "broyden":
        frozen["scc_mixer"] = args.scc_mixer
    if args.soscf_onset != 0.25:
        frozen["soscf_onset"] = args.soscf_onset
    if args.maxit != 200:
        frozen["maxit"] = args.maxit
    if args.cnvdmp != 50:
        frozen["cnvdmp"] = args.cnvdmp
    if args.dampex != 4.0:
        frozen["dampex"] = args.dampex
    if args.repair_clashes:
        frozen["repair_clashes"] = True
    if args.fix_ligand_heavy:
        frozen["fix_ligand_heavy"] = True
    if args.minimization_threads is not None:
        frozen["minimization_threads"] = args.minimization_threads
    if not args.use_pcm:
        frozen["pcm"] = False
    if args.poses_root is not None:
        frozen["poses_root"] = str(args.poses_root)
        frozen["poses_protocol_sha256"] = digest(args.poses_root / "frozen_protocol.json")
    if args.prepared_protein is not None:
        frozen["prepared_protein_override"] = str(args.prepared_protein)
        frozen["prepared_protein_override_sha256"] = digest(args.prepared_protein)
    if args.reference_water_hydrogens is not None:
        if args.prepared_protein is not None:
            parser.error("source water-hydrogen transfer requires a newly prepared receptor")
        frozen["reference_water_hydrogens"] = str(args.reference_water_hydrogens)
        frozen["reference_water_hydrogens_sha256"] = digest(args.reference_water_hydrogens)
    if args.ligand_sdf is not None:
        frozen["ligand_sdf_override"] = str(args.ligand_sdf)
        frozen["ligand_sdf_override_sha256"] = digest(args.ligand_sdf)
    if args.reference_pocket_poses is not None:
        frozen["reference_pocket_root"] = str(args.reference_pocket_poses)
        frozen["reference_pocket_protocol_sha256"] = digest(
            args.reference_pocket_poses / "frozen_protocol.json")
        if args.allow_reference_pocket_delta:
            frozen["reference_pocket_delta_diagnostic"] = True
    if args.fixed_pocket_membership is not None:
        frozen["fixed_pocket_membership"] = str(args.fixed_pocket_membership.resolve())
        frozen["fixed_pocket_membership_sha256"] = digest(args.fixed_pocket_membership)
    if manifest_path.exists() and not equivalent_frozen_protocol(
            json.loads(manifest_path.read_text()), frozen):
        raise ValueError("holdout structure or frozen protocol has changed")
    if not manifest_path.exists():
        manifest_path.write_text(json.dumps(frozen, indent=2) + "\n")
    if args.stage == "prepare":
        if args.prepared_protein is not None:
            copy_prepared_receptor(args.prepared_protein, protein_pdb)
        if not protein_pdb.exists():
            if args.prepared_protein is not None:
                shutil.copyfile(args.prepared_protein, protein_pdb)
            else:
                prepare_protein_pdb2pqr(source_protein, protein_pdb, ph=7.4,
                                        forcefield="AMBER", keep_water=True)
                if args.reference_water_hydrogens is not None:
                    preserve_source_waters(args.reference_water_hydrogens, protein_pdb)
        if args.reference_water_hydrogens is not None:
            audit_file = protein_pdb.with_suffix(".water_hydrogens.json")
            if not audit_file.is_file():
                raise ValueError("prepared receptor has no source-water transfer audit")
            audit = json.loads(audit_file.read_text())
            if (audit.get("source_sha256") != digest(args.reference_water_hydrogens) or
                    audit.get("prepared_sha256") != digest(protein_pdb)):
                raise ValueError("prepared source-water audit fingerprint changed")
        if args.prepared_protein is not None and digest(protein_pdb) != digest(args.prepared_protein):
            raise ValueError("prepared receptor does not match frozen reference")
        print(f"Prepared {len(names)} label-free molecules and {protein_pdb}", flush=True)
        return
    if not protein_pdb.exists():
        raise ValueError("first run --stage prepare")
    if args.stage == "minimize":
        from fmotpiefinal.core.minimization import minimize_complex

        if args.minimization_threads is not None:
            os.environ["OPENMM_CPU_THREADS"] = str(args.minimization_threads)
        for index, name in enumerate(names, start=1):
            folder = output / "minimized" / name
            record_path = folder / "minimization.json"
            if record_path.exists():
                record = json.loads(record_path.read_text())
                protein_input = record.get("protein_input")
                ligand_input = record.get("ligand_sdf")
                protein_matches = (protein_input is not None and
                                   digest(Path(protein_input)) == digest(protein_pdb))
                ligand_matches = (ligand_input is not None and
                                  digest(Path(ligand_input)) == digest(clean_sdf))
                protocol_matches = (
                    record.get("protein_coordinate_policy") == (
                        "free" if args.repair_clashes else
                        "fixed_prepared_receptor_and_waters") and
                    record.get("restraint_k_kj_mol_angstrom2") == args.restraint_k and
                    record.get("ligand_coordinate_policy") == (
                        "fixed_heavy_atoms_hydrogens_only" if args.fix_ligand_heavy else None) and
                    (folder / "minimized_complex.pdb").exists())
                if not (protein_matches and ligand_matches and protocol_matches):
                    record_path.unlink()
                    (folder / "minimized_complex.pdb").unlink(missing_ok=True)
                    minimize_complex(protein_pdb, clean_sdf, name, folder,
                                     restraint_k_kj_mol_a2=args.restraint_k,
                                     move_protein=args.repair_clashes,
                                     fix_ligand_heavy=args.fix_ligand_heavy)
            else:
                minimize_complex(protein_pdb, clean_sdf, name, folder,
                                 restraint_k_kj_mol_a2=args.restraint_k,
                                 move_protein=args.repair_clashes,
                                 fix_ligand_heavy=args.fix_ligand_heavy)
            print(f"[{index}/{args.expected}] {name} minimized", flush=True)
        return
    poses_root = args.poses_root or output
    if not all((poses_root / "minimized" / name / "minimization.json").exists()
               for name in names):
        raise ValueError(f"all {args.expected} blind minimized structures are required before scoring")
    if args.fix_ligand_heavy:
        from fmotpiefinal.pipeline.geometry_audit import audit as audit_geometry

        if args.poses_root is not None or args.repair_clashes:
            raise ValueError("hydrogen-only scoring requires locally minimized fixed receptor/ligand poses")
        geometry = audit_geometry(protein_pdb, clean_sdf, poses_root / "minimized")
        if (len(geometry["cases"]) != len(names) or
                {case["ligand"] for case in geometry["cases"]} != set(names) or
                any(case.get("minimized_pose") is None or
                    case["minimized_pose"]["n_severe_contacts"] != 0
                    for case in geometry["cases"])):
            raise ValueError("hydrogen-only candidate has unresolved pose contacts")
        for name in names:
            record = json.loads((poses_root / "minimized" / name / "minimization.json").read_text())
            if (record.get("ligand_coordinate_policy") != "fixed_heavy_atoms_hydrogens_only" or
                    record.get("max_ligand_heavy_displacement_angstrom", 1.0) > 0.0001):
                raise ValueError(f"hydrogen-only minimization changed ligand heavy atoms: {name}")
        report_path = output / "post_k10_geometry.json"
        content = json.dumps(geometry, indent=2) + "\n"
        if report_path.exists() and report_path.read_text() != content:
            raise ValueError("previous frozen hydrogen-only geometry report differs")
        report_path.write_text(content)
    if not all((args.gamess_dir, args.library_dir, args.dftb_parameter_dir)):
        parser.error("score requires --gamess-dir --library-dir --dftb-parameter-dir")
    args.dftb_parameter_dir = args.dftb_parameter_dir.resolve(strict=True)
    protein = read_pdb(protein_pdb)
    ligands = [load_bound_ligand(clean_sdf, name, poses_root / "minimized" / name /
                                 "minimized_complex.pdb") for name in names]
    component_region = validate_retained_component_region(protein_pdb, protein, ligands, args.cutoff)
    if component_region is not None:
        (output / "retained_component_qm_region.json").write_text(json.dumps(component_region, indent=2) + "\n")
    if args.repair_clashes:
        if len(names) != 1 or args.reference_pocket_poses is not None:
            raise ValueError("free-receptor QM requires a single case and no fixed reference pocket")
        protein = receptor_at_minimized_geometry(
            protein, poses_root / "minimized" / names[0] / "minimized_complex.pdb")
    if args.fixed_pocket_membership is not None:
        reference = json.loads(args.fixed_pocket_membership.read_text(encoding="utf-8"))
        pocket = select_fixed_pocket(protein, reference, args.target_name, names, args.cutoff)
    elif args.reference_pocket_poses is not None:
        reference = args.reference_pocket_poses
        reference_protein = reference / f"{args.target_name.lower()}_prepared.pdb"
        reference_matches = digest(protein_pdb) == digest(reference_protein)
        if not reference_matches and not args.allow_reference_pocket_delta:
            raise ValueError("reference pocket and scoring protein differ")
        reference_protocol = json.loads((reference / "frozen_protocol.json").read_text())
        reference_poses = Path(reference_protocol.get("poses_root", str(reference)))
        reference_ligands = [load_bound_ligand(reference / "ligands_no_properties.sdf", name,
                             reference_poses / "minimized" / name / "minimized_complex.pdb")
                             for name in names]
        pocket = shared_pocket(protein, reference_ligands, args.cutoff, include_waters=True,
                               heavy_only=args.pocket_heavy_only)
    else:
        pocket = shared_pocket(protein, ligands, args.cutoff, include_waters=True,
                               heavy_only=args.pocket_heavy_only)
    component_manifest = protein_pdb.with_suffix(".components.json")
    if component_manifest.exists():
        components = json.loads(component_manifest.read_text())
        if components.get("prepared_sha256") != digest(protein_pdb):
            raise ValueError("registered special-component receptor fingerprint changed")
        special_ids = {model["id"] for model in components["models"]}
        if special_ids & {residue_key(atom) for atom in pocket.atoms}:
            raise ValueError("retained MM component intersects actual QM pocket; validated special QM model required")
    if args.exclude_pocket_residue is not None:
        key = args.exclude_pocket_residue
        included = {residue_key(atom) for atom in pocket.atoms}
        if key not in included:
            raise ValueError(f"registered failing residue is not in the frozen shared pocket: {key}")
        chain, _, index = key.split(":")
        if any(f"{chain}:{name}:{int(index) + offset}" in included
               for name in {atom.residue_name for atom in pocket.atoms if atom.record == "ATOM"}
               for offset in (-1, 1)):
            raise ValueError(f"registered failing residue has contiguous selected protein neighbor: {key}")
        pocket = Protein(tuple(atom for atom in pocket.atoms if residue_key(atom) != key))
    prepared, fragments, metadata = make_facio_like_fragments(
        pocket, ligands[0], include_waters=True, reference_protein=protein, use_hop=True,
        merge_disulfides=args.merge_disulfides,
        repair_peptide_numbering=args.repair_peptide_numbering)
    signature = (prepared, fragments[:-1], metadata["hop_pairs_gamess_indices"])
    if args.reference_pocket_poses is not None:
        if not args.allow_reference_pocket_delta:
            original = (args.reference_pocket_poses / "qm" / directory_key(names[0]) / "input.inp").read_text()
            original_xyz = original.split(" $FMOXYZ\n", 1)[1].split(" $END", 1)[0].splitlines()
            reference_deck = render_fmo2_pieda(
                prepared, ligands[0], fragments, title=f"blind {args.target_name} {names[0]}",
                dftb_parameter_dir="../params", fragment_bonds=[tuple(pair) for pair in signature[2]],
                deck_options=FMODeckOptions(use_pcm=True, dampex=4.0))
            new_xyz = reference_deck.split(" $FMOXYZ\n", 1)[1].split(" $END", 1)[0].splitlines()
            if original_xyz[:len(prepared.atoms)] != new_xyz[:len(prepared.atoms)]:
                raise ValueError("reference receptor and cap coordinates changed")
    jobs = [(index, name, ligand) for index, (name, ligand) in
            enumerate(zip(names, ligands, strict=True), start=1)
            if args.only_ligand is None or name == args.only_ligand]
    if args.only_ligand is not None and len(jobs) != 1:
        raise ValueError(f"pilot ligand missing from frozen cohort: {args.only_ligand}")
    from fmotpiefinal.core.cohort_failure_gate import repeated_scf_failure

    def failure_gate(failed_names):
        evidence = repeated_scf_failure(output, [directory_key(name) for name in failed_names])
        if evidence is not None:
            (output / "cohort_failure_gate.json").write_text(json.dumps(evidence, indent=2) + "\n")
        return evidence

    rows = collect_cohort_scores(jobs, args.workers, args.expected,
                                 lambda name, ligand: score_one(name, protein, pocket, ligand,
                                                                args, output, signature),
                                 failure_gate=failure_gate)
    if args.only_ligand is not None:
        return
    rows.sort(key=lambda row: names.index(row["ligand"]))
    result = {"n_requested": args.expected, "n_succeeded": len(rows),
              "frozen_protocol_sha256": digest(manifest_path),
              "protein_sha256": digest(protein_pdb), "pocket_atoms": len(pocket.atoms),
              "predictions": rows}
    (output / "blind_predictions.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
