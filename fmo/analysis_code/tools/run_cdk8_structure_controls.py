"""Frozen full-cohort CDK8 receptor-flexibility and water-H controls."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import io
import json
from pathlib import Path
import re
import sys
import time

import numpy as np

from run_cdk8_solution_cycle import array, digest, group, metrics, write_json
from run_cdk8_bound_qm import replace_coordinates


NAMES = ("12", "13", "14", "15", "16", "17", "19", "20", "21", "22")
ARMS = ("baseline_reproduction", "local_sidechains", "water_hydrogens")


def topology_signature(deck):
    content = group(deck, "FMO")[1]
    bonds = re.search(r"\$FMOBND\s+(.*?)\$END", deck, re.S | re.I)
    return {"charges": array(content, "ICHARG", "FRGNAM"),
            "names": array(content, "FRGNAM", "INDAT"),
            "membership": array(content, "INDAT", r"(?!)"),
            "bonds": bonds[1].split() if bonds else []}


def coordinate_array(deck):
    rows = [line.split() for line in group(deck, "FMOXYZ")[1].strip().splitlines()]
    return [row[0].upper() for row in rows], np.asarray([row[2:] for row in rows], float)


def minimize_job(payload):
    from fmotpiefinal.core.minimization import minimize_complex
    started = time.monotonic()
    arm, name, prepared, output, selected = payload
    prepared, output = Path(prepared), Path(output)
    destination = output / arm / name / "minimization"
    arguments = {"restraint_k_kj_mol_a2": 10.0, "move_protein": False}
    if arm == "local_sidechains":
        arguments.update(local_flexible_residues=tuple(map(tuple, selected)),
                         receptor_restraint_k_kj_mol_a2=2.5)
    elif arm == "water_hydrogens":
        arguments.update(relax_water_hydrogens=True, fix_ligand_all=True,
            starting_ligand_complex=prepared / f"minimized/{name}/minimized_complex.pdb")
    try:
        result = minimize_complex(prepared / "cdk8_prepared.pdb",
            prepared / "ligands_no_properties.sdf", name, destination, **arguments)
        status = {"status": "SUCCEEDED", "complex": str(result.complex_pdb),
                  "wall_seconds": time.monotonic() - started}
    except Exception as error:
        status = {"status": "FAILED", "error": f"{type(error).__name__}: {error}",
                  "wall_seconds": time.monotonic() - started}
    write_json(output / arm / name / "preparation_status.json", status)
    return arm, name, status


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--prepare-workers", type=int, default=4)
    parser.add_argument("--qm-workers", type=int, default=5)
    args = parser.parse_args()
    if min(args.prepare_workers, args.qm_workers) < 1:
        parser.error("worker counts must be positive")
    workspace, reference, prepared, output = (value.resolve() for value in
        (args.workspace, args.reference, args.prepared, args.output))
    sys.path.insert(0, str(workspace / "src"))
    from fmotpiefinal.core.structure import Protein, read_pdb, residue_key
    from fmotpiefinal.core.bound_geometry import load_bound_ligand, shared_pocket
    from fmotpiefinal.core.fragmentation import make_facio_like_fragments
    from fmotpiefinal.core.fmo_input import render_fmo2_pieda
    from fmotpiefinal.core.minimization import _water_compatible_pdb
    from fmotpiefinal.core.fmo_convergence import final_fmo_convergence
    from fmotpiefinal.core.solution_cycle import pcm_energy_inventory
    from fmotpiefinal.pipeline.target_stages import receptor_at_minimized_geometry
    from sophosqm_baseline.relative_total_energy import execute_deck
    from sophosqm_baseline.parser import parse_pie
    protein = read_pdb(prepared / "cdk8_prepared.pdb")
    ligands = {name: load_bound_ligand(prepared / "ligands_no_properties.sdf", name,
        prepared / f"minimized/{name}/minimized_complex.pdb") for name in NAMES}
    pocket = shared_pocket(protein, list(ligands.values()), 5.0)
    keys = {residue_key(atom) for atom in pocket.atoms}
    ligand_heavy = np.asarray([atom.xyz for ligand in ligands.values()
                              for atom in ligand.atoms if atom.element.upper() != "H"])
    selected = sorted({(atom.chain, atom.residue_number, atom.residue_name)
        for atom in protein.atoms if atom.record == "ATOM" and atom.element.upper() != "H"
        and float(np.min(np.linalg.norm(ligand_heavy - atom.xyz, axis=1))) <= 4.0})
    parameters = (reference / f"pcm/{NAMES[0]}/params").resolve(strict=True)

    def render(receptor, ligand):
        fixed_pocket = Protein(tuple(atom for atom in receptor.atoms if residue_key(atom) in keys))
        capped, fragments, metadata = make_facio_like_fragments(fixed_pocket, ligand,
            reference_protein=receptor, merge_disulfides=True, repair_peptide_numbering=True)
        return render_fmo2_pieda(capped, ligand, fragments, dftb_parameter_dir="../params",
            fragment_bonds=metadata["hop_pairs_gamess_indices"], full_energy=True)

    for name in NAMES:
        regenerated = render(protein, ligands[name])
        source = (reference / f"pcm/{name}/complex/input.inp").read_text()
        if topology_signature(regenerated) != topology_signature(source):
            raise ValueError(f"{name}: regenerated baseline topology differs")
        elements, positions = coordinate_array(regenerated)
        original_elements, original_positions = coordinate_array(source)
        if elements != original_elements or not np.array_equal(positions, original_positions):
            raise ValueError(f"{name}: regenerated baseline caps/coordinates differ")
    plan = {"names": list(NAMES), "arms": list(ARMS), "labels_used_for_generation": False,
        "ligand_restraint_k_kj_mol_angstrom2": 10.0,
        "local_receptor_restraint_k_kj_mol_angstrom2": 2.5,
        "local_selection": "4A heavy-atom union over all ten baseline k10 poses",
        "local_residues": selected, "pocket_residues": sorted(keys), "pocket_cutoff_angstrom": 5.0,
        "receptor_sha256": digest(prepared / "cdk8_prepared.pdb"),
        "sdf_sha256": digest(prepared / "ligands_no_properties.sdf"),
        "reference_inputs": {name: digest(reference / f"pcm/{name}/complex/input.inp") for name in NAMES},
        "starting_complexes": {name: digest(prepared / f"minimized/{name}/minimized_complex.pdb") for name in NAMES},
        "driver_sha256": digest(__file__),
        "minimization_sha256": digest(Path(__file__).parents[1] / "src/fmotpiefinal/core/minimization.py"),
        "baseline_deck_regeneration": "all ten coordinates, elements, fragments, charges and HOPs exact",
        "water_control": "all water H mobile; ligand, protein and water O fixed to baseline k10",
        "local_control": "selected sidechains mobile; backbone, waters, remote atoms fixed",
        "parameters": str(parameters), "qm_ncores_per_job": 2,
        "production_default_changed": False}
    output.mkdir(parents=True, exist_ok=True)
    existing = output / "frozen_protocol.json"
    if existing.exists() and json.loads(existing.read_text()) != json.loads(json.dumps(plan)):
        raise ValueError("frozen structure control plan changed")
    write_json(existing, plan)
    _water_compatible_pdb(prepared / "cdk8_prepared.pdb")
    started = time.monotonic()
    statuses = {}
    pending = []
    for arm in ARMS:
        for name in NAMES:
            cached = output / arm / name / "preparation_status.json"
            if cached.exists() and json.loads(cached.read_text()).get("status") == "SUCCEEDED":
                statuses[f"{arm}/{name}"] = json.loads(cached.read_text())
            else:
                pending.append((arm, name, str(prepared), str(output), selected))
    with ProcessPoolExecutor(max_workers=args.prepare_workers) as executor:
        futures = [executor.submit(minimize_job, payload) for payload in pending]
        for future in as_completed(futures):
            arm, name, status = future.result()
            statuses[f"{arm}/{name}"] = status
            write_json(output / "preparation_progress.json", statuses)
            print("PREP", arm, name, status["status"], status.get("error", ""), flush=True)

    def score(arm, name):
        parent = output / arm / name
        try:
            if statuses[f"{arm}/{name}"]["status"] != "SUCCEEDED":
                raise ValueError("preparation failed")
            complex_path = parent / "minimization/minimized_complex.pdb"
            receptor = receptor_at_minimized_geometry(protein, complex_path)
            ligand = load_bound_ligand(prepared / "ligands_no_properties.sdf", name, complex_path)
            regenerated = render(receptor, ligand)
            source = (reference / f"pcm/{name}/complex/input.inp").read_text()
            if topology_signature(regenerated) != topology_signature(source):
                raise ValueError("prepared topology differs from frozen reference")
            rows = group(regenerated, "FMOXYZ")[1].strip().splitlines()
            deck = replace_coordinates(source, rows).replace("IPIEDA=0", "IPIEDA=1")
            (parent / "params").symlink_to(parameters, target_is_directory=True) if not (parent / "params").exists() else None
            directory = parent / "qm"
            directory.mkdir(exist_ok=True)
            (directory / "input.inp").write_text(deck)
            runtime = json.loads((reference / f"pcm/{name}/complex/result.json").read_text())
            result = execute_deck(deck, directory, gamess_dir=runtime["gamess_dir"],
                library_dir=runtime["library_dir"], ncores=2, timeout=1800)
            if result["status"] != "SUCCEEDED":
                raise ValueError(f"QM failed: {result}")
            log = (directory / "gamess.log").read_text()
            convergence = final_fmo_convergence(deck, io.StringIO(log))
            inventory = pcm_energy_inventory(log)
            count = int(re.search(r"NFRAG=(\d+)", deck)[1])
            pairs = [pair for pair in parse_pie(directory / "gamess.log")
                     if count in (pair.fragment_i, pair.fragment_j)]
            if len(pairs) != count - 1:
                raise ValueError("incomplete ligand pair inventory")
            result = {"status": "SUCCEEDED", "tpie": sum(pair.value_kcal_mol for pair in pairs),
                "pair_count": len(pairs), "convergence": convergence, "pcm_inventory": inventory,
                "input_sha256": digest(directory / "input.inp"),
                "receptor_max_displacement_angstrom": float(np.max(np.linalg.norm(
                    np.asarray([atom.xyz for atom in receptor.atoms]) -
                    np.asarray([atom.xyz for atom in protein.atoms]), axis=1)))}
        except Exception as error:
            result = {"status": "FAILED", "error": f"{type(error).__name__}: {error}"}
        write_json(parent / "summary.json", result)
        return arm, name, result

    results = {}
    with ThreadPoolExecutor(max_workers=args.qm_workers) as executor:
        futures = [executor.submit(score, arm, name) for arm in ARMS for name in NAMES]
        for future in as_completed(futures):
            arm, name, result = future.result()
            results.setdefault(arm, {})[name] = result
            write_json(output / "qm_progress.json", results)
            print("QM", arm, name, result["status"], result.get("tpie", result.get("error")), flush=True)
    from rdkit import Chem
    sdf = workspace / "vendor/public_binding_free_energy_benchmark/fep_benchmark_inputs/structure_inputs/misc/cdk8_koehler_ligands.sdf"
    molecules = {molecule.GetProp("_Name").strip(): molecule
                 for molecule in Chem.SDMolSupplier(str(sdf), removeHs=False) if molecule is not None}
    truth = [float(molecules[name].GetProp("r_exp_dg")) for name in NAMES]
    evaluation = {}
    for arm in ARMS:
        successful = [name for name in NAMES if results[arm][name]["status"] == "SUCCEEDED"]
        evaluation[arm] = {"n_succeeded": len(successful), "n_expected": len(NAMES)}
        if len(successful) == len(NAMES):
            evaluation[arm]["metrics"] = metrics(truth, [results[arm][name]["tpie"] for name in NAMES])
    baseline_verified = False
    historical_path = prepared.parents[2] / "independent_verification.json"
    if historical_path.exists() and evaluation[ARMS[0]]["n_succeeded"] == len(NAMES):
        historical = {record["ligand"]: record["tpie"] for record in json.loads(historical_path.read_text())["records"]}
        difference = max(abs(results[ARMS[0]][name]["tpie"] - historical[name]) for name in NAMES)
        evaluation[ARMS[0]]["max_historical_score_difference_kcal_mol"] = difference
        baseline_verified = difference <= 0.05
    write_json(output / "evaluation.json", {"arms": evaluation, "baseline_matched_within_0_05_kcal": baseline_verified,
        "wall_seconds": time.monotonic() - started, "results": results,
        "target_spearman": 0.6, "experimental_labels_used_only_after_generation": True})
    print(json.dumps(evaluation, indent=2), flush=True)


if __name__ == "__main__":
    main()
