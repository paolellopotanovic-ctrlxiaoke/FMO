"""Boundary-free ligand/water energy-model diagnostic; not a binding free energy."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import time
import traceback

import numpy as np

from run_matched_openfe import digest, write_json


HARTREE_TO_KCAL_MOL = 627.5094740631


def closest_complete_water(protein, ligand):
    from fmotpiefinal.core.structure import residue_key

    waters = {}
    for atom in protein.atoms:
        if atom.residue_name in {"HOH", "WAT", "TIP"}:
            waters.setdefault(residue_key(atom), []).append(atom)
    heavy = np.asarray([atom.xyz for atom in ligand.atoms if atom.element.upper() != "H"])
    candidates = []
    for key, atoms in waters.items():
        if sorted(atom.element.upper() for atom in atoms) != ["H", "H", "O"]:
            raise ValueError(f"incomplete water must be repaired before reference: {key}")
        oxygen = next(atom for atom in atoms if atom.element.upper() == "O")
        hydrogens = [atom for atom in atoms if atom.element.upper() == "H"]
        vectors = np.asarray([np.asarray(atom.xyz) - oxygen.xyz for atom in hydrogens])
        distances = np.linalg.norm(vectors, axis=1)
        angle = float(np.degrees(np.arccos(np.clip(np.dot(*vectors) / np.prod(distances), -1, 1))))
        if not np.all((distances >= .7) & (distances <= 1.25)) or not 80 <= angle <= 130:
            raise ValueError(f"invalid water geometry requires preparation diagnosis: {key}")
        nearest = float(np.linalg.norm(heavy - oxygen.xyz, axis=1).min())
        candidates.append((nearest, key, atoms))
    if not candidates:
        raise ValueError("no retained complete structural water; no artificial water added")
    distance, key, atoms = min(candidates, key=lambda item: (item[0], item[1]))
    if distance < 1.5:
        raise ValueError("severe ligand/water overlap")
    return key, distance, atoms


def reference_energy(atoms, charge, basis, ghost_indices, directory, threads):
    from pyscf import gto, lib, mp, scf
    import pyscf

    directory.mkdir(parents=True, exist_ok=False)
    lib.num_threads(threads)
    symbols = [("ghost-" if index in ghost_indices else "") + atom.element.capitalize()
               for index, atom in enumerate(atoms)]
    molecule = gto.M(atom=list(zip(symbols, [atom.xyz for atom in atoms], strict=True)),
                     basis=basis, unit="Angstrom", charge=charge, spin=0,
                     output=str(directory / "pyscf.log"), verbose=4, max_memory=6000)
    mean_field = scf.RHF(molecule).density_fit()
    mean_field.conv_tol = 1e-10
    mean_field.max_cycle = 200
    mean_field.chkfile = str(directory / "scf.chk")
    hartree_fock = float(mean_field.kernel())
    if not mean_field.converged:
        raise ValueError("reference RHF did not converge; no relaxed tolerance substitution")
    perturbation = mp.MP2(mean_field)
    perturbation.kernel(with_t2=False)
    total = float(perturbation.e_tot)
    if not np.isfinite([hartree_fock, total]).all():
        raise ValueError("nonfinite correlated reference energy")
    record = {"pyscf": pyscf.__version__, "basis": basis, "method": "density-fitted all-electron RHF/MP2",
        "hf_hartree": hartree_fock, "mp2_hartree": total, "charge": charge,
        "ghost_indices_zero_based": sorted(ghost_indices), "n_basis": molecule.nao_nr(),
        "scf_converged": True, "scf_tolerance": mean_field.conv_tol,
        "scf_stability_test_performed": False, "threads": threads}
    write_json(directory / "result.json", record)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--target", choices=("BACE", "CDK8"), required=True)
    parser.add_argument("--case-index", type=int, default=0)
    parser.add_argument("--basis", choices=("def2-svpd", "def2-tzvpd"), default="def2-svpd")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    if args.case_index < 0 or args.threads < 1:
        parser.error("nonnegative case index and positive threads required")
    workspace = args.workspace.resolve()
    sys.path.insert(0, str(workspace / "src"))
    from fmotpiefinal.core.bound_geometry import load_bound_ligand
    from fmotpiefinal.core.fmo_input import FMODeckOptions, render_fmo2_pieda
    from fmotpiefinal.core.fragmentation import Fragment
    from fmotpiefinal.core.structure import Protein, read_pdb
    from sophosqm_baseline.relative_total_energy import execute_deck

    plan = json.loads(args.protocol.read_text())
    cohort = plan["targets"][args.target]
    for path, expected in cohort["input_sha256"].items():
        if digest(path) != expected:
            raise ValueError("frozen geometry input changed")
    name = cohort["names"][args.case_index]
    ligand = load_bound_ligand(cohort["sdf"], name, cohort["poses"][name])
    water_key, nearest, water_atoms = closest_complete_water(read_pdb(cohort["protein"]), ligand)
    atoms = list(ligand.atoms) + water_atoms
    directory = args.output.resolve() / args.target / f"case_{args.case_index:03d}" / args.basis
    directory.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    inventory = {"schema": "matched-ligand-water-energy-reference-v1", "target": args.target, "ligand": name,
        "case_selection": "registered cohort index, independent of experimental values", "water": water_key,
        "nearest_ligand_heavy_to_water_O_angstrom": nearest, "basis": args.basis,
        "ligand_charge": ligand.formal_charge, "water_charge": 0,
        "atoms": [{"element": atom.element, "xyz_angstrom": atom.xyz} for atom in atoms],
        "n_ligand_atoms": len(ligand.atoms), "protocol_sha256": digest(args.protocol),
        "driver_sha256": digest(__file__), "thermodynamic_free_energy": False,
        "boundary_policy": "entire ligand plus intact water; no truncated receptor covalent bonds",
        "interpretation": "Gas-phase local model comparison only; not PCM tPIE replacement or final binding reference. Basis, SCF stability and omitted environment remain limitations."}
    write_json(directory / "input_manifest.json", inventory)
    write_json(directory / "status.json", {"status": "RUNNING", "pid": os.getpid()})
    try:
        tasks = {"complex": (atoms, ligand.formal_charge, set()),
            "ligand_cp": (atoms, ligand.formal_charge, set(range(len(ligand.atoms), len(atoms)))),
            "water_cp": (atoms, 0, set(range(len(ligand.atoms)))),
            "ligand": (list(ligand.atoms), ligand.formal_charge, set()),
            "water": (water_atoms, 0, set())}
        references = {state: reference_energy(*values[:2], args.basis, values[2], directory / state, args.threads)
                      for state, values in tasks.items()}
        energies = {}
        for method in ("hf_hartree", "mp2_hartree"):
            energies[method.replace("_hartree", "")] = {
                "interaction_uncorrected_kcal_mol": (references["complex"][method] - references["ligand"][method] - references["water"][method]) * HARTREE_TO_KCAL_MOL,
                "interaction_counterpoise_kcal_mol": (references["complex"][method] - references["ligand_cp"][method] - references["water_cp"][method]) * HARTREE_TO_KCAL_MOL}
        registration = workspace / "FMOtpiefinal_HIV1lab/runs/graph_sidechain_transfer_20261021/frozen_protocol.json"
        original = json.loads(registration.read_text())["cohorts"][args.target]
        runtime = original["runtime"][name]
        (directory / "params").symlink_to(Path(original["parameters"]).resolve(), target_is_directory=True)
        receptor = Protein(tuple(water_atoms))
        water_fragment = Fragment(1, "WATER", "water", tuple(atom.index for atom in water_atoms))
        ligand_fragment = Fragment(2, "LIGAND", "ligand", tuple(atom.index for atom in ligand.atoms), charge=ligand.formal_charge)
        ligand_only = Fragment(1, "LIGAND", "ligand", ligand_fragment.atom_indices, charge=ligand.formal_charge)
        dftb_tasks = {"complex": (receptor, ligand, [water_fragment, ligand_fragment]),
            "ligand": (Protein(()), ligand, [ligand_only]), "water": (receptor, None, [water_fragment])}
        dftb = {}
        for state, values in dftb_tasks.items():
            deck = render_fmo2_pieda(*values, dftb_parameter_dir="../params",
                full_energy=True, deck_options=FMODeckOptions(use_pcm=False))
            deck += " $DFT DC=.TRUE. IDCVER=4 $END\n"
            result = execute_deck(deck, directory / ("dftb_" + state), gamess_dir=runtime["gamess_dir"],
                library_dir=runtime["library_dir"], timeout=1200, ncores=2)
            if result["status"] != "SUCCEEDED":
                raise ValueError(f"matched DFTB reference failed: {state}")
            dftb[state] = result["energy_hartree"]
        energies["dftb3_dispersion"] = {"interaction_kcal_mol": (dftb["complex"] - dftb["ligand"] - dftb["water"]) * HARTREE_TO_KCAL_MOL}
        write_json(directory / "status.json", {"status": "SUCCEEDED", "energies": energies,
            "wall_seconds": time.monotonic() - started, "affinity_metric_claim_allowed": False})
    except Exception as error:
        write_json(directory / "status.json", {"status": "FAILED", "error": f"{type(error).__name__}: {error}",
            "traceback": traceback.format_exc(), "wall_seconds": time.monotonic() - started})
        raise


if __name__ == "__main__":
    main()
