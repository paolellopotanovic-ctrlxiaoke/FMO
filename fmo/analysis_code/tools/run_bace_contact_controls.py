"""Full-cohort mechanistic BACE controls; never choose a branch using labels."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from copy import deepcopy
import json
import multiprocessing
from pathlib import Path

import numpy as np

from fmotpiefinal.core.bound_geometry import load_bound_ligand
from fmotpiefinal.core.structure import read_pdb
from run_cdk8_solution_cycle import digest, write_json
from run_local_sidechain_transfer import directory_key, evaluate, score_job


def execute(payload):
    cohort, name, output, arm = payload
    from fmotpiefinal.core.minimization import minimize_complex
    parent = Path(output) / arm / "BACE" / directory_key(name)
    preparation_path = parent / "preparation_status.json"
    try:
        if not preparation_path.exists():
            arguments = {"restraint_k_kj_mol_a2": 10.0, "move_protein": False}
            if arm == "contact_carboxylates":
                arguments.update(local_flexible_residues=tuple(map(tuple, cohort["local_residues"])),
                                 receptor_restraint_k_kj_mol_a2=2.5)
            elif arm == "water_hydrogens_only":
                arguments.update(relax_water_hydrogens=True, fix_ligand_all=True,
                                 starting_ligand_complex=cohort["complexes"][name])
            else:
                raise ValueError("unknown frozen control")
            result = minimize_complex(cohort["protein"], cohort["sdf"], name,
                                      parent / "minimization", **arguments)
            write_json(preparation_path, {"status": "SUCCEEDED", "complex": str(result.complex_pdb)})
        return arm, name, score_job(cohort, name, Path(output) / arm, 2)[2]
    except Exception as error:
        record = {"status": "FAILED", "error": f"{type(error).__name__}: {error}"}
        write_json(parent / "summary.json", record)
        return arm, name, record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("positive workers required")
    source, output, workspace = args.source.resolve(), args.output.resolve(), args.workspace.resolve()
    original = json.loads((source / "frozen_protocol.json").read_text())["cohorts"]["BACE"]
    for path, expected in original["input_sha256"].items():
        if digest(path) != expected:
            raise ValueError("verified baseline evidence changed")
    protein = read_pdb(original["protein"])
    ligand_points = np.asarray([atom.xyz for name in original["names"] for atom in
        load_bound_ligand(original["sdf"], name, original["complexes"][name]).atoms if atom.element != "H"])
    selected = sorted({(atom.chain, atom.residue_number, atom.residue_name) for atom in protein.atoms
        if atom.residue_name in {"ASP", "ASH", "GLU", "GLH"} and atom.atom_name in {"OD1", "OD2", "OE1", "OE2"}
        and np.linalg.norm(ligand_points - np.asarray(atom.xyz), axis=1).min() <= 3.5})
    if not selected:
        raise ValueError("no contacting carboxylate mechanism identified")
    cohorts = {arm: deepcopy(original) for arm in ("contact_carboxylates", "water_hydrogens_only")}
    cohorts["contact_carboxylates"]["local_residues"] = selected
    cohorts["water_hydrogens_only"]["local_residues"] = []
    files = [Path(__file__), Path(__file__).with_name("run_local_sidechain_transfer.py"),
             Path(__file__).parents[1] / "src/fmotpiefinal/core/minimization.py"]
    frozen = {"schema": "bace-contact-mechanism-controls-v1", "cohorts": cohorts,
        "source_protocol_sha256": digest(source / "frozen_protocol.json"),
        "code_sha256": {str(path.resolve()): digest(path) for path in files},
        "labels_used_for_generation": False, "production_default_changed": False,
        "scope": "Matched prepared-geometry controls, not a fresh raw preparation or free-energy calculation",
        "selection": "whole-cohort carboxylate oxygen within 3.5A of ligand heavy atoms; fixed receptor backbone and water; k10 ligand / k2.5 selected sidechains",
        "water_control": "baseline k10 ligand fully fixed, protein and water oxygens fixed, only existing water H can move; no water deletion, proton transfer or state change",
        "primary_readout": "complete 36-case tPIE; no best-arm automatic adoption"}
    path = output / "frozen_protocol.json"
    if path.exists() and json.loads(path.read_text()) != frozen:
        raise ValueError("refusing changed registered experiment")
    write_json(path, frozen)
    results = {arm: {"BACE": {}} for arm in cohorts}
    baseline = workspace / "FMOtpiefinal_HIV1lab/runs/all_accelerated_robust_20261020"
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context("spawn")) as pool:
        futures = [pool.submit(execute, (cohort, name, str(output), arm))
                   for arm, cohort in cohorts.items() for name in cohort["names"]]
        for future in as_completed(futures):
            arm, name, record = future.result()
            results[arm]["BACE"][name] = record
            write_json(output / "progress.json", results)
            print(arm, name, record["status"], record.get("tpie", record.get("error")), flush=True)
    for arm in cohorts:
        evaluation = evaluate(workspace, baseline, results[arm], output / arm)
        succeeded = evaluation["groups"]["BACE"]["n_succeeded"]
        write_json(output / arm / "final_status.json", {"status": "COMPLETE" if succeeded == 36 else "INCOMPLETE",
            "n_requested": 36, "n_succeeded": succeeded, "evaluation_sha256": digest(output / arm / "evaluation.json")})


if __name__ == "__main__":
    main()
