"""Graph-triggered guanidino NH2 reconstruction with full-cohort QM controls."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import multiprocessing
from pathlib import Path
import time

import numpy as np
from rdkit import Chem

from fmotpiefinal.core.conjugated_amine_geometry import induced_conjugated_nh2, reconstruct_flagged_nh2
from fmotpiefinal.core.solution_cycle import HARTREE_TO_KCAL_MOL
from report_local_sidechain_transfer import load_truth
from run_cdk8_bound_qm import replace_coordinates
from run_cdk8_solution_cycle import digest, metrics, write_json
from run_direction_preserving_hydrogen_control import worker as execute_geometry_pair
from run_hydrogen_normalized_strain import fixed_heavy_displacement, transplant_ligand, xyz_rows
from run_local_sidechain_transfer import directory_key


def scan_registry(registry):
    cohorts, molecules = {}, {}
    for target, cohort in registry["cohorts"].items():
        matches = {molecule.GetProp("_Name").strip(): molecule for molecule in Chem.SDMolSupplier(cohort["sdf"], removeHs=False) if molecule is not None}
        if set(matches) != set(cohort["names"]):
            raise ValueError("complete source chemical inventory required")
        cases = {}
        for name in cohort["names"]:
            molecule = matches[name]
            bound_rows = xyz_rows(Path(cohort["sources"][name]).read_text())[-molecule.GetNumAtoms():]
            if [row.split()[0].upper() for row in bound_rows] != [atom.GetSymbol().upper() for atom in molecule.GetAtoms()]:
                raise ValueError("source ligand atom mapping mismatch")
            bound = np.asarray([row.split()[2:] for row in bound_rows], float)
            diagnostics = induced_conjugated_nh2(molecule, molecule.GetConformer().GetPositions(), bound)
            repaired = reconstruct_flagged_nh2(molecule, bound, diagnostics)
            if not any(row["new_deformation_flag"] for row in diagnostics) and not np.array_equal(bound, repaired):
                raise ValueError("untriggered geometry was changed")
            cases[name] = {"groups": diagnostics, "triggered": any(row["new_deformation_flag"] for row in diagnostics)}
        cohorts[target] = {"n": len(cases), "n_triggered": sum(row["triggered"] for row in cases.values()),
                           "n_existing_input_flags": sum(any(group["already_flagged_input"] for group in row["groups"]) for row in cases.values()), "cases": cases}
        molecules[target] = matches
    return {"schema": "induced-conjugated-nh2-all-cohort-audit-v1", "cohorts": cohorts,
            "labels_used": False, "departure_threshold_degrees": 30., "increase_threshold_degrees": 20.,
            "n": sum(row["n"] for row in cohorts.values()), "scope": "Chemical group and pre/post geometry screen; not proof of a force-field bug or affinity prediction accuracy"}, molecules


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--prior", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=16)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("positive workers required")
    start = time.monotonic()
    registry = json.loads(args.registry.read_text())
    source = json.loads(args.protocol.read_text())
    prior = json.loads(args.prior.read_text())
    target, names = source["target"], source["names"]
    if prior["status"] != "COMPLETE_VERIFIED" or set(prior["scores"]) != set(names):
        raise ValueError("full matched cycle reference required")
    audit, molecules = scan_registry(registry)
    write_json(args.output/"graph_audit.json", audit)
    cohort = registry["cohorts"][target]
    for path, fingerprint in cohort["input_sha256"].items():
        if digest(path) != fingerprint:
            raise ValueError("source registry evidence changed")
    tasks = {}
    for name in names:
        ligand_key = source["case_tasks"][name]["ligand"]
        complex_key = source["case_tasks"][name]["complex"]
        ligand_path, complex_path = [Path(source["tasks"][key]["input"]) for key in (ligand_key, complex_key)]
        for key, path in ((ligand_key, ligand_path), (complex_key, complex_path)):
            if digest(path) != prior["rows"][key]["input_sha256"]:
                raise ValueError("source cycle input changed")
        ligand_deck = ligand_path.read_text()
        before = np.asarray([row.split()[2:] for row in xyz_rows(ligand_deck)], float)
        molecule = molecules[target][name]
        diagnostics = induced_conjugated_nh2(molecule, molecule.GetConformer().GetPositions(), before)
        points = reconstruct_flagged_nh2(molecule, before, diagnostics)
        rows = [f"{row.split()[0]} {row.split()[1]} {point[0]:.10f} {point[1]:.10f} {point[2]:.10f}"
                for row, point in zip(xyz_rows(ligand_deck), points, strict=True)]
        fixed_heavy_displacement(xyz_rows(ligand_deck), rows)
        parent = args.output.resolve()/target/directory_key(name)
        parent.mkdir(parents=True, exist_ok=True)
        link = parent/"params"
        if link.exists() and link.resolve() != Path(source["parameters"]).resolve():
            raise ValueError("parameter link changed")
        if not link.exists():
            link.symlink_to(Path(source["parameters"]).resolve(), target_is_directory=True)
        tasks[name] = {"name": name, "runtime": source["tasks"][ligand_key]["runtime"],
            "molblock": Chem.MolToMolBlock(molecule), "bound_ligand_deck": ligand_deck,
            "radial_coordinates": rows, "candidate_kind": "graph-flagged conjugated NH2 only, not radial normalization",
            "diagnostics": diagnostics}
        for state, deck in {"ligand": replace_coordinates(ligand_deck, rows), "complex": transplant_ligand(complex_path.read_text(), ligand_deck, rows)}.items():
            path = parent/state/"input.inp"
            path.parent.mkdir(parents=True, exist_ok=True)
            if path.exists() and path.read_text() != deck:
                raise ValueError("refusing changed graph repair deck")
            path.write_text(deck)
            tasks[name][state+"_input"] = str(path)
            tasks[name][state+"_input_sha256"] = digest(path)
    frozen = {"schema": "graph-induced-conjugated-nh2-repair-control-v1", "target": target, "names": names,
        "tasks": tasks, "source_case_tasks": source["case_tasks"], "registry_sha256": digest(args.registry),
        "source_protocol_sha256": digest(args.protocol), "prior_sha256": digest(args.prior), "driver_sha256": digest(__file__),
        "executor_sha256": digest(Path(__file__).with_name("run_direction_preserving_hydrogen_control.py")),
        "graph_core_sha256": digest(Path(__file__).parents[1]/"src/fmotpiefinal/core/conjugated_amine_geometry.py"),
        "audit_sha256": digest(args.output/"graph_audit.json"), "labels_used_for_trigger": False,
        "controls": "No heavy, receptor/water, charge, bond-order or QM model changes; only newly distorted guanidino NH2 H rebuilt in conjugated plane at 120 degrees and original N-H lengths. H assignment by minimum displacement, not score.",
        "arms": ["matched_tpie", "original_total_strain_score", "graph_repaired_tpie", "graph_repaired_tpie_plus_local_strain"],
        "warning": "Development-only geometry candidate, not default; idealized NH2 geometry is not a QM stationary point or a proof of improved affinity prediction"}
    path = args.output/"frozen_protocol.json"
    if path.exists() and json.loads(path.read_text()) != frozen:
        raise ValueError("frozen graph experiment changed")
    write_json(path, frozen)
    results = {}
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context("spawn")) as pool:
        futures = [pool.submit(execute_geometry_pair, (task, str(args.workspace.resolve()))) for task in tasks.values()]
        for future in as_completed(futures):
            name, row = future.result()
            results[name] = row
            write_json(args.output/"progress.json", {"n_expected": len(tasks), "n_finished": len(results),
                "n_succeeded": sum(value["status"] == "SUCCEEDED" for value in results.values()), "rows": results})
            print(name, row["status"], row.get("error", ""), flush=True)
    complete = len(results) == len(names) and all(row["status"] == "SUCCEEDED" for row in results.values())
    report = {"status": "COMPLETE_VERIFIED" if complete else "INCOMPLETE", "rows": results, "scores": {}, "metrics": {},
              "n_expected": len(names), "n_succeeded": sum(row["status"] == "SUCCEEDED" for row in results.values()), "production_default_changed": False,
              "wall_seconds": time.monotonic()-start}
    if complete:
        truth = load_truth(args.workspace.resolve(), target)
        if set(truth) != set(names):
            raise ValueError("complete evaluation identity inventory required")
        for name, row in results.items():
            old = prior["scores"][name]
            bound = prior["rows"][source["case_tasks"][name]["ligand"]]["energy_hartree"]
            strain = old["ligand_strain"]-(bound-row["ligand"]["energy_hartree"])*HARTREE_TO_KCAL_MOL
            report["scores"][name] = {"matched_tpie": old["matched_tpie"], "original_total_strain_score": old["tpie_plus_ligand_strain"],
                "graph_repaired_tpie": row["complex"]["tpie"], "graph_repaired_local_strain": strain,
                "graph_repaired_tpie_plus_local_strain": row["complex"]["tpie"]+strain}
        report["metrics"] = {arm: metrics([truth[name] for name in sorted(names)], [-report["scores"][name][arm] for name in sorted(names)]) for arm in frozen["arms"]}
    write_json(args.output/"evaluation.json", report)
    print(report["status"], report["metrics"], flush=True)


if __name__ == "__main__":
    main()
