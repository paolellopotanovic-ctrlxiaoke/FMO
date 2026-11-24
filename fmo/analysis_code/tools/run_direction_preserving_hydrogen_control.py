"""Preserve bound X-H directions while transplanting QM reference bond lengths."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import multiprocessing
from pathlib import Path
import sys
import time

import numpy as np
from rdkit import Chem

from audit_cycle_ligand_chemistry import geometry_screen, stereochemistry
from fmotpiefinal.core.parser import parse_pie
from fmotpiefinal.core.runner import charge_aware_pair_outliers
from fmotpiefinal.core.solution_cycle import HARTREE_TO_KCAL_MOL
from report_local_sidechain_transfer import load_truth
from run_cdk8_bound_qm import replace_coordinates
from run_cdk8_solution_cycle import digest, metrics, write_json
from run_cohort_solution_cycle import execute
from run_hydrogen_normalized_strain import fixed_heavy_displacement, transplant_ligand, xyz_rows
from run_local_sidechain_transfer import directory_key


def radial_coordinates(molecule, initial_rows, qm_rows):
    fixed_heavy_displacement(initial_rows, qm_rows)
    before = np.asarray([row.split()[2:] for row in initial_rows], float)
    after = np.asarray([row.split()[2:] for row in qm_rows], float)
    if molecule.GetNumAtoms() != len(before):
        raise ValueError("chemical atom inventory mismatch")
    points = before.copy()
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() != 1:
            continue
        if atom.GetDegree() != 1 or atom.GetNeighbors()[0].GetAtomicNum() == 1:
            raise ValueError("expected a hydrogen bound to one heavy atom")
        index, parent = atom.GetIdx(), atom.GetNeighbors()[0].GetIdx()
        direction = before[index]-before[parent]
        distance = float(np.linalg.norm(direction))
        target = float(np.linalg.norm(after[index]-after[parent]))
        if min(distance, target) < .6 or max(distance, target) > 1.45:
            raise ValueError("unphysical hydrogen bond length")
        points[index] = before[parent]+direction*(target/distance)
    rows = [f"{row.split()[0]} {row.split()[1]} {point[0]:.10f} {point[1]:.10f} {point[2]:.10f}"
            for row, point in zip(initial_rows, points, strict=True)]
    fixed_heavy_displacement(initial_rows, rows)
    return rows


def worker(payload):
    task, workspace = payload
    parent = Path(task["ligand_input"]).parent.parent
    try:
        molecule = Chem.MolFromMolBlock(task["molblock"], removeHs=False)
        initial = np.asarray([row.split()[2:] for row in xyz_rows(task["bound_ligand_deck"])], float)
        final = np.asarray([row.split()[2:] for row in task["radial_coordinates"]], float)
        chemistry = geometry_screen(molecule, final, stereochemistry(molecule, initial))
        if not chemistry["passed_screen"]:
            raise ValueError("direction-preserving geometry failed chemical screening")
        summaries = {}
        for state in ("ligand", "complex"):
            path = Path(task[state+"_input"])
            if digest(path) != task[state+"_input_sha256"]:
                raise ValueError("frozen input changed")
            actual = {"key": task["name"]+"/"+state, "state": state, "input": str(path),
                      "input_sha256": digest(path), "runtime": task["runtime"]}
            _, summary = execute((actual, workspace))
            if summary["status"] != "SUCCEEDED":
                raise ValueError(f"{state}: {summary.get('error')}")
            summaries[state] = summary
        path = Path(task["complex_input"])
        fragment = task["runtime"]["ligand_fragment"]
        pairs = [pair for pair in parse_pie(path.parent/"gamess.log") if fragment in (pair.fragment_i, pair.fragment_j)]
        outliers, _ = charge_aware_pair_outliers(pairs, path.read_text(), 100.)
        if len(pairs) != fragment-1 or outliers:
            raise ValueError("original charge-aware ligand PIEDA gate failed")
        result = {"status": "SUCCEEDED", **summaries, "chemistry": chemistry}
    except Exception as error:
        result = {"status": "FAILED", "error": f"{type(error).__name__}: {error}"}
    write_json(parent/"summary.json", result)
    return task["name"], result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hydrogen-run", type=Path, required=True)
    parser.add_argument("--prior", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=16)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("positive workers required")
    start = time.monotonic()
    source = json.loads((args.hydrogen_run/"frozen_protocol.json").read_text())
    hydrogen = json.loads((args.hydrogen_run/"evaluation.json").read_text())
    prior = json.loads(args.prior.read_text())
    if hydrogen["status"] != "COMPLETE_VERIFIED" or digest(args.prior) != source["evaluation_sha256"]:
        raise ValueError("complete verified hydrogen and prior references required")
    tasks = {}
    for name in source["names"]:
        old = source["tasks"][name]
        path = Path(old["input"])
        if digest(path.parent/"gamess.log") != hydrogen["rows"][name]["hydrogen"]["log_sha256"]:
            raise ValueError("source constrained-optimization log changed")
        molecule = Chem.MolFromMolBlock(old["molblock"], removeHs=False)
        rows = radial_coordinates(molecule, xyz_rows(old["bound_ligand_deck"]), hydrogen["rows"][name]["hydrogen"]["optimized_coordinates"])
        template = Path(old["complex_template"])
        if digest(template) != old["complex_template_sha256"]:
            raise ValueError("bound complex template changed")
        decks = {"ligand": replace_coordinates(old["bound_ligand_deck"], rows),
                 "complex": transplant_ligand(template.read_text(), old["bound_ligand_deck"], rows)}
        parent = args.output.resolve()/source["target"]/directory_key(name)
        parent.mkdir(parents=True, exist_ok=True)
        link = parent/"params"
        parameters = path.parent.parent/"params"
        if link.exists() and link.resolve() != parameters.resolve():
            raise ValueError("parameter directory changed")
        if not link.exists():
            link.symlink_to(parameters.resolve(), target_is_directory=True)
        tasks[name] = {"name": name, "runtime": old["runtime"], "molblock": old["molblock"],
                       "bound_ligand_deck": old["bound_ligand_deck"], "radial_coordinates": rows}
        for state, deck in decks.items():
            destination = parent/state/"input.inp"
            destination.parent.mkdir(parents=True, exist_ok=True)
            if destination.exists() and destination.read_text() != deck:
                raise ValueError("refusing changed radial-geometry input")
            destination.write_text(deck)
            tasks[name][state+"_input"] = str(destination)
            tasks[name][state+"_input_sha256"] = digest(destination)
    frozen = {"schema": "direction-preserving-hydrogen-bondlength-control-v1", "target": source["target"],
        "names": source["names"], "tasks": tasks, "driver_sha256": digest(__file__),
        "hydrogen_protocol_sha256": digest(args.hydrogen_run/"frozen_protocol.json"),
        "hydrogen_evaluation_sha256": digest(args.hydrogen_run/"evaluation.json"), "prior_sha256": digest(args.prior),
        "source_case_tasks": source["source_case_tasks"], "affinity_labels_used_for_geometry": False,
        "dependencies": {name: digest(Path(__file__).with_name(name)) for name in ("run_cohort_solution_cycle.py", "run_hydrogen_normalized_strain.py", "audit_cycle_ligand_chemistry.py")},
        "arms": ["matched_tpie", "original_total_strain_score", "radial_tpie", "radial_tpie_plus_residual_strain",
                 "original_tpie_plus_radial_component", "original_tpie_plus_radial_residual"],
        "protocol": "All heavy/receptor/water coordinates and X-H directions fixed. Per-bond length from same-model constrained-QM reference, no fitted or cohort-average bond lengths. One L energy and one full AB energy per case. All coefficients one.",
        "warning": "Post-hoc development diagnostic; reconstructed radial geometry is not claimed to be a stationary point or thermodynamic state"}
    path = args.output/"frozen_protocol.json"
    if path.exists() and json.loads(path.read_text()) != frozen:
        raise ValueError("frozen radial experiment changed")
    write_json(path, frozen)
    results = {}
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context("spawn")) as pool:
        futures = [pool.submit(worker, (task, str(args.workspace.resolve()))) for task in tasks.values()]
        for future in as_completed(futures):
            name, row = future.result()
            results[name] = row
            write_json(args.output/"progress.json", {"n_expected": len(tasks), "n_finished": len(results),
                "n_succeeded": sum(item["status"] == "SUCCEEDED" for item in results.values()), "rows": results})
            print(name, row["status"], row.get("error", ""), flush=True)
    complete = set(results) == set(source["names"]) and all(row["status"] == "SUCCEEDED" for row in results.values())
    report = {"status": "COMPLETE_VERIFIED" if complete else "INCOMPLETE", "rows": results,
              "n_expected": len(tasks), "n_succeeded": sum(row["status"] == "SUCCEEDED" for row in results.values()),
              "scores": {}, "metrics": {}, "production_default_changed": False, "wall_seconds": time.monotonic()-start}
    if complete:
        truth = load_truth(args.workspace.resolve(), source["target"])
        names = sorted(results)
        if set(truth) != set(names):
            raise ValueError("complete evaluation identity inventory required")
        for name, row in results.items():
            old = prior["scores"][name]
            bound = prior["rows"][source["source_case_tasks"][name]["ligand"]]["energy_hartree"]
            change = (bound-row["ligand"]["energy_hartree"])*HARTREE_TO_KCAL_MOL
            residual = old["ligand_strain"]-change
            report["scores"][name] = {"matched_tpie": old["matched_tpie"], "original_total_strain_score": old["tpie_plus_ligand_strain"],
                "radial_component": change, "radial_residual": residual, "radial_tpie": row["complex"]["tpie"],
                "radial_tpie_plus_residual_strain": row["complex"]["tpie"]+residual,
                "original_tpie_plus_radial_component": old["matched_tpie"]+change,
                "original_tpie_plus_radial_residual": old["matched_tpie"]+residual}
        report["metrics"] = {arm: metrics([truth[name] for name in names], [-report["scores"][name][arm] for name in names]) for arm in frozen["arms"]}
    write_json(args.output/"evaluation.json", report)
    print(report["status"], report["metrics"], flush=True)


if __name__ == "__main__":
    main()
