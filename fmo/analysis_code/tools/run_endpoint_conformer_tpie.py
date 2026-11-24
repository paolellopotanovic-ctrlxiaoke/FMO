"""Fixed-receptor, pose-removed ligand conformer control; not a free energy."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import re

import numpy as np

from analyze_openfe_endpoints import fit_coordinates
from run_cdk8_solution_cycle import group
from run_matched_openfe import digest, write_json


def conformer_deck(source, elements, coordinates):
    coordinates = np.asarray(coordinates, dtype=float)
    if coordinates.shape != (len(elements), 3) or not np.isfinite(coordinates).all():
        raise ValueError("finite complete ligand coordinates required")
    section = group(source, "FMOXYZ")
    rows = section[1].strip().splitlines()
    ligand = rows[-len(elements):]
    if len(rows) <= len(elements) or [row.split()[0].upper() for row in ligand] != [value.upper() for value in elements]:
        raise ValueError("ligand-last atom inventory mismatch")
    original = np.asarray([row.split()[2:] for row in ligand], dtype=float)
    heavy = np.asarray([element.upper() != "H" for element in elements])
    if heavy.sum() < 3:
        raise ValueError("three heavy atoms required for alignment")
    aligned = fit_coordinates(coordinates[heavy], original[heavy], coordinates)
    changed = [f"{row.split()[0]} {row.split()[1]} " + " ".join(f"{value:.10f}" for value in xyz)
               for row, xyz in zip(ligand, aligned, strict=True)]
    replacement = "$FMOXYZ\n" + "\n".join(rows[:-len(elements)] + changed) + "\n $END"
    deck = source[:section.start()] + replacement + source[section.end():]
    receptor = np.asarray([row.split()[2:] for row in rows[:-len(elements)]], dtype=float)
    receptor_heavy = np.asarray([row.split()[0].upper() != "H" for row in rows[:-len(elements)]])
    nearest = float(np.linalg.norm(receptor[receptor_heavy, None] - aligned[heavy][None], axis=-1).min())
    return deck, {"minimum_heavy_contact_angstrom": nearest,
                  "severe_heavy_overlap": nearest < 1.2,
                  "ligand_internal_rmsd_angstrom": float(np.sqrt(np.mean(np.sum((aligned[heavy] - original[heavy]) ** 2, axis=1))))}


def execute(task, timeout):
    from fmotpiefinal.core.runner import run_gamess

    directory = Path(task["directory"])
    directory.mkdir(parents=True, exist_ok=False)
    source = Path(task["source"])
    if digest(source) != task["source_sha256"]:
        raise ValueError("frozen reference deck changed")
    text = source.read_text()
    if task["kind"] == "baseline_replay":
        deck, geometry = text, {}
    else:
        deck, geometry = conformer_deck(text, task["elements"], task["xyz_angstrom"])
    path = directory / "input.inp"
    path.write_text(deck)
    if geometry.get("severe_heavy_overlap"):
        result = {"status": "FAILED_GEOMETRY_CONTROL", "reason": "transplanted conformer clashes with fixed receptor"}
    else:
        count = int(re.search(r"\bNFRAG=(\d+)", deck)[1])
        result = run_gamess(path, directory / "run", gamess_dir=task["runtime"]["gamess_dir"],
            library_dir=task["runtime"]["library_dir"], ligand_fragment=count,
            parameter_dir=task["parameters"], ncores=1, timeout_seconds=timeout)
    row = {key: task[key] for key in ("target", "name", "kind", "directory")}
    row.update(iteration=task.get("iteration"), trajectory=task.get("trajectory"),
               geometry=geometry, result=result)
    write_json(directory / "summary.json", row)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout", type=int, default=1800)
    args = parser.parse_args()
    if min(args.workers, args.timeout) < 1:
        parser.error("positive budgets required")
    root = args.workspace.resolve() / "FMOtpiefinal_HIV1lab/runs"
    run, output = args.run.resolve(), args.output.resolve()
    registry = root / "graph_sidechain_transfer_20261021/frozen_protocol.json"
    cohorts = json.loads(registry.read_text())["cohorts"]
    endpoint_plan = run / "final_endpoint_audit_20261023/qm/frozen_protocol.json"
    plan = json.loads(endpoint_plan.read_text())
    if plan["labels_used"]:
        raise ValueError("label-independent frame plan required")
    tasks = []
    baseline = {}
    for original in plan["tasks"]:
        if original["leg"] != "bound" or original["iteration"] not in (400, 1200, 2000):
            continue
        target, name = original["target"], original["name"]
        source = (root / f"main_response_regression_20261021/CDK8/sidechain_response/qm/{name}/input.inp"
                  if target == "CDK8" else Path(cohorts[target]["sources"][name]))
        task = {**original, "source": str(source), "source_sha256": digest(source), "kind": "conformer",
                "directory": str(output / target / name / Path(original["trajectory"]).parent.name / f'iteration_{original["iteration"]}')}
        tasks.append(task)
        baseline[target, name] = {**task, "kind": "baseline_replay", "directory": str(output / target / name / "baseline_replay")}
    tasks = list(baseline.values()) + tasks
    if len(tasks) != 40:
        raise ValueError("four baseline replays and 36 fixed conformer tasks required")
    output.mkdir(parents=True, exist_ok=False)
    write_json(output / "frozen_protocol.json", {"tasks": tasks, "driver_sha256": digest(__file__),
        "endpoint_plan_sha256": digest(endpoint_plan), "registry_sha256": digest(registry),
        "labels_used": False, "main_pipeline_changed": False,
        "scope": "Ligand internal conformers only: heavy-atom self-fit removes MD pose; source receptor, caps, charges, PCM and fragments unchanged. Not same-frame full-complex FMO or a free energy."})
    rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(execute, task, args.timeout): task for task in tasks}
        for future in as_completed(futures):
            task = futures[future]
            try:
                row = future.result()
            except Exception as error:
                row = {"target": task["target"], "name": task["name"], "kind": task["kind"],
                       "directory": task["directory"], "result": {"status": "FAILED", "error": f"{type(error).__name__}: {error}"}}
            rows.append(row)
            write_json(output / "progress.json", rows)
            print(len(rows), len(tasks), task["target"], task["name"], row["result"]["status"], flush=True)
    succeeded = sum(row["result"]["status"] == "SUCCEEDED" for row in rows)
    write_json(output / "status.json", {"status": "COMPLETE" if succeeded == len(tasks) else "INCOMPLETE_WITH_FAILURES",
        "n_requested": len(tasks), "n_succeeded": succeeded, "rows": rows, "cohort_ranking_available": False})


if __name__ == "__main__":
    main()
