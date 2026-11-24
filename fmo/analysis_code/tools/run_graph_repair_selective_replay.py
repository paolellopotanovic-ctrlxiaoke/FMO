"""Replay graph NH2 coordinates with the unchanged production tPIE Hamiltonian."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import io
import json
import multiprocessing
from pathlib import Path

from fmotpiefinal.core.fmo_convergence import final_fmo_convergence
from fmotpiefinal.core.parser import parse_pie
from fmotpiefinal.core.runner import run_gamess
from report_local_sidechain_transfer import load_truth
from run_cdk8_solution_cycle import digest, metrics, write_json
from run_hydrogen_normalized_strain import transplant_ligand
from run_local_sidechain_transfer import directory_key


def worker(payload):
    task = payload
    path = Path(task["input"])
    try:
        if digest(path) != task["input_sha256"]:
            raise ValueError("frozen selective deck changed")
        cached = path.parent/"run/result.json"
        runtime = task["runtime"]
        if cached.exists():
            result = json.loads(cached.read_text())
        else:
            result = run_gamess(path, path.parent/"run", gamess_dir=runtime["gamess_dir"],
                library_dir=runtime["library_dir"], ligand_fragment=runtime["ligand_fragment"],
                ncores=1, timeout_seconds=2400, parameter_dir=task["parameters"])
        if result["status"] != "SUCCEEDED" or result["input_sha256"] != digest(path):
            raise ValueError(f"selective production QM did not pass: {result.get('status')}")
        if result["executable_sha256"] != digest(Path(runtime["gamess_dir"])/"gamess.00.x"):
            raise ValueError("executable changed")
        log = path.parent/"run/gamess.log"
        convergence = final_fmo_convergence(path.read_text(), io.StringIO(log.read_text()))
        count = runtime["ligand_fragment"]
        pairs = [pair for pair in parse_pie(log) if count in (pair.fragment_i, pair.fragment_j)]
        if len(pairs) != count-1 or abs(sum(pair.value_kcal_mol for pair in pairs)-result["tpie_kcal_mol"]) > 1.e-8:
            raise ValueError("independent ligand pair inventory/sum failed")
        row = {"status": "SUCCEEDED", "tpie": result["tpie_kcal_mol"], "input_sha256": digest(path),
            "log_sha256": digest(log), "convergence": convergence, "executable_sha256": result["executable_sha256"]}
    except Exception as error:
        row = {"status": "FAILED", "error": f"{type(error).__name__}: {error}"}
    write_json(path.parent/"summary.json", row)
    return task["name"], row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph-run", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=18)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("positive workers required")
    graph = json.loads((args.graph_run/"frozen_protocol.json").read_text())
    registry = json.loads(args.registry.read_text())
    cohort = registry["cohorts"][graph["target"]]
    if set(graph["names"]) != set(cohort["names"]):
        raise ValueError("complete matched graph/production identities required")
    tasks = {}
    for name in graph["names"]:
        source = Path(cohort["sources"][name])
        if digest(source) != cohort["input_sha256"][str(source)]:
            raise ValueError("original production input changed")
        candidate = graph["tasks"][name]
        deck = transplant_ligand(source.read_text(), candidate["bound_ligand_deck"], candidate["radial_coordinates"])
        path = args.output.resolve()/cohort["target"]/directory_key(name)/"input.inp"
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists() and path.read_text() != deck:
            raise ValueError("refusing changed production replay")
        path.write_text(deck)
        tasks[name] = {"name": name, "input": str(path), "input_sha256": digest(path),
            "runtime": cohort["runtime"][name], "parameters": cohort["parameters"], "original_sha256": digest(source)}
    frozen = {"schema": "graph-nh2-original-production-selective-replay-v1", "target": graph["target"], "names": graph["names"],
        "tasks": tasks, "registry_sha256": digest(args.registry), "graph_protocol_sha256": digest(args.graph_run/"frozen_protocol.json"),
        "driver_sha256": digest(__file__), "labels_used_for_geometry": False,
        "scope": "Only graph-flagged NH2 H coordinates changed. Original selective FMO/IDISP/PCM and all noncoordinate groups preserved, no strain added."}
    path = args.output/"frozen_protocol.json"
    if path.exists() and json.loads(path.read_text()) != frozen:
        raise ValueError("frozen selective replay changed")
    write_json(path, frozen)
    rows = {}
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context("spawn")) as pool:
        futures = [pool.submit(worker, task) for task in tasks.values()]
        for future in as_completed(futures):
            name, row = future.result()
            rows[name] = row
            write_json(args.output/"progress.json", {"n_expected": len(tasks), "n_finished": len(rows),
                "n_succeeded": sum(value["status"] == "SUCCEEDED" for value in rows.values()), "rows": rows})
            print(name, row["status"], row.get("error", ""), flush=True)
    complete = len(rows) == len(tasks) and all(row["status"] == "SUCCEEDED" for row in rows.values())
    report = {"status": "COMPLETE_VERIFIED" if complete else "INCOMPLETE", "rows": rows, "metrics": {}, "production_default_changed": False}
    if complete:
        truth = load_truth(args.workspace.resolve(), graph["target"])
        if set(truth) != set(tasks):
            raise ValueError("complete truth inventory required")
        names = sorted(tasks)
        report["baseline_metrics"] = metrics([truth[name] for name in names], [-cohort["baseline"][name]["tpie"] for name in names])
        report["metrics"] = metrics([truth[name] for name in names], [-rows[name]["tpie"] for name in names])
    write_json(args.output/"evaluation.json", report)
    print(report["status"], report.get("baseline_metrics"), report["metrics"], flush=True)


if __name__ == "__main__":
    main()
