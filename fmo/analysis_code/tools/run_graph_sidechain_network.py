"""Freeze and run exactly one contact-shell extension, without score-guided search."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from copy import deepcopy
import json
import multiprocessing
from pathlib import Path
import time

from fmotpiefinal.core.sidechain_network import extend_one_contact_shell
from fmotpiefinal.core.structure import read_pdb
from run_cdk8_solution_cycle import digest, write_json
from run_graph_sidechain_transfer import execute
from run_local_sidechain_transfer import evaluate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--seed-run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("positive workers required")
    workspace, seed_run, output = (path.resolve() for path in (args.workspace, args.seed_run, args.output))
    seed_protocol = json.loads((seed_run / "frozen_protocol.json").read_text())
    seed_final = json.loads((seed_run / "final_status.json").read_text())
    if seed_final["status"] != "COMPLETE":
        raise ValueError("complete frozen seed control required")
    baseline = workspace / "FMOtpiefinal_HIV1lab/runs/all_accelerated_robust_20261020"
    output.mkdir(parents=True, exist_ok=True)
    plan_path = output / "frozen_protocol.json"
    if plan_path.exists():
        plan = json.loads(plan_path.read_text())
    else:
        cohorts = deepcopy(seed_protocol["cohorts"])
        for target, cohort in cohorts.items():
            seed = json.loads((seed_run / target / "sidechain_decision_graph.json").read_text())
            decision = extend_one_contact_shell(read_pdb(cohort["protein"]), seed)
            write_json(output / target / "sidechain_decision_graph.json", decision)
            cohort["local_residues"] = decision["selected_residues"]
            cohort["decision_plan_sha256"] = decision["plan_sha256"]
            print("PLAN", target, cohort["local_residues"], flush=True)
        source = Path(__file__).parents[1] / "src/fmotpiefinal"
        files = [Path(__file__), Path(__file__).with_name("run_graph_sidechain_transfer.py"),
            Path(__file__).with_name("run_local_sidechain_transfer.py"), source / "core/sidechain_network.py",
            source / "core/minimization.py", source / "pipeline/target_stages.py"]
        plan = {"schema": "one-hop-sidechain-full-cohort-control-v1", "targets": list(cohorts),
            "cohorts": cohorts, "seed_protocol_sha256": digest(seed_run / "frozen_protocol.json"),
            "selection_uses_scores": False, "labels_used_for_selection": False,
            "scope": "13 Ross cohorts, no OpenBind", "production_default_changed": False,
            "network_control": "exactly one 4A geometric sidechain contact shell; only baseline binding-region nodes",
            "code_sha256": {str(path.resolve()): digest(path) for path in files}}
        write_json(plan_path, plan)
    for path, expected in plan["code_sha256"].items():
        if digest(path) != expected:
            raise ValueError("network control code changed")
    for cohort in plan["cohorts"].values():
        for path, expected in cohort["input_sha256"].items():
            if digest(path) != expected:
                raise ValueError("source changed")
    results = {target: {} for target in plan["targets"]}
    started = time.monotonic()
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context("spawn")) as pool:
        futures = [pool.submit(execute, (cohort, name, str(output)))
            for cohort in plan["cohorts"].values() for name in cohort["names"]]
        for future in as_completed(futures):
            target, name, result = future.result()
            results[target][name] = result
            print("CASE", target, name, result["status"], result.get("tpie", result.get("error")), flush=True)
            write_json(output / "progress.json", results)
            if len(results[target]) == len(plan["cohorts"][target]["names"]):
                evaluate(workspace, baseline, results, output)
    evaluation = evaluate(workspace, baseline, results, output)
    successful = sum(group["n_succeeded"] for group in evaluation["groups"].values())
    expected = sum(len(cohort["names"]) for cohort in plan["cohorts"].values())
    write_json(output / "final_status.json", {"status": "COMPLETE" if successful == expected else "INCOMPLETE_WITH_FAILURES",
        "n_requested": expected, "n_succeeded": successful, "wall_seconds": time.monotonic()-started,
        "evaluation_sha256": digest(output / "evaluation.json")})


if __name__ == "__main__":
    main()
