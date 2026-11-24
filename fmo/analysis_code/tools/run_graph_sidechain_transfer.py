"""Execute one frozen geometry-triggered sidechain policy on complete cohorts."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from copy import deepcopy
import json
import multiprocessing
from pathlib import Path
import time

from fmotpiefinal.sidechains import audit_sidechains
from run_cdk8_solution_cycle import digest, write_json
from run_local_sidechain_transfer import directory_key, evaluate, prepare_job, register_cohort, score_job


def execute(payload):
    cohort, name, output = payload
    target = cohort["target"]
    parent = Path(output) / target / directory_key(name)
    if not cohort["local_residues"]:
        baseline = cohort["baseline"][name]
        result = {"status": "SUCCEEDED", "tpie": baseline["tpie"], "baseline_tpie": baseline["tpie"],
            "execution_path": "UNCHANGED_VERIFIED_BASELINE_NO_STRUCTURAL_TRIGGER",
            "new_qm_executed": False, "input_sha256": baseline["input_sha256"],
            "log_sha256": baseline["log_sha256"], "pair_count": baseline["n_pairs"],
            "sidechain_plan_sha256": cohort["decision_plan_sha256"]}
        write_json(parent / "summary.json", result)
        return target, name, result
    summary = parent / "summary.json"
    if summary.exists():
        previous = json.loads(summary.read_text())
        if previous.get("status") == "SUCCEEDED" and previous.get("sidechain_plan_sha256") == cohort["decision_plan_sha256"]:
            return target, name, previous
        raise ValueError("existing failed or unregistered output requires explicit diagnostic continuation")
    preparation = parent / "preparation_status.json"
    if not preparation.exists():
        _, _, prepared = prepare_job((cohort, name, output))
        if prepared["status"] != "SUCCEEDED":
            result = {**prepared, "execution_path": "GRAPH_SELECTED_PREPARATION_FAILED",
                "sidechain_plan_sha256": cohort["decision_plan_sha256"]}
            write_json(summary, result)
            return target, name, result
    target, name, result = score_job(cohort, name, Path(output), 2)
    result.update(execution_path="GRAPH_SELECTED_LOCAL_RELAX_AND_QM", new_qm_executed=True,
        sidechain_plan_sha256=cohort["decision_plan_sha256"])
    write_json(summary, result)
    return target, name, result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--registered-transfer", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--register-only", action="store_true")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("positive worker count required")
    workspace, baseline, registered, output = (path.resolve() for path in
        (args.workspace, args.baseline, args.registered_transfer, args.output))
    output.mkdir(parents=True, exist_ok=True)
    plan_path = output / "frozen_protocol.json"
    if plan_path.exists():
        plan = json.loads(plan_path.read_text())
    else:
        original = json.loads((registered / "frozen_protocol.json").read_text())
        cohorts = deepcopy(original["cohorts"])
        cohorts["CDK8"] = register_cohort(baseline, "CDK8")
        summaries = {}
        for target, cohort in cohorts.items():
            decision = audit_sidechains(cohort["protein"], cohort["sdf"], cohort["complexes"])
            directory = output / target
            directory.mkdir(exist_ok=True)
            write_json(directory / "sidechain_decision_graph.json", decision)
            cohort["uniform_4A_residues"] = cohort["local_residues"]
            cohort["local_residues"] = decision["selected_residues"]
            cohort["decision_plan_sha256"] = decision["plan_sha256"]
            summaries[target] = {"n": len(cohort["names"]), "uniform_4A_count": len(cohort["uniform_4A_residues"]),
                "selected_count": len(cohort["local_residues"]), "selected_residues": cohort["local_residues"],
                "action": decision["action"], "protected_with_trigger": [node["id"] for node in decision["nodes"]
                    if node["action"] == "REVIEW_PROTECTED_NETWORK_KEEP_FIXED"]}
            print("PLAN", target, json.dumps(summaries[target]), flush=True)
        code = [Path(__file__), Path(__file__).with_name("run_local_sidechain_transfer.py")]
        source = Path(__file__).parents[1] / "src/fmotpiefinal"
        code += [source / "core/sidechain_policy.py", source / "sidechains.py",
                 source / "core/minimization.py", source / "pipeline/target_stages.py"]
        plan = {"schema": "graph-sidechain-full-cohort-control-v1", "cohorts": cohorts,
            "targets": list(cohorts), "selection_summary": summaries, "code_sha256": {str(path.resolve()): digest(path) for path in code},
            "source_registration_sha256": digest(registered / "frozen_protocol.json"),
            "labels_used_for_selection": False, "weights_fitted": False,
            "selection_uses_scores": False, "production_default_changed": False,
            "water_H_optimization_added": False, "scope": "13 Ross cohorts; OpenBind explicitly excluded",
            "conditional_execution": "No trigger: verified baseline unchanged. Trigger: same selected sidechains and k2.5/k10 for entire target.",
            "policy_development_status": "retrospective method development; not independent prospective validation"}
        write_json(plan_path, plan)
    for path, expected in plan["code_sha256"].items():
        if digest(path) != expected:
            raise ValueError("frozen graph experiment implementation changed")
    for target, cohort in plan["cohorts"].items():
        for path, expected in cohort["input_sha256"].items():
            if digest(path) != expected:
                raise ValueError(f"frozen source changed: {target}")
        decision = json.loads((output / target / "sidechain_decision_graph.json").read_text())
        if decision["plan_sha256"] != cohort["decision_plan_sha256"]:
            raise ValueError("decision plan identity changed")
    if args.register_only:
        return
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
    expected = sum(len(cohort["names"]) for cohort in plan["cohorts"].values())
    successful = sum(group["n_succeeded"] for group in evaluation["groups"].values())
    write_json(output / "final_status.json", {"status": "COMPLETE" if successful == expected else "INCOMPLETE_WITH_FAILURES",
        "n_requested": expected, "n_succeeded": successful, "wall_seconds": time.monotonic()-started,
        "evaluation_sha256": digest(output / "evaluation.json")})


if __name__ == "__main__":
    main()
