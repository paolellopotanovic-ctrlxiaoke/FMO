"""Two preregistered broader controls, no target-wise or score-guided tuning."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from copy import deepcopy
import json
import multiprocessing
from pathlib import Path
import time

import numpy as np
from rdkit import Chem

from fmotpiefinal.core.bound_geometry import load_bound_ligand
from fmotpiefinal.core.sidechain_motion import add_displaced_clash_seeds
from fmotpiefinal.core.sidechain_network import extend_one_contact_shell
from fmotpiefinal.core.structure import read_pdb
from run_cdk8_solution_cycle import digest, write_json
from run_graph_sidechain_transfer import execute
from run_local_sidechain_transfer import evaluate


ARMS = ("strict_seeds_two_hops", "motion_seeds_two_hops")


def run_case(payload):
    arm, cohort, name, root = payload
    target, name, result = execute((cohort, name, str(Path(root) / arm)))
    return arm, target, name, result


def register(seed_run, output):
    source = json.loads((seed_run / "frozen_protocol.json").read_text())
    source_final = json.loads((seed_run / "final_status.json").read_text())
    if source_final["status"] != "COMPLETE":
        raise ValueError("complete seed experiment required")
    cohorts = {arm: deepcopy(source["cohorts"]) for arm in ARMS}
    summaries = {}
    for target, cohort in source["cohorts"].items():
        protein = read_pdb(cohort["protein"])
        seed = json.loads((seed_run / target / "sidechain_decision_graph.json").read_text())
        molecules = {molecule.GetProp("_Name").strip(): molecule for molecule in
            Chem.SDMolSupplier(cohort["sdf"], removeHs=False) if molecule is not None}
        if set(molecules) != set(cohort["names"]):
            raise ValueError("full ligand inventory mismatch")
        bound = {name: np.asarray([atom.xyz for atom in load_bound_ligand(cohort["sdf"], name, cohort["complexes"][name]).atoms])
                 for name in cohort["names"]}
        motion = add_displaced_clash_seeds(protein, molecules, bound, seed)
        write_json(output / target / "motion_seed_graph.json", motion)
        summaries[target] = {"n": len(cohort["names"]), "strict_seeds": seed["selected_residues"],
            "motion_seeds": motion["selected_residues"], "n_motion_events": len(motion["motion_events"])}
        for arm, initial in zip(ARMS, (seed, motion), strict=True):
            decision = extend_one_contact_shell(protein, initial, n_hops=2)
            write_json(output / arm / target / "sidechain_decision_graph.json", decision)
            cohorts[arm][target]["local_residues"] = decision["selected_residues"]
            cohorts[arm][target]["decision_plan_sha256"] = decision["plan_sha256"]
            summaries[target][arm] = decision["selected_residues"]
        print("PLAN", target, json.dumps(summaries[target]), flush=True)
    files = [Path(__file__), Path(__file__).with_name("run_graph_sidechain_transfer.py"),
        Path(__file__).with_name("run_local_sidechain_transfer.py")]
    source_root = Path(__file__).parents[1] / "src/fmotpiefinal"
    files += [source_root / name for name in ("core/sidechain_motion.py", "core/sidechain_network.py",
        "core/sidechain_policy.py", "core/minimization.py", "pipeline/target_stages.py")]
    code = {str(path.resolve()): digest(path) for path in files}
    plan = {"schema": "adaptive-sidechain-two-control-v1", "arms": list(ARMS), "cohorts": cohorts,
        "targets": list(source["cohorts"]), "selection_summary": summaries, "code_sha256": code,
        "seed_protocol_sha256": digest(seed_run / "frozen_protocol.json"),
        "thresholds": {"additional_trigger": "original heavy overlap >=0.6A plus ligand heavy displacement >=0.3A",
            "network_hops": 2, "network_contact_cutoff_A": 4.0, "polar_original_distance_below_A": 2.2},
        "design": "Control A changes only network extent. Control B additionally changes seed evidence. No further parameter search.",
        "scores_used_for_selection": False, "labels_used_for_selection": False,
        "development_status": "retrospective mechanism-driven development, not a pristine external validation",
        "production_default_changed": False, "openbind_requested": False}
    write_json(output / "frozen_protocol.json", plan)
    for arm in ARMS:
        write_json(output / arm / "frozen_protocol.json", {"schema": "adaptive-sidechain-arm-v1",
            "arm": arm, "targets": plan["targets"], "cohorts": cohorts[arm], "code_sha256": code,
            "parent_protocol_sha256": digest(output / "frozen_protocol.json"),
            "production_default_changed": False})
    return plan


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", required=True, type=Path)
    parser.add_argument("--seed-run", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--register-only", action="store_true")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("workers must be positive")
    workspace, seed_run, output = (path.resolve() for path in (args.workspace, args.seed_run, args.output))
    output.mkdir(parents=True, exist_ok=True)
    plan_path = output / "frozen_protocol.json"
    plan = json.loads(plan_path.read_text()) if plan_path.exists() else register(seed_run, output)
    for path, expected in plan["code_sha256"].items():
        if digest(path) != expected:
            raise ValueError("frozen adaptive implementation changed")
    for cohort in plan["cohorts"][ARMS[0]].values():
        for path, expected in cohort["input_sha256"].items():
            if digest(path) != expected:
                raise ValueError("frozen input changed")
    if args.register_only:
        return
    results = {arm: {target: {} for target in plan["targets"]} for arm in ARMS}
    started = time.monotonic()
    baseline = workspace / "FMOtpiefinal_HIV1lab/runs/all_accelerated_robust_20261020"
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context("spawn")) as pool:
        futures = [pool.submit(run_case, (arm, cohort, name, str(output)))
            for arm in ARMS for cohort in plan["cohorts"][arm].values() for name in cohort["names"]]
        for future in as_completed(futures):
            arm, target, name, result = future.result()
            results[arm][target][name] = result
            print("CASE", arm, target, name, result["status"], result.get("tpie", result.get("error")), flush=True)
            write_json(output / arm / "progress.json", results[arm])
            if len(results[arm][target]) == len(plan["cohorts"][arm][target]["names"]):
                evaluate(workspace, baseline, results[arm], output / arm)
    final = {}
    for arm in ARMS:
        evaluation = evaluate(workspace, baseline, results[arm], output / arm)
        expected = sum(len(cohort["names"]) for cohort in plan["cohorts"][arm].values())
        successful = sum(group["n_succeeded"] for group in evaluation["groups"].values())
        final[arm] = {"status": "COMPLETE" if successful == expected else "INCOMPLETE_WITH_FAILURES",
            "n_requested": expected, "n_succeeded": successful, "wall_seconds": time.monotonic()-started,
            "evaluation_sha256": digest(output / arm / "evaluation.json")}
        write_json(output / arm / "final_status.json", final[arm])
    write_json(output / "final_status.json", final)


if __name__ == "__main__":
    main()
