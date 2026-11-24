"""Matched full-cohort input-heavy-fixed hydrogen optimization control."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from copy import deepcopy
import json
import multiprocessing
from pathlib import Path
import time

from run_cdk8_solution_cycle import digest, write_json
from run_local_sidechain_transfer import directory_key, evaluate, score_job


TARGETS = ("CHK1", "BACE", "Galectin3", "CDK8")


def minimization_arguments(control):
    arguments = {"restraint_k_kj_mol_a2": 10.0, "move_protein": False}
    if control == "input_heavy_fixed":
        arguments["fix_ligand_heavy"] = True
    elif control == "hydrogens_unrestrained":
        arguments["restrain_ligand_hydrogens"] = False
    else:
        raise ValueError("unknown frozen ligand geometry control")
    return arguments


def execute(payload):
    cohort, name, output, control = payload
    from fmotpiefinal.core.minimization import minimize_complex
    parent = Path(output) / cohort["target"] / directory_key(name)
    started = time.monotonic()
    try:
        preparation_path = parent / "preparation_status.json"
        if preparation_path.exists():
            preparation = json.loads(preparation_path.read_text())
            if preparation["status"] != "SUCCEEDED":
                raise ValueError("previous preparation failed; no silent retry")
        else:
            minimized = minimize_complex(cohort["protein"], cohort["sdf"], name, parent / "minimization",
                **minimization_arguments(control))
            write_json(preparation_path, {"status": "SUCCEEDED", "complex": str(minimized.complex_pdb),
                "wall_seconds": time.monotonic()-started})
        _, _, result = score_job(cohort, name, Path(output), 2)
    except Exception as error:
        result = {"status": "FAILED", "error": f"{type(error).__name__}: {error}"}
        write_json(parent / "summary.json", result)
    return cohort["target"], name, result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", required=True, type=Path)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--control", choices=("input_heavy_fixed", "hydrogens_unrestrained"), default="input_heavy_fixed")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("workers must be positive")
    workspace, source, output = args.workspace.resolve(), args.source.resolve(), args.output.resolve()
    parent = json.loads((source / "frozen_protocol.json").read_text())
    cohorts = {target: deepcopy(parent["cohorts"][target]) for target in TARGETS}
    for cohort in cohorts.values():
        cohort["local_residues"] = []
        for path, fingerprint in cohort["input_sha256"].items():
            if digest(path) != fingerprint:
                raise ValueError("source fingerprint changed")
    files = [Path(__file__), Path(__file__).with_name("run_local_sidechain_transfer.py"),
        Path(__file__).parents[1] / "src/fmotpiefinal/core/minimization.py"]
    plan = {"schema": "matched-ligand-geometry-control-v2", "targets": list(TARGETS),
        "cohorts": cohorts, "code_sha256": {str(path.resolve()): digest(path) for path in files},
        "source_protocol_sha256": digest(source / "frozen_protocol.json"),
        "control": args.control, "minimization_arguments": minimization_arguments(args.control),
        "description": "Receptor fixed. input_heavy_fixed fixes ligand heavy atoms and retains H k10 restraints; hydrogens_unrestrained keeps mobile heavy atoms under k10 but frees H from positional restraints",
        "qm_model_changed": False, "fragment_membership_changed": False,
        "labels_used_for_control_selection": False, "production_default_changed": False}
    frozen = output / "frozen_protocol.json"
    if frozen.exists() and json.loads(frozen.read_text()) != plan:
        raise ValueError("refusing changed frozen pose control")
    write_json(frozen, plan)
    started = time.monotonic()
    results = {target: {} for target in TARGETS}
    baseline = workspace / "FMOtpiefinal_HIV1lab/runs/all_accelerated_robust_20261020"
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context("spawn")) as pool:
        futures = [pool.submit(execute, (cohort, name, str(output), args.control)) for cohort in cohorts.values() for name in cohort["names"]]
        for future in as_completed(futures):
            target, name, result = future.result()
            results[target][name] = result
            write_json(output / "progress.json", results)
            print(target, name, result["status"], result.get("tpie", result.get("error")), flush=True)
            if len(results[target]) == len(cohorts[target]["names"]):
                evaluate(workspace, baseline, results, output)
    evaluation = evaluate(workspace, baseline, results, output)
    failed = sum(row["status"] != "SUCCEEDED" for group in results.values() for row in group.values())
    write_json(output / "final_status.json", {"status": "COMPLETE" if not failed else "INCOMPLETE_WITH_FAILURES",
        "n_requested": sum(len(cohort["names"]) for cohort in cohorts.values()), "n_failed": failed,
        "wall_seconds": time.monotonic()-started, "evaluation_sha256": digest(output / "evaluation.json")})
    print(json.dumps(evaluation, indent=2), flush=True)


if __name__ == "__main__":
    main()
