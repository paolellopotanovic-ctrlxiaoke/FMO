"""Complete the fixed sidechain-flexibility by H-restraint factorial control."""
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


def execute(payload):
    cohort, name, output = payload
    from fmotpiefinal.core.minimization import minimize_complex
    parent = Path(output) / cohort["target"] / directory_key(name)
    try:
        preparation_path = parent / "preparation_status.json"
        if preparation_path.exists():
            if json.loads(preparation_path.read_text())["status"] != "SUCCEEDED":
                raise ValueError("previous preparation failed; no silent retry")
        else:
            minimized = minimize_complex(cohort["protein"], cohort["sdf"], name, parent / "minimization",
                restraint_k_kj_mol_a2=10.0, move_protein=False,
                local_flexible_residues=tuple(map(tuple, cohort["local_residues"])),
                receptor_restraint_k_kj_mol_a2=2.5, restrain_ligand_hydrogens=False)
            write_json(preparation_path, {"status": "SUCCEEDED", "complex": str(minimized.complex_pdb)})
        _, _, result = score_job(cohort, name, Path(output), 2)
    except Exception as error:
        result = {"status": "FAILED", "error": f"{type(error).__name__}: {error}"}
        write_json(parent / "summary.json", result)
    return cohort["target"], name, result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=10)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("workers must be positive")
    workspace, source, output = args.workspace.resolve(), args.source.resolve(), args.output.resolve()
    original = json.loads((source / "frozen_protocol.json").read_text())
    cohorts = {target: deepcopy(original["cohorts"][target]) for target in TARGETS}
    for cohort in cohorts.values():
        for path, fingerprint in cohort["input_sha256"].items():
            if digest(path) != fingerprint:
                raise ValueError("source input fingerprint changed")
    files = [Path(__file__), Path(__file__).with_name("run_local_sidechain_transfer.py"),
        Path(__file__).parents[1] / "src/fmotpiefinal/core/minimization.py"]
    plan = {"schema": "sidechain-hydrogen-factorial-cell-v1", "cohorts": cohorts,
        "control": "sidechain_hydrogens_unrestrained", "code_sha256": {str(path.resolve()): digest(path) for path in files},
        "source_protocol_sha256": digest(source / "frozen_protocol.json"),
        "change_from_sidechain_B": "ligand hydrogen positional restraints removed; same frozen sidechain graph and receptor k2.5, ligand heavy k10",
        "labels_used_for_selection": False, "production_default_changed": False}
    frozen = output / "frozen_protocol.json"
    if frozen.exists() and json.loads(frozen.read_text()) != plan:
        raise ValueError("refusing changed frozen factorial control")
    write_json(frozen, plan)
    started = time.monotonic()
    rows = {target: {} for target in TARGETS}
    baseline = workspace / "FMOtpiefinal_HIV1lab/runs/all_accelerated_robust_20261020"
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context("spawn")) as pool:
        futures = [pool.submit(execute, (cohort, name, str(output))) for cohort in cohorts.values() for name in cohort["names"]]
        for future in as_completed(futures):
            target, name, result = future.result()
            rows[target][name] = result
            write_json(output / "progress.json", rows)
            print(target, name, result["status"], result.get("tpie", result.get("error")), flush=True)
    evaluate(workspace, baseline, rows, output)
    failed = sum(row["status"] != "SUCCEEDED" for group in rows.values() for row in group.values())
    write_json(output / "final_status.json", {"status": "COMPLETE" if not failed else "INCOMPLETE_WITH_FAILURES",
        "n_requested": sum(len(cohort["names"]) for cohort in cohorts.values()), "n_failed": failed,
        "wall_seconds": time.monotonic()-started, "evaluation_sha256": digest(output / "evaluation.json")})


if __name__ == "__main__":
    main()
