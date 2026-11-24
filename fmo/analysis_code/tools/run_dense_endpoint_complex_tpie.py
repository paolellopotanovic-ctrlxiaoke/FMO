"""Run all predefined saved bound frames, with exact-input validated reuse."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import time

import run_endpoint_complex_tpie as complex_driver
from run_matched_openfe import digest, write_json


def task_identity(task):
    return tuple(task[key] for key in ("target", "name", "trajectory", "iteration", "replica", "endpoint"))


def reuse_result(task, cached):
    if task_identity(task) != task_identity(cached):
        raise ValueError("cache identity mismatch")
    for key in ("input_sha256", "source_sha256", "runtime", "parameters", "ligand_fragment"):
        if task[key] != cached[key]:
            raise ValueError(f"cache chemistry mismatch: {key}")
    for entry in (task, cached):
        if digest(Path(entry["directory"]) / "input.inp") != entry["input_sha256"]:
            raise ValueError("cache input digest mismatch")
        if digest(entry["source"]) != entry["source_sha256"]:
            raise ValueError("source digest mismatch")
    summary_path = Path(cached["directory"]) / "summary.json"
    result = json.loads(summary_path.read_text())["result"]
    if result != cached["result"] or result["status"] != "SUCCEEDED":
        raise ValueError("cache must be verified successful and unchanged")
    row = {**task, "result": result, "reused_from": str(summary_path),
           "reused_summary_sha256": digest(summary_path)}
    write_json(Path(task["directory"]) / "summary.json", row)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("workspace", "run", "output", "reuse"):
        parser.add_argument(f"--{name}", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--timeout", type=int, default=1800)
    args = parser.parse_args()
    if min(args.workers, args.timeout) < 1:
        parser.error("positive budgets required")
    started = time.monotonic()
    output = args.output.resolve()
    plan = complex_driver.freeze(args, iterations=tuple(range(400, 2001, 40)), driver_path=__file__)
    reuse_path = args.reuse.resolve() / "progress.json"
    cached = json.loads(reuse_path.read_text())
    identities = [task_identity(row) for row in cached]
    if len(set(identities)) != len(identities):
        raise ValueError("duplicate cached task")
    lookup = dict(zip(identities, cached, strict=True))
    plan["source_sha256"][str(reuse_path)] = digest(reuse_path)
    plan["source_sha256"][str(Path(complex_driver.__file__).resolve())] = digest(complex_driver.__file__)
    plan["workers"] = args.workers
    plan["schedule_reason"] = "All saved 0.1 ns frames from 1 to 5 ns; no label or score selection."
    write_json(output / "frozen_protocol.json", plan)
    rows, pending = [], []
    for task in plan["tasks"]:
        candidate = lookup.get(task_identity(task))
        if candidate is not None and candidate["result"]["status"] == "SUCCEEDED":
            rows.append(reuse_result(task, candidate))
        else:
            pending.append(task)
    reused = len(rows)
    write_json(output / "progress.json", rows)
    write_json(output / "status.json", {"status": "RUNNING", "n_requested": len(plan["tasks"]),
        "n_reused": reused, "n_succeeded": reused, "workers": args.workers})
    print("FROZEN", len(plan["tasks"]), "REUSED", reused, "NEW", len(pending), flush=True)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(complex_driver.execute, task, args.timeout): task for task in pending}
        for future in as_completed(futures):
            task = futures[future]
            try:
                row = future.result()
            except Exception as error:
                row = {**task, "result": {"status": "FAILED", "error": f"{type(error).__name__}: {error}"}}
                write_json(Path(task["directory"]) / "summary.json", row)
            rows.append(row)
            write_json(output / "progress.json", rows)
            print(len(rows), len(plan["tasks"]), task["target"], task["name"], row["result"]["status"], flush=True)
    succeeded = sum(row["result"]["status"] == "SUCCEEDED" for row in rows)
    write_json(output / "status.json", {"status": "COMPLETE" if succeeded == len(plan["tasks"]) else "INCOMPLETE_WITH_FAILURES",
        "n_requested": len(plan["tasks"]), "n_succeeded": succeeded, "n_reused": reused,
        "n_new_executed": len(pending), "wall_seconds_including_preflight": time.monotonic() - started,
        "cohort_ranking_available": False, "binding_free_energy_claim_allowed": False,
        "rows": rows})


if __name__ == "__main__":
    main()
