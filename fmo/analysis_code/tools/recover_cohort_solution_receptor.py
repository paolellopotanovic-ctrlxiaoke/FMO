"""Repair only a failed receptor DATA title; preserve all other frozen work."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from run_cdk8_solution_cycle import digest, group, write_json
from run_cohort_solution_cycle import evaluate, execute


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--evaluate", action="store_true")
    args = parser.parse_args()
    root, workspace = args.run.resolve(), args.workspace.resolve()
    source = root / "frozen_protocol.json"
    original = json.loads(source.read_text())
    recovered = json.loads(source.read_text())
    repaired = {}
    for key, task in original["tasks"].items():
        if task["state"] != "receptor":
            continue
        path = Path(task["input"])
        if digest(path) != task["input_sha256"]:
            raise ValueError("original input changed")
        deck = path.read_text()
        if "$DATA\n\nMatched common receptor\nC1" not in deck:
            raise ValueError("recovery does not match the one registered formatting error")
        corrected = deck.replace("$DATA\n\nMatched common receptor\nC1", "$DATA\nMatched common receptor\nC1", 1)
        if group(corrected, "FMOXYZ")[1] != group(deck, "FMOXYZ")[1]:
            raise ValueError("recovery changed coordinates")
        parent = root / "receptor_data_title_repair" / key
        parent.mkdir(parents=True, exist_ok=True)
        link = parent.parent / "params"
        if not link.exists():
            link.symlink_to(Path(original["parameters"]), target_is_directory=True)
        target = parent / "input.inp"
        if target.exists() and target.read_text() != corrected:
            raise ValueError("refusing changed repair")
        target.write_text(corrected)
        replacement = dict(task, input=str(target), input_sha256=digest(target))
        recovered["tasks"][key] = replacement
        repaired[key] = replacement
    plan = {"schema": "data-title-blank-line-only-recovery-v1", "source_sha256": digest(source),
        "driver_sha256": digest(__file__), "repaired_tasks": repaired,
        "change": "Remove one blank line after DATA; no atom/model/charge/numerical change",
        "labels_used": False, "old_failed_logs_preserved": True}
    plan_path = root / "receptor_data_title_repair/frozen_recovery.json"
    if plan_path.exists() and json.loads(plan_path.read_text()) != plan:
        raise ValueError("recovery plan changed")
    write_json(plan_path, plan)
    repaired_rows = {key: execute((task, str(workspace)))[1] for key, task in repaired.items()}
    write_json(root / "receptor_data_title_repair/results.json", repaired_rows)
    if args.evaluate:
        progress = json.loads((root / "progress.json").read_text())
        if progress["n_finished"] != progress["n_expected"]:
            raise ValueError("original jobs still running; do not evaluate a partial cohort")
        rows = dict(progress["rows"], **repaired_rows)
        report = evaluate(rows, recovered, workspace)
        report["recovery"] = plan
        write_json(root / "receptor_data_title_repair/evaluation.json", report)
        write_json(root / "receptor_data_title_repair/recovered_protocol.json", recovered)


if __name__ == "__main__":
    main()
