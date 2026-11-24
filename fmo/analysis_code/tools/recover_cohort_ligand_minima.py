"""Apply a label-blind RFO optimization rescue without weakening convergence."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from run_cdk8_solution_cycle import digest, group, write_json
from run_cohort_solution_cycle import evaluate, execute
from recover_cdk8_solution_cycle import recovery_deck


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    args = parser.parse_args()
    root, workspace = args.run.resolve(), args.workspace.resolve()
    source = root / "receptor_data_title_repair"
    protocol = source / "recovered_protocol.json"
    frozen = json.loads(protocol.read_text())
    rows = json.loads((source / "evaluation.json").read_text())["rows"]
    selected = {key: frozen["tasks"][key] for key, row in rows.items()
                if row["status"] != "SUCCEEDED" and frozen["tasks"][key]["state"] == "free_ligand_minimum"
                and "isolated ligand geometry optimization did not converge" in row.get("error", "")}
    recovery = root / "ligand_minimum_rfo_recovery"
    tasks = {}
    for key, task in selected.items():
        old = Path(task["input"])
        if digest(old) != task["input_sha256"]:
            raise ValueError("original optimizer input changed")
        deck = old.read_text()
        corrected = recovery_deck(deck, True, "charge_diis_rfo")
        if corrected == deck or group(corrected, "FMOXYZ")[1] != group(deck, "FMOXYZ")[1]:
            raise ValueError("unrecognized optimization budget or changed starting coordinates")
        parent = recovery / key
        parent.mkdir(parents=True, exist_ok=True)
        link = parent.parent / "params"
        if not link.exists():
            link.symlink_to(Path(frozen["parameters"]), target_is_directory=True)
        path = parent / "input.inp"
        if path.exists() and path.read_text() != corrected:
            raise ValueError("refusing changed numerical recovery")
        path.write_text(corrected)
        tasks[key] = dict(task, input=str(path), input_sha256=digest(path))
        frozen["tasks"][key] = tasks[key]
    plan = {"schema": "strict-unchanged-opt-tolerance-rfo-rescue-v1", "source_sha256": digest(protocol),
        "driver_sha256": digest(__file__), "tasks": tasks, "labels_used": False,
        "profile": "All and only failed 200-step ligand optimizations: RFO, NSTEP=500, TRMIN=0.001; unchanged OPTTOL=0.000100 and all SCF thresholds, starting geometry, charge, energy and solvent model"}
    path = recovery / "frozen_recovery.json"
    if path.exists() and json.loads(path.read_text()) != plan:
        raise ValueError("changed numerical recovery plan")
    write_json(path, plan)
    for key, task in tasks.items():
        rows[key] = execute((task, str(workspace)))[1]
        print(key, rows[key]["status"], rows[key].get("error", ""), flush=True)
    report = evaluate(rows, frozen, workspace)
    report["numerical_recovery"] = plan
    write_json(recovery / "evaluation.json", report)
    write_json(recovery / "recovered_protocol.json", frozen)


if __name__ == "__main__":
    main()
