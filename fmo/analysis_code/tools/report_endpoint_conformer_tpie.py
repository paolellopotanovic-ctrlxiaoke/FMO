"""Summarize a complete fixed-receptor conformer control without fitting."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from run_matched_openfe import digest, write_json


def summarize_control(plan, rows):
    expected = {task["directory"] for task in plan["tasks"]}
    actual = [row["directory"] for row in rows]
    if len(set(actual)) != len(actual) or not set(actual) <= expected:
        raise ValueError("duplicate or unexpected task results")
    failed = [row["directory"] for row in rows if row["result"]["status"] != "SUCCEEDED"]
    complete = set(actual) == expected and not failed
    report = {"status": "COMPLETE" if complete else "INCOMPLETE", "n_expected": len(expected),
              "n_received": len(rows), "failed": failed, "missing": sorted(expected - set(actual)),
              "cohort_ranking_available": False, "binding_free_energy_claim_allowed": False,
              "scope": plan["scope"], "targets": {}}
    if not complete:
        return report
    for target in sorted({row["target"] for row in rows}):
        selected = [row for row in rows if row["target"] == target]
        ligands = {}
        for name in sorted({row["name"] for row in selected}):
            baseline = [row for row in selected if row["name"] == name and row["kind"] == "baseline_replay"]
            conformers = [row for row in selected if row["name"] == name and row["kind"] == "conformer"]
            if len(baseline) != 1 or len(conformers) != 9:
                raise ValueError("one baseline and nine conformers per ligand required")
            repeats = {}
            for trajectory in sorted({row["trajectory"] for row in conformers}):
                frames = [row for row in conformers if row["trajectory"] == trajectory]
                if sorted(row["iteration"] for row in frames) != [400, 1200, 2000]:
                    raise ValueError("fixed three-frame schedule required per repeat")
                repeats[trajectory] = float(np.mean([row["result"]["tpie_kcal_mol"] for row in frames]))
            if len(repeats) != 3:
                raise ValueError("three independent repeats required")
            means = list(repeats.values())
            baseline_tpie = baseline[0]["result"]["tpie_kcal_mol"]
            ligands[name] = {"baseline_tpie_kcal_mol": baseline_tpie,
                "conformer_mean_tpie_kcal_mol": float(np.mean(means)),
                "repeat_std_kcal_mol": float(np.std(means, ddof=1)),
                "repeat_means_kcal_mol": repeats,
                "mean_shift_from_baseline_kcal_mol": float(np.mean(means) - baseline_tpie),
                "max_internal_rmsd_angstrom": max(row["geometry"]["ligand_internal_rmsd_angstrom"] for row in conformers)}
        report["targets"][target] = {"ligands": ligands}
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True, type=Path)
    args = parser.parse_args()
    root = args.run.resolve()
    plan = json.loads((root / "frozen_protocol.json").read_text())
    rows = json.loads((root / "progress.json").read_text())
    report = summarize_control(plan, rows)
    report["report_driver_sha256"] = digest(__file__)
    write_json(root / "evaluation.json", report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
