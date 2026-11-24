"""Evaluate a registered QM-gated input-pose control, retaining every cohort member."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from run_cdk8_solution_cycle import digest, metrics, write_json
from report_bace_extensions import audit, paired_rank_summary
from report_local_sidechain_transfer import load_truth


def select_candidate(candidate, baseline):
    if candidate["status"] == "SUCCEEDED":
        return candidate["tpie"], "new_input_heavy_fixed"
    if candidate.get("error") == "ValueError: QM status: FAILED_UNPHYSICAL_PIEDA":
        if baseline.get("status") != "SUCCEEDED":
            raise ValueError("fallback requires independently successful baseline")
        return baseline["tpie"], "baseline_retained_after_unphysical_input_pose"
    raise ValueError("unregistered failure must not silently fall back")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--control", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source, control, output = args.source.resolve(), args.control.resolve(), args.output.resolve()
    frozen = {"schema": "qm-gated-input-heavy-retention-v1", "source": str(source),
        "source_sha256": digest(source), "control": str(control),
        "control_protocol_sha256": digest(control / "frozen_protocol.json"),
        "driver_sha256": digest(__file__), "labels_used_for_selection": False,
        "selection": "Use new fixed-input-heavy tPIE if strict QM succeeds; retain independently verified current baseline only for FAILED_UNPHYSICAL_PIEDA; all other errors fail closed",
        "timing": "Registered after structural rejection of CAT-13i, before evaluating this guarded complete-cohort score",
        "production_default_changed": False}
    plan = output / "frozen_protocol.json"
    if plan.exists() and json.loads(plan.read_text()) != frozen:
        raise ValueError("guarded protocol changed")
    write_json(plan, frozen)
    cohort = json.loads(source.read_text())["cohorts"]["BACE"]
    evaluation = json.loads((control / "evaluation.json").read_text())
    names = sorted(cohort["names"])
    if set(evaluation["rows"]) != set(names):
        raise ValueError("missing cohort members")
    for path, expected in cohort["input_sha256"].items():
        if digest(path) != expected:
            raise ValueError("baseline evidence changed")
    rows = {}
    for name in names:
        candidate = evaluation["rows"][name]
        baseline = dict(cohort["baseline"][name], status=cohort["runtime"][name]["status"])
        score, selection = select_candidate(candidate, baseline)
        if selection == "new_input_heavy_fixed":
            directory = control / "BACE" / name
            if digest(directory / "input.inp") != candidate["input_sha256"] or digest(directory / "run/gamess.log") != candidate["log_sha256"]:
                raise ValueError("new candidate evidence changed")
        else:
            failure = json.loads((control / "BACE" / name / "run/result.json").read_text())
            if failure["status"] != "FAILED_UNPHYSICAL_PIEDA":
                raise ValueError("fallback trigger not substantiated")
        rows[name] = {"tpie": score, "selection": selection, "candidate": candidate}
    truth = load_truth(args.workspace.resolve(), "BACE")
    if set(truth) != set(names):
        raise ValueError("evaluation identity mismatch")
    observed = [truth[name] for name in names]
    baseline = [-cohort["baseline"][name]["tpie"] for name in names]
    predicted = [-rows[name]["tpie"] for name in names]
    result = {"status": "COMPLETE_GUARDED_DEVELOPMENT_CONTROL", "n": len(names), "rows": rows,
        "n_new_geometry": sum(row["selection"] == "new_input_heavy_fixed" for row in rows.values()),
        "n_baseline_retained": sum(row["selection"] != "new_input_heavy_fixed" for row in rows.values()),
        "metrics": metrics(observed, predicted), "baseline_metrics": metrics(observed, baseline),
        "paired_comparison": paired_rank_summary(observed, baseline, predicted),
        "warning": "Not the pure fixed-heavy arm, which failed 1/36; fallback explicitly preserves the missing case. Not a new blind evaluation.",
        "production_default_changed": False}
    write_json(output / "evaluation.json", result)


if __name__ == "__main__":
    main()
