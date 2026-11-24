"""Audit independent-repeat and time stability without treating frames as repeats."""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np

from report_endpoint_energy_comparison import HARTREE_TO_KCAL_MOL
from run_matched_openfe import digest, write_json


def paired_repeats(frames):
    grouped = defaultdict(dict)
    for frame in frames:
        key = (frame["leg"], frame["trajectory"])
        endpoint_frames = grouped[key].setdefault(frame["endpoint"], {})
        if frame["iteration"] in endpoint_frames:
            raise ValueError("duplicate endpoint frame")
        endpoint_frames[frame["iteration"]] = frame
    result = []
    for (leg, trajectory), endpoints in sorted(grouped.items()):
        if set(endpoints) != {0, 10} or set(endpoints[0]) != set(endpoints[10]):
            raise ValueError("both physical endpoints with matching times required")
        paired = []
        for iteration in sorted(endpoints[0]):
            first, second = endpoints[0][iteration], endpoints[10][iteration]
            mm = second["MM_total_kcal_mol"] - first["MM_total_kcal_mol"]
            qm = (second["DFTB_hartree"] - first["DFTB_hartree"]) * HARTREE_TO_KCAL_MOL
            paired.append({"iteration": iteration, "MM": mm, "DFTB": qm, "model_difference": qm - mm})
        means = {model: float(np.mean([frame[model] for frame in paired])) for model in ("MM", "DFTB", "model_difference")}
        early = [frame for frame in paired if frame["iteration"] <= 800]
        late = [frame for frame in paired if frame["iteration"] >= 1600]
        if not early or not late:
            raise ValueError("registered early and late windows required")
        result.append({"leg": leg, "trajectory": trajectory, "n_paired_frames": len(paired),
            "mean_B_minus_A_potential_kcal_mol": means,
            "late_minus_early_B_minus_A_potential_kcal_mol": {
                model: float(np.mean([frame[model] for frame in late]) - np.mean([frame[model] for frame in early]))
                for model in means}})
    return result


def stability_summary(frames, expected_repeats=3):
    repeats = paired_repeats(frames)
    legs = {leg: [row for row in repeats if row["leg"] == leg] for leg in ("bound", "solution")}
    if any(not rows for rows in legs.values()):
        raise ValueError("both legs required")
    reference = legs["solution"][0]["mean_B_minus_A_potential_kcal_mol"]
    result = {"repeat_counts": {leg: len(rows) for leg, rows in legs.items()}, "expected_repeats_per_leg": expected_repeats,
        "all_planned_repeats_present": all(len(rows) == expected_repeats for rows in legs.values()),
        "n_independent_bound_solution_combinations": None,
        "shared_reference_warning": "Bound repeats minus the same solution mean are not independent bound/solution combinations.",
        "convergence_proven": False, "models": {}}
    for model in reference:
        centered = {leg: np.array([row["mean_B_minus_A_potential_kcal_mol"][model] - reference[model] for row in rows])
                    for leg, rows in legs.items()}
        average = {leg: float(values.mean()) for leg, values in centered.items()}
        deviation = {leg: float(values.std(ddof=1)) if len(values) > 1 else None for leg, values in centered.items()}
        standard_error = (float(np.sqrt(sum(values.var(ddof=1) / len(values) for values in centered.values())))
                          if all(len(values) > 1 for values in centered.values()) else None)
        leave_out = []
        for leg, values in centered.items():
            if len(values) <= 1:
                continue
            for index, row in enumerate(legs[leg]):
                changed = dict(average)
                changed[leg] = float(np.delete(values, index).mean())
                leave_out.append({"removed_leg": leg, "removed_trajectory": row["trajectory"],
                    "double_difference_kcal_mol": changed["bound"] - changed["solution"]})
        result["models"][model] = {
            "internal_energy_double_difference_kcal_mol": average["bound"] - average["solution"],
            "leg_repeat_std_kcal_mol": deviation,
            "repeat_standard_error_diagnostic_kcal_mol": standard_error,
            "single_bound_repeat_minus_shared_solution_mean_kcal_mol": [float(value - average["solution"]) for value in centered["bound"]],
            "leave_one_repeat_out": leave_out,
            "centered_repeat_means_kcal_mol": {leg: values.tolist() for leg, values in centered.items()},
            "late_minus_early_by_repeat_kcal_mol": {leg: [row["late_minus_early_B_minus_A_potential_kcal_mol"][model] for row in rows]
                                                    for leg, rows in legs.items()}}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    comparison = json.loads(args.comparison.read_text())
    report = {"schema": "endpoint-independent-repeat-stability-v1", "targets": {},
        "comparison_sha256": digest(args.comparison), "driver_sha256": digest(__file__),
        "labels_used": False, "main_pipeline_changed": False,
        "limitations": ["Within-repeat frames and the two endpoint observations are not independent repeats.",
            "Repeat standard error is only a small-sample diagnostic, not total uncertainty or a confidence interval.",
            "One bound repeat cannot estimate bound repeat variability; no zero uncertainty substitution.",
            "Neither a binding free energy nor a valid MM-to-QM free-energy correction is estimated."]}
    for target, data in comparison["targets"].items():
        report["targets"][target] = stability_summary(data["matched_frames"])
    write_json(args.output, report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
