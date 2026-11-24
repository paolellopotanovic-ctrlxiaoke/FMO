"""Compare matched single-frame isolated-ligand potentials, never free energies."""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np

from run_endpoint_dftb_internal import FIXED_ITERATIONS
from run_matched_openfe import digest, write_json


HARTREE_TO_KCAL_MOL = 627.5094740631


def grouped_summary(frames, iterations=FIXED_ITERATIONS):
    grouped = defaultdict(list)
    for frame in frames:
        grouped[(frame["leg"], frame["endpoint"], frame["trajectory"])].append(frame)
    summaries = []
    for (leg, endpoint, trajectory), rows in sorted(grouped.items()):
        if tuple(sorted(row["iteration"] for row in rows)) != tuple(iterations):
            raise ValueError("matched fixed-time frame set incomplete")
        summaries.append({"leg": leg, "endpoint": endpoint, "trajectory": trajectory,
            "n_frames": len(rows), "MM_mean_kcal_mol": float(np.mean([row["MM_total_kcal_mol"] for row in rows])),
            "DFTB_mean_hartree": float(np.mean([row["DFTB_hartree"] for row in rows]))})
    within = {}
    for endpoint in (0, 10):
        legs = {leg: [row for row in summaries if row["leg"] == leg and row["endpoint"] == endpoint]
                for leg in ("bound", "solution")}
        if any(not rows for rows in legs.values()):
            raise ValueError("matched bound and solution data required")
        within[str(endpoint)] = {}
        for model, key, conversion in (("MM", "MM_mean_kcal_mol", 1), ("DFTB", "DFTB_mean_hartree", HARTREE_TO_KCAL_MOL)):
            average = {leg: float(np.mean([row[key] for row in rows])) for leg, rows in legs.items()}
            within[str(endpoint)][model] = (average["bound"] - average["solution"]) * conversion
    double = {model: within["10"][model] - within["0"][model] for model in ("MM", "DFTB")}
    return {"repeat_means": summaries, "bound_minus_solution_mean_internal_energy_kcal_mol": within,
            "B_minus_A_internal_energy_double_difference_kcal_mol": double,
            "DFTB_minus_MM_double_difference_kcal_mol": double["DFTB"] - double["MM"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mm", type=Path, required=True)
    parser.add_argument("--qm", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    mm = json.loads(args.mm.read_text())
    qm = json.loads(args.qm.read_text())
    plan_path = args.qm.parent / "frozen_protocol.json"
    plan = json.loads(plan_path.read_text())
    if (qm["status"] != "COMPLETE" or qm["n_requested"] != qm["n_succeeded"]
            or qm["n_requested"] != len(qm["results"]) or mm["status"] != "COMPLETE"):
        raise ValueError("complete MM and QM diagnostics required; no dropping failed frames")
    if plan["diagnostics_sha256"] != mm["diagnostics_sha256"]:
        raise ValueError("MM and QM frame provenance differs")
    tasks = {(task["target"], task["directory"]): task for task in plan["tasks"]}
    result_keys = [(row["target"], row["directory"]) for row in qm["results"]]
    if len(set(result_keys)) != len(result_keys) or set(result_keys) != set(tasks):
        raise ValueError("QM results must cover each frozen task exactly once")
    lookup = {}
    for target, data in mm["targets"].items():
        for leg, rows in data["legs"].items():
            for repeat in rows["repeats"]:
                for frame in repeat["frames"]:
                    lookup[(target, leg, repeat["endpoint"], repeat["trajectory"], frame["iteration"])] = frame
    grouped = defaultdict(list)
    for row in qm["results"]:
        task = tasks[(row["target"], row["directory"])]
        key = (row["target"], row["leg"], row["endpoint"], task["trajectory"], row["iteration"])
        frame = lookup[key]
        if frame["replica"] != row["replica"]:
            raise ValueError("MM/QM replica mismatch")
        grouped[row["target"]].append({"leg": row["leg"], "endpoint": row["endpoint"],
            "trajectory": task["trajectory"], "iteration": row["iteration"],
            "MM_total_kcal_mol": frame["energies_kcal_mol"]["total"], "DFTB_hartree": row["energy_hartree"]})
    report = {"schema": "matched-endpoint-MM-DFTB-internal-v1", "targets": {},
              "sources": {str(path): digest(path) for path in (args.mm, args.qm, plan_path)},
              "labels_used": False, "main_pipeline_changed": False,
              "interpretation": "Comparisons use identical fixed frames, molecules and gas-phase isolation, but different MM/DFTB Hamiltonians.",
              "warning": "Mean isolated potential differences are NOT free-energy corrections; do not subtract them from FEP or add to tPIE.",
              "QM_wall_seconds": qm["wall_seconds"], "QM_singlepoint_wall_seconds_sum": qm["sum_singlepoint_wall_seconds"]}
    for target, rows in grouped.items():
        report["targets"][target] = {**grouped_summary(rows, plan.get("fixed_iterations", FIXED_ITERATIONS)), "matched_frames": rows,
            "all_41_frame_MM_double_difference_kcal_mol": mm["targets"][target]["B_minus_A_internal_energy_double_difference_kcal_mol"]}
        iterations = plan.get("fixed_iterations", list(FIXED_ITERATIONS))
        if len(iterations) == 41:
            subsets = {"integer_ns_only": list(FIXED_ITERATIONS),
                       "early_1_to_2ns": [value for value in iterations if value <= 800],
                       "late_4_to_5ns": [value for value in iterations if value >= 1600],
                       "alternating_grid_0": iterations[::2], "alternating_grid_1": iterations[1::2]}
            report["targets"][target]["fixed_subwindow_sensitivity"] = {
                name: grouped_summary([row for row in rows if row["iteration"] in selected], selected)
                for name, selected in subsets.items()}
        print(target, json.dumps(report["targets"][target]["B_minus_A_internal_energy_double_difference_kcal_mol"]))
    write_json(args.output, report)


if __name__ == "__main__":
    main()
