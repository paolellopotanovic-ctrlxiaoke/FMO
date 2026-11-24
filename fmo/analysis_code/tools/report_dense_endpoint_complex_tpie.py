"""Assess fixed-window and independent-repeat sensitivity without fitting scores."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import subprocess
import sys

import numpy as np

from run_matched_openfe import digest, write_json


def is_water_fragment(name):
    return re.fullmatch(r".?(?:HOH|WAT|TIP)[+-]?\d+[A-Za-z]?", name.upper()) is not None


def repeat_diagnostics(differences):
    values = np.asarray(differences, float)
    if values.shape != (3, 41) or not np.isfinite(values).all():
        raise ValueError("three complete repeats of 41 fixed paired frames required")
    means = values.mean(axis=1)
    windows = {}
    for name, selected in (("early_1_to_2p9_ns", slice(0, 20)), ("late_3_to_5_ns", slice(20, 41))):
        per_repeat = values[:, selected].mean(axis=1)
        windows[name] = {"paired_repeat_means": per_repeat.tolist(), "mean": float(per_repeat.mean()),
                         "repeat_std": float(per_repeat.std(ddof=1))}
    prefixes = {}
    for count in (11, 21, 31, 41):
        per_repeat = values[:, :count].mean(axis=1)
        prefixes[str(count)] = {"ending_ns": 1 + (count - 1) * .1,
                               "paired_repeat_means": per_repeat.tolist(),
                               "mean": float(per_repeat.mean()), "repeat_std": float(per_repeat.std(ddof=1))}
    return {"paired_repeat_means": means.tolist(), "mean": float(means.mean()),
        "repeat_std": float(means.std(ddof=1)),
        "all_three_repeats_same_nonzero_sign": bool(np.all(means > 0) or np.all(means < 0)),
        "leave_one_repeat_out_means": [float(np.delete(means, index).mean()) for index in range(3)],
        "fixed_windows": windows, "prefix_sensitivity": prefixes,
        "late_minus_early_by_repeat": (values[:, 20:].mean(axis=1) - values[:, :20].mean(axis=1)).tolist(),
        "convergence_proven": False,
        "warning": "Saved frames are correlated; 123 paired observations are not 123 independent repeats."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("run", "baseline", "openfe"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    root = args.run.resolve()
    with (root / "report.stdout").open("w") as handle:
        subprocess.run([sys.executable, str(Path(__file__).with_name("report_endpoint_complex_tpie.py")),
                        "--run", str(root), "--baseline", str(args.baseline), "--openfe", str(args.openfe)],
                       stdout=handle, check=True)
    report = json.loads((root / "evaluation.json").read_text())
    plan = json.loads((root / "frozen_protocol.json").read_text())
    if plan["fixed_iterations"] != list(range(400, 2001, 40)):
        raise ValueError("fixed dense schedule required")
    rows = json.loads((root / "progress.json").read_text())
    if report["status"] == "COMPLETE":
        for target, data in report["targets"].items():
            selected = [row for row in rows if row["target"] == target]
            matrices = {}
            water_fragments = {row["fragment_id"] for row in data["receptor_fragment_B_minus_A_kcal_mol"]
                               if is_water_fragment(row["fragment_name"])}
            for quantity in ("tPIE", "retained_water_pair_sum"):
                differences = []
                for trajectory in sorted({row["trajectory"] for row in selected}):
                    frame_values = {name: [] for name in (data["first"], data["second"])}
                    for name in frame_values:
                        frames = sorted((row for row in selected if row["trajectory"] == trajectory and row["name"] == name),
                                        key=lambda row: row["iteration"])
                        for row in frames:
                            value = row["result"]["tpie_kcal_mol"]
                            if quantity == "retained_water_pair_sum":
                                ligand = row["ligand_fragment"]
                                value = sum(pair["value_kcal_mol"] for pair in row["result"]["pairs"]
                                            if ligand in (pair["fragment_i"], pair["fragment_j"])
                                            and (pair["fragment_j"] if pair["fragment_i"] == ligand else pair["fragment_i"]) in water_fragments)
                            frame_values[name].append(value)
                    differences.append(np.asarray(frame_values[data["second"]]) - frame_values[data["first"]])
                matrices[quantity] = repeat_diagnostics(differences)
            data["dense_stability"] = matrices
            data["water_decomposition_warning"] = "Diagnostic subset of the same tPIE, not a solvent free energy, new correction or default water deletion."
    report["dense_report_driver_sha256"] = digest(__file__)
    write_json(root / "dense_evaluation.json", report)
    compact = {"status": report["status"], "n_expected": report["n_expected"],
               "n_received": report["n_received"], "failed": report["failed"],
               "targets": {target: data["dense_stability"] for target, data in report["targets"].items()}}
    write_json(root / "dense_summary.json", compact)
    print(json.dumps(compact, indent=2))


if __name__ == "__main__":
    main()
