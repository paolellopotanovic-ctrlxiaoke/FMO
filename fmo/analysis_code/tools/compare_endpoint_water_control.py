"""Validate water-only paired inputs and compare all predefined repeats."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from analyze_cdk8_geometry_evidence import fragments, xyz
from report_dense_endpoint_complex_tpie import is_water_fragment
from run_cdk8_solution_cycle import group
from run_dense_endpoint_complex_tpie import task_identity
from run_matched_openfe import digest, write_json


def validate_water_only_change(original, updated):
    first, second = group(original, "FMOXYZ"), group(updated, "FMOXYZ")
    if original[:first.start()] + original[first.end():] != updated[:second.start()] + updated[second.end():]:
        raise ValueError("non-coordinate controls changed")
    labels, membership = fragments(original)
    symbols, coordinates = xyz(original)
    other_symbols, other_coordinates = xyz(updated)
    if symbols != other_symbols:
        raise ValueError("elements or atom count changed")
    water_atoms = {index - 1 for label, members in zip(labels, membership, strict=True)
                   if is_water_fragment(label) for index in members}
    fixed = sorted(set(range(len(symbols))) - water_atoms)
    if not np.allclose(coordinates[fixed], other_coordinates[fixed], atol=1e-9, rtol=0):
        raise ValueError("protein, ligand or caps changed")
    return len(water_atoms)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("original", "control", "coverage"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    original = json.loads((args.original / "evaluation.json").read_text())
    control = json.loads((args.control / "evaluation.json").read_text())
    report = {"status": "INCOMPLETE", "labels_used": False, "main_pipeline_changed": False,
        "cohort_ranking_available": False, "binding_free_energy_claim_allowed": False,
        "targets": {}, "driver_sha256": digest(__file__),
        "warning": "Fixed-count nearest-water replacement is a diagnostic control, not a validated default or solvent free energy."}
    if original["status"] == control["status"] == "COMPLETE":
        plan = json.loads((args.control / "frozen_protocol.json").read_text())
        coverage = json.loads(args.coverage.read_text())
        if coverage["status"] != "COMPLETE":
            raise ValueError("complete independent water coverage audit required")
        original_plan = json.loads((args.original / "frozen_protocol.json").read_text())
        original_lookup = {task_identity(task): task for task in original_plan["tasks"]}
        checks = []
        for task in plan["tasks"]:
            reference = original_lookup[task_identity(task)]
            path = Path(task["directory"]) / "input.inp"
            reference_path = Path(reference["directory"]) / "input.inp"
            if digest(path) != task["input_sha256"] or digest(reference_path) != reference["input_sha256"]:
                raise ValueError("matched input changed")
            changed_count = validate_water_only_change(reference_path.read_text(), path.read_text())
            mappings = json.loads(Path(task["directory"], "water_mapping.json").read_text())
            checks.append({"directory": task["directory"], "water_atom_slots": changed_count,
                           "n_selected_nearby_water_oxygens": sum(water["minimum_polar_distance_angstrom"] < 3.5 for water in mappings)})
        for target, data in control["targets"].items():
            previous = original["targets"][target]
            if set(data["repeats"]) != set(previous["repeats"]):
                raise ValueError("matched repeat identities differ")
            changes = []
            for trajectory in sorted(data["repeats"]):
                new = data["repeats"][trajectory]["B_minus_A_kcal_mol"]["tPIE"]
                old = previous["repeats"][trajectory]["B_minus_A_kcal_mol"]["tPIE"]
                changes.append({"trajectory": trajectory, "original_B_minus_A": old, "control_B_minus_A": new,
                                "control_minus_original": new - old})
            coverage_lookup = {(row["name"], row["trajectory"], row["iteration"]): row
                               for row in coverage["targets"][target]["frames"]}
            contact_checks = []
            for task, check in zip(plan["tasks"], checks, strict=True):
                if task["target"] != target:
                    continue
                frame = coverage_lookup[task["name"], task["trajectory"], task["iteration"]]
                missing = frame["n_nearby_water_oxygens"] - check["n_selected_nearby_water_oxygens"]
                if missing < 0:
                    raise ValueError("selected nearby waters exceed independently measured total")
                contact_checks.append({"name": task["name"], "iteration": task["iteration"], "trajectory": task["trajectory"],
                                       "omitted_before": frame["n_nearby_omitted_water_oxygens"], "omitted_after": missing})
            report["targets"][target] = {"paired_repeats": changes,
                "original_mean_B_minus_A": previous["mean_B_minus_A_kcal_mol"]["tPIE"],
                "control_mean_B_minus_A": data["mean_B_minus_A_kcal_mol"]["tPIE"],
                "original_repeat_std": previous["repeat_std_B_minus_A_kcal_mol"]["tPIE"],
                "control_repeat_std": data["repeat_std_B_minus_A_kcal_mol"]["tPIE"],
                "mean_control_minus_original": float(np.mean([row["control_minus_original"] for row in changes])),
                "water_coverage_by_frame": contact_checks}
        report["input_checks"] = checks
        report["status"] = "COMPLETE"
    report["source_sha256"] = {str(path): digest(path) for path in
                               (args.original / "evaluation.json", args.control / "evaluation.json", args.coverage)}
    write_json(args.control / "paired_water_comparison.json", report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
