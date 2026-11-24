"""Read completed/running pilot diagnostics without claiming cohort accuracy."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml

from report_local_sidechain_transfer import load_truth
from run_matched_openfe import digest, write_json


def leg_snapshot(path, thermal_energy, solution_estimate):
    rows = yaml.safe_load(path.read_text())
    if not isinstance(rows, list) or not rows:
        raise ValueError("nonempty diagnostic history required")
    last = rows[-1]
    analysis = last["mbar_analysis"]
    return {"source": str(path), "sha256": digest(path), "iteration": last["iteration"],
        "production_complete": last["percent_complete"] == 100,
        "bound_estimate_kcal_mol": analysis["free_energy_in_kT"] * thermal_energy,
        "binding_B_minus_A_kcal_mol": analysis["free_energy_in_kT"] * thermal_energy - solution_estimate,
        "bound_internal_mbar_se_kcal_mol": analysis["standard_error_in_kT"] * thermal_energy,
        "reported_uncorrelated_samples": analysis["number_of_uncorrelated_samples"],
        "reported_equilibrium_discard_iterations": analysis["n_equilibrium_iterations"],
        "history": [{"iteration": row["iteration"],
            "binding_B_minus_A_kcal_mol": row["mbar_analysis"]["free_energy_in_kT"] * thermal_energy - solution_estimate}
            for row in rows if row["iteration"] in (400, 800, 1200, 1600, 2000)],
        "full_uncertainty_estimated": False,
        "warning": "Single-history MBAR error is not independent-repeat, chemical-state or model uncertainty."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads((args.run / "frozen_protocol.json").read_text())
    report = {"schema": "partial-openfe-pilot-diagnostics-v1", "targets": {},
        "temperature_kelvin": 298.15, "full_cohort_metrics_available": False,
        "scope": "first predefined edge per target only; not validation of a QM correction or cohort ranking",
        "network_expansion_approved": False}
    for target in ("CDK8", "BACE"):
        root = args.run / "production" / target / "edge_000"
        solution_path = root / "solution/status.json"
        solution = json.loads(solution_path.read_text())
        if solution["status"] != "SUCCEEDED":
            raise ValueError("completed solution leg required for partial difference")
        first, second = plan["targets"][target]["edges"][0]
        truth = load_truth(args.workspace.resolve(), target)
        snapshots = [leg_snapshot(path, 0.00198720425864083 * 298.15, solution["estimate_kcal_mol"])
                     for path in sorted((root / "bound/shared").glob("*/simulation_real_time_analysis.yaml"))]
        report["targets"][target] = {"first": first, "second": second,
            "experimental_dG_B_minus_A_kcal_mol": truth[first] - truth[second],
            "solution_estimate_kcal_mol": solution["estimate_kcal_mol"],
            "solution_repeat_std_kcal_mol": solution["repeat_std_kcal_mol"],
            "solution_status_sha256": digest(solution_path), "bound_snapshots": snapshots}
    write_json(args.output, report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
