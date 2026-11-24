"""Audit an isolated-ligand QM/MM internal-potential replacement cycle."""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp


HARTREE_TO_KCAL_MOL = 627.5094740631
GAS_CONSTANT_KCAL_MOL_K = 0.00198720425864083
FIXED_ITERATIONS = tuple(range(400, 2001, 40))


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def component_perturbation_summary(delta_u_kcal_mol, *, temperature_kelvin):
    """One-sided FEP for replacing one explicitly named Hamiltonian component."""
    values = np.asarray(delta_u_kcal_mol, dtype=float)
    if values.ndim != 1 or len(values) < 2 or not np.isfinite(values).all():
        raise ValueError("at least two finite component-energy differences required")
    if not math.isfinite(temperature_kelvin) or temperature_kelvin <= 0:
        raise ValueError("positive finite temperature required")
    thermal_energy = GAS_CONSTANT_KCAL_MOL_K * temperature_kelvin
    log_weights = -values / thermal_energy
    log_normalizer = logsumexp(log_weights)
    normalized_weights = np.exp(log_weights - log_normalizer)
    effective_samples = float(1.0 / np.sum(normalized_weights ** 2))
    return {
        "n_frames": len(values),
        "delta_f_kcal_mol": float(-thermal_energy * (log_normalizer - math.log(len(values)))),
        "mean_delta_u_kcal_mol": float(values.mean()),
        "sd_delta_u_kcal_mol": float(values.std(ddof=1)),
        "weight_ess": effective_samples,
        "weight_ess_fraction": effective_samples / len(values),
        "largest_normalized_weight": float(normalized_weights.max()),
        "convergence_proven": False,
    }


def relative_cycle(bound_a, bound_b, solution_a, solution_b):
    return (bound_b - bound_a) - (solution_b - solution_a)


def saved_temperature(trajectory):
    from netCDF4 import Dataset, chartostring
    import yaml

    temperatures = []
    with Dataset(trajectory) as data:
        for state_index in (1, 10):
            serialized = chartostring(data[f"thermodynamic_states/state{state_index}"][:]).item()
            state = yaml.load(serialized, Loader=yaml.BaseLoader)
            temperature = state["thermodynamic_state"]["temperature"]
            if temperature["unit"] != "kelvin":
                raise ValueError("saved thermodynamic-state temperature must use kelvin")
            temperatures.append(float(temperature["value"]))
    if not math.isclose(temperatures[0], temperatures[1], rel_tol=0, abs_tol=1e-10):
        raise ValueError("endpoint thermodynamic-state temperatures differ")
    return temperatures[0]


def _leg_summary(rows, *, temperature_kelvin):
    rows = sorted(rows, key=lambda row: (row["trajectory"], row["iteration"]))
    values = np.asarray([row["delta_u_centered_kcal_mol"] for row in rows], dtype=float)
    summary = component_perturbation_summary(values, temperature_kelvin=temperature_kelvin)
    return {
        **summary,
        "n_trajectories": len({row["trajectory"] for row in rows}),
        "replicates": [],
        "early_1_to_3ns": None,
        "late_3p1_to_5ns": None,
    }


def build_report(mm_path, qm_path, protocol_path, openfe_run, trajectory_base):
    mm = json.loads(Path(mm_path).read_text())
    qm = json.loads(Path(qm_path).read_text())
    plan = json.loads(Path(protocol_path).read_text())
    if mm.get("status") != "COMPLETE":
        raise ValueError("complete matched MM internal-energy report required")
    if (qm.get("status") != "COMPLETE" or qm.get("n_requested") != qm.get("n_succeeded")
            or qm.get("n_requested") != len(qm.get("results", []))):
        raise ValueError("complete DFTB endpoint set required; failed frames cannot be dropped")
    if plan.get("diagnostics_sha256") != mm.get("diagnostics_sha256"):
        raise ValueError("QM and MM reports do not share frozen trajectory diagnostics")
    if plan.get("method") != "isolated complete ligand FMO1-DFTB3 gas phase, DC=.TRUE. IDCVER=4":
        raise ValueError("unexpected high-level component definition")
    if qm["n_requested"] != len(plan.get("tasks", [])):
        raise ValueError("QM task count differs from frozen schedule")

    task_by_key = {(task["target"], task["directory"]): task for task in plan["tasks"]}
    if len(task_by_key) != len(plan["tasks"]):
        raise ValueError("duplicate frozen QM task identity")
    result_by_key = {(row["target"], row["directory"]): row for row in qm["results"]}
    if (len(result_by_key) != len(qm["results"]) or set(result_by_key) != set(task_by_key)
            or any(row.get("status") != "SUCCEEDED" for row in qm["results"])):
        raise ValueError("QM results must match every frozen task exactly once")

    mm_by_frame = {}
    for target, target_data in mm["targets"].items():
        for leg in ("bound", "solution"):
            for repeat in target_data["legs"][leg]["repeats"]:
                for frame in repeat["frames"]:
                    key = (target, leg, repeat["endpoint"], repeat["trajectory"], frame["iteration"])
                    if key in mm_by_frame:
                        raise ValueError("duplicate MM frame identity")
                    mm_by_frame[key] = frame

    joined = defaultdict(list)
    trajectory_paths = set()
    for key, task in task_by_key.items():
        row = result_by_key[key]
        frame_key = (row["target"], row["leg"], row["endpoint"], task["trajectory"], row["iteration"])
        frame = mm_by_frame.get(frame_key)
        if frame is None or frame["replica"] != row["replica"]:
            raise ValueError(f"QM/MM physical-frame identity mismatch: {frame_key}")
        q_energy = float(row["energy_hartree"]) * HARTREE_TO_KCAL_MOL
        mm_internal = float(frame["energies_kcal_mol"]["total"])
        trajectory_path = (Path(trajectory_base) / task["trajectory"]).resolve(strict=True)
        trajectory_paths.add(trajectory_path)
        repeat_key = (row["target"], row["leg"], row["endpoint"], task["trajectory"])
        joined[repeat_key].append({
            "iteration": row["iteration"],
            "replica": row["replica"],
            "trajectory": task["trajectory"],
            "delta_u_kcal_mol": q_energy - mm_internal,
        })

    temperatures = {str(path): saved_temperature(path) for path in sorted(trajectory_paths)}
    unique_temperatures = {round(value, 8) for value in temperatures.values()}
    if len(unique_temperatures) != 1:
        raise ValueError(f"trajectory temperatures differ: {sorted(unique_temperatures)}")
    temperature = unique_temperatures.pop()

    target_reports = {}
    for target, mm_target in mm["targets"].items():
        endpoint_names = {
            0: mm_target["endpoints"]["0"]["name"],
            10: mm_target["endpoints"]["10"]["name"],
        }
        centers = {}
        for endpoint in (0, 10):
            endpoint_values = [item["delta_u_kcal_mol"]
                               for (name, _, ep, _), rows in joined.items()
                               if name == target and ep == endpoint for item in rows]
            centers[endpoint] = float(np.median(endpoint_values))

        endpoint_data = {}
        leg_dg = {}
        for leg in ("bound", "solution"):
            for endpoint in (0, 10):
                repeat_groups = [
                    (trajectory, rows)
                    for (name, current_leg, ep, trajectory), rows in joined.items()
                    if name == target and current_leg == leg and ep == endpoint
                ]
                repeat_reports = []
                leg_rows = []
                for trajectory, rows in sorted(repeat_groups):
                    ordered = sorted(rows, key=lambda item: item["iteration"])
                    if tuple(item["iteration"] for item in ordered) != FIXED_ITERATIONS:
                        raise ValueError(f"incomplete fixed QM/MM frame grid: {target}/{leg}/{endpoint}/{trajectory}")
                    centered = [{**item, "delta_u_centered_kcal_mol": item["delta_u_kcal_mol"] - centers[endpoint]}
                                for item in ordered]
                    leg_rows.extend(centered)
                    repeat_summary = component_perturbation_summary(
                        [item["delta_u_centered_kcal_mol"] for item in centered],
                        temperature_kelvin=temperature)
                    repeat_reports.append({
                        "trajectory": trajectory,
                        "n_frames": len(centered),
                        "delta_f_kcal_mol": repeat_summary["delta_f_kcal_mol"],
                        "weight_ess_fraction": repeat_summary["weight_ess_fraction"],
                        "largest_normalized_weight": repeat_summary["largest_normalized_weight"],
                    })
                if not repeat_reports:
                    raise ValueError(f"missing endpoint ensemble: {target}/{leg}/{endpoint}")
                pooled = _leg_summary(leg_rows, temperature_kelvin=temperature)
                pooled["replicates"] = repeat_reports
                for name, lower, upper in (("early_1_to_3ns", 400, 1200), ("late_3p1_to_5ns", 1240, 2000)):
                    selected = [item for item in leg_rows if lower <= item["iteration"] <= upper]
                    pooled[name] = component_perturbation_summary(
                        [item["delta_u_centered_kcal_mol"] for item in selected],
                        temperature_kelvin=temperature)
                endpoint_data[f"{leg}_{'A' if endpoint == 0 else 'B'}"] = {
                    "ligand": endpoint_names[endpoint],
                    "median_reference_shift_kcal_mol": centers[endpoint],
                    **pooled,
                }
                leg_dg[(leg, endpoint)] = pooled["delta_f_kcal_mol"]

        correction = relative_cycle(
            leg_dg[("bound", 0)], leg_dg[("bound", 10)],
            leg_dg[("solution", 0)], leg_dg[("solution", 10)],
        )
        phase_corrections = {}
        for phase_key, phase_name in (("early_1_to_3ns", "early_1_to_3ns"),
                                      ("late_3p1_to_5ns", "late_3p1_to_5ns")):
            phase_dg = {
                ("bound", 0): endpoint_data["bound_A"][phase_key]["delta_f_kcal_mol"],
                ("bound", 10): endpoint_data["bound_B"][phase_key]["delta_f_kcal_mol"],
                ("solution", 0): endpoint_data["solution_A"][phase_key]["delta_f_kcal_mol"],
                ("solution", 10): endpoint_data["solution_B"][phase_key]["delta_f_kcal_mol"],
            }
            phase_corrections[phase_name] = relative_cycle(
                phase_dg[("bound", 0)], phase_dg[("bound", 10)],
                phase_dg[("solution", 0)], phase_dg[("solution", 10)],
            )
        bound_status = json.loads((Path(openfe_run) / "production" / target / "edge_000" / "bound" / "status.json").read_text())
        solution_status = json.loads((Path(openfe_run) / "production" / target / "edge_000" / "solution" / "status.json").read_text())
        mm_relative = float(bound_status["estimate_kcal_mol"]) - float(solution_status["estimate_kcal_mol"])
        target_reports[target] = {
            "edge": [endpoint_names[0], endpoint_names[10]],
            "component_definition": "gas-phase isolated-ligand DFTB3 internal potential minus the matched isolated-ligand OpenFF MM internal potential; added as a replacement of only that intramolecular term while retaining the original MM environment Hamiltonian",
            "endpoint_component_corrections": endpoint_data,
            "relative_component_correction_kcal_mol": correction,
            "time_window_relative_component_correction_kcal_mol": phase_corrections,
            "time_window_drift_kcal_mol": phase_corrections["late_3p1_to_5ns"] - phase_corrections["early_1_to_3ns"],
            "openfe_mm_edge_estimate_bound_minus_solution_kcal_mol": mm_relative,
            "partial_component_adjusted_edge_estimate_kcal_mol": mm_relative + correction,
            "openfe_bound_status_claim_allowed": bound_status.get("equilibrium_affinity_claim_allowed"),
            "openfe_solution_status_claim_allowed": solution_status.get("equilibrium_affinity_claim_allowed"),
            "convergence_proven": False,
            "ranking_eligible": False,
        }

    return {
        "schema": "isolated-ligand-intramolecular-dftb-mm-four-leg-correction-v1",
        "temperature_kelvin": temperature,
        "n_qm_successful": len(qm["results"]),
        "n_mm_qm_frames_matched": sum(len(rows) for rows in joined.values()),
        "n_unique_simulation_trajectories": len({(key[0], key[1], key[3]) for key in joined}),
        "n_endpoint_trajectory_groups": len(joined),
        "trajectory_temperature_audit_kelvin": temperatures,
        "sources": {str(Path(path).resolve()): sha256(path) for path in (mm_path, qm_path, protocol_path)},
        "selection": "All existing completed OpenFE repeats and all 41 pre-registered frames from 1.0 to 5.0 ns at 0.1 ns spacing; no affinity-guided frame selection; no new MD.",
        "labels_used": False,
        "main_pipeline_changed": False,
        "interpretation": "This is a four-leg free-energy perturbation of the ligand intramolecular potential component only. It is not an FMO interaction score and not a full QM/MM total-potential correction.",
        "limitations": [
            "The high-level term is isolated gas-phase ligand FMO1-DFTB3; receptor/water electronic polarization, ligand-environment QM/MM coupling, and environment Hamiltonian corrections are absent.",
            "One-sided exponential averaging has poor weight ESS in several endpoint ensembles and does not establish overlap or convergence.",
            "Saved frames are time-correlated; weight ESS is not an independent-sample count.",
            "CDK8 has only one complete bound repeat in this saved production set; one additional repeat was incomplete and was not included.",
            "The pre-existing OpenFE edge estimates are marked as not eligible for equilibrium-affinity claims; adjusted edge values are diagnostic only.",
            "A two-ligand edge cannot yield a meaningful cohort ranking statistic such as Spearman correlation.",
        ],
        "targets": target_reports,
        "convergence_proven": False,
        "ranking_eligible": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mm", type=Path, required=True)
    parser.add_argument("--qm", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--openfe-run", type=Path, required=True)
    parser.add_argument("--trajectory-base", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = build_report(args.mm, args.qm, args.protocol, args.openfe_run, args.trajectory_base)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({target: data["relative_component_correction_kcal_mol"]
                      for target, data in report["targets"].items()}, sort_keys=True))


if __name__ == "__main__":
    main()
