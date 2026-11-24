from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from report_endpoint_complex_tpie import summarize


def fixture():
    rows = []
    for repeat in range(3):
        for name in ("A", "B"):
            for iteration in (400, 1200, 2000):
                value = 1 + repeat if name == "B" else 0
                pair = {"fragment_i": 2, "fragment_j": 1, "value_kcal_mol": value, "electrostatic_kcal_mol": value,
                        "exchange_repulsion_kcal_mol": 0, "charge_transfer_kcal_mol": 0,
                        "dispersion_kcal_mol": 0, "solvation_kcal_mol": 0}
                rows.append({"directory": f"{repeat}/{name}/{iteration}", "target": "T", "name": name,
                    "trajectory": str(repeat), "iteration": iteration, "ligand_fragment": 2,
                    "result": {"status": "SUCCEEDED", "tpie_kcal_mol": value, "pairs": [pair]}})
    return {"tasks": rows, "scope": "test"}, rows


def test_repeat_differences_not_independent_frames():
    plan, rows = fixture()
    report = summarize(plan, rows, {"T": ["A", "B"]}, {"T": {"A": 0, "B": 4}})
    target = report["targets"]["T"]
    assert target["n_independent_repeats"] == 3
    assert target["mean_B_minus_A_kcal_mol"]["tPIE"] == pytest.approx(2)
    assert target["repeat_std_B_minus_A_kcal_mol"]["tPIE"] == pytest.approx(1)
    assert target["receptor_fragment_B_minus_A_kcal_mol"][0]["mean"] == pytest.approx(2)
    assert not report["binding_free_energy_claim_allowed"]


def test_incomplete_never_publishes_partial_average():
    plan, rows = fixture()
    report = summarize(plan, rows[:-1], {"T": ["A", "B"]}, {})
    assert report["status"] == "INCOMPLETE"
    assert report["targets"] == {}


def test_duplicate_results_rejected():
    plan, rows = fixture()
    with pytest.raises(ValueError, match="duplicate"):
        summarize(plan, rows + rows[:1], {}, {})


def test_dense_schedule_supported_without_assuming_three_frames():
    plan, rows = fixture()
    plan["fixed_iterations"] = [400, 1200]
    rows = [row for row in rows if row["iteration"] != 2000]
    plan["tasks"] = rows
    report = summarize(plan, rows, {"T": ["A", "B"]}, {"T": {"A": 0, "B": 4}})
    assert report["status"] == "COMPLETE"
    assert len(report["targets"]["T"]["fixed_time_slice_sensitivity"]) == 2


def test_result_identity_must_match_frozen_task():
    plan, rows = fixture()
    changed = [{**row} for row in rows]
    changed[0]["iteration"] = 800
    with pytest.raises(ValueError, match="identity"):
        summarize(plan, changed, {}, {})
