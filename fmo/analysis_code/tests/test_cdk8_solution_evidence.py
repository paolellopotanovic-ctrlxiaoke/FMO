"""The published development result must retain failures and unmet goals."""
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_cdk8_report_keeps_complete_cohort_and_negative_result():
    report = json.loads((ROOT / "docs/CDK8_SOLUTION_CYCLE_RESULTS.json").read_text())
    assert report["n"] == len(report["rows"]) == 10
    assert {row["ligand"] for row in report["rows"]} == {"12", "13", "14", "15", "16", "17", "19", "20", "21", "22"}
    assert not report["goal_achieved"]
    assert not report["production_changed"]
    assert len(report["failed_gas_jobs"]) == 10
    assert report["native_max_abs_reconciliation_kcal_mol"] < 0.001
    assert report["metrics"]["tpie_pcm"]["spearman"] < 0.24


def test_multistart_evidence_is_full_not_success_subset():
    report = json.loads((ROOT / "docs/CDK8_SOLUTION_CYCLE_RESULTS.json").read_text())
    sampling = report["sampling"]
    assert sampling["n_jobs"] == 40
    assert not sampling["failed_starts"]
    assert all(row["n_successful_minima_including_bound_start"] == 5 for row in sampling["rows"])
    assert sampling["metrics"]["tpie_plus_multistart_strain"]["spearman"] < 0.6
