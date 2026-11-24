"""Saved model comparisons preserve complete cohorts and honest missingness."""
import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("build_model_comparison", ROOT / "tools/build_model_comparison.py")
COMPARISON = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(COMPARISON)


def test_generated_model_report_matches_metrics():
    result = json.loads((ROOT / "docs/MODEL_COMPARISON.json").read_text())
    assert (ROOT / "docs/MODEL_COMPARISON.zh-CN.md").read_text() == COMPARISON.render_report(result)
    assert result["new_QM_or_Boltz_or_FEP_inference"] is False
    assert len(result["groups"]) == 14
    assert sum(group["n"] for group in result["groups"].values()) == 760
    assert sum(group["n"] for group in result["groups"].values() if group["FEPplus_21_4"] is not None) == 270


def test_missing_boltz_is_not_filled_with_zero_or_fep():
    result = json.loads((ROOT / "docs/MODEL_COMPARISON.json").read_text())
    missing = {"CHK1", "CDK8", "BACE", "HSP90", "PTP1B", "FAAH"}
    assert set(result["missing_complete_Boltz2"]) in (missing, set())
    for name in missing:
        group = result["groups"][name]
        if result["missing_complete_Boltz2"]:
            assert group["Boltz2"] is None
            assert group["Boltz2_status"] == "NO_COMPLETE_RESULT_LOCATED"
        else:
            assert group["Boltz2"] is not None
            assert group["Boltz2_status"] == "COMPLETE_SAVED_RESULTS"
        assert group["FEPplus_21_4"] is not None
    assert sum(result["groups"][name]["n"] for name in missing) == 117
    openbind = result["groups"]["OpenBind"]
    assert openbind["FEPplus_21_4"] is None
    assert openbind["Boltz2_status"] == "PUBLISHED_BENCHMARK_SCORES_NOT_LOCAL_RERUN"


def test_comparison_uses_current_fmo_and_weighted_metrics():
    result = json.loads((ROOT / "docs/MODEL_COMPARISON.json").read_text())
    history = json.loads((ROOT / "docs/HISTORY_METRICS.json").read_text())
    for name, group in result["groups"].items():
        assert group["FMO_raw"] == history["latest_2004_groups"][name]["raw"]
        assert group["FMO_softcap"] == history["latest_2004_groups"][name]["softcap"]
    for model in ("FMO_raw", "FMO_softcap", "Boltz2", "FEPplus_21_4"):
        for metric in ("pearson", "spearman", "kendall"):
            expected = sum(result["groups"][name]["n"] * result["groups"][name][model][metric]
                           for name in ("CDK2", "JNK1", "P38", "TYK2")) / 87
            assert result["original87_weighted"][model][metric] == expected


def test_duplicate_identities_are_rejected():
    with pytest.raises(ValueError, match="duplicate identities"):
        COMPARISON.unique_rows([{"name": "a"}, {"name": "a"}], "name")
    assert set(COMPARISON.unique_rows([{"name": 1}, {"name": 2}], "name")) == {"1", "2"}


def test_metrics_reject_partial_nonfinite_and_constant_arrays():
    for observed, predicted in (([1, 2, 3], [1, 2]), ([1, 2, 3], [1, 2, float("nan")])):
        with pytest.raises(ValueError):
            COMPARISON.metrics(observed, predicted)
    assert COMPARISON.metrics([1, 2, 3], [3, 2, 1])["pearson"] == pytest.approx(-1)
    with pytest.raises(ValueError, match="saved metric differs"):
        COMPARISON.check_metrics({"pearson": 0.5}, {"pearson": 0.6})
