import pytest

from fmotpiefinal.core.interaction_audit import audit_rows


def _row(**overrides):
    row = {"case": "x", "tpie": -10.0, "ES": -8.0, "EX": 2.0,
           "CT": -1.0, "DI": -3.0, "SOLV": 0.0,
           "three_state": -4.0, "pb_polar": 2.0}
    row.update(overrides)
    return row


def test_closure_and_overlap_are_reported_without_fitting():
    report = audit_rows([_row()])
    assert report["status"] == "PASS_BOOKKEEPING"
    assert report["component_closure"]["n_pass"] == 1
    assert report["ranking_eligible"] is False
    assert report["label_free_diagnostics"]["optional_observable_coverage"]["three_state"] == 1
    assert any("double-counts" in warning for warning in report["overlap_warnings"])


def test_component_closure_failure_is_fail_closed():
    report = audit_rows([_row(tpie=-9.8)], closure_tolerance=0.05)
    assert report["status"] == "REVIEW_COMPONENT_CLOSURE"
    assert report["component_closure"]["n_fail"] == 1


def test_missing_tpie_is_allowed_only_with_complete_components():
    report = audit_rows([_row(tpie=None)])
    assert report["component_closure"]["n_complete"] == 0
    assert report["status"] == "REVIEW_INCOMPLETE_COMPONENT_COVERAGE"
    with pytest.raises(ValueError, match="tPIE or all five"):
        audit_rows([{"case": "bad", "ES": -1.0}])
