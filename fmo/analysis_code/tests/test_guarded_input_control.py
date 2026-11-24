import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from report_guarded_input_control import select_candidate


def test_guard_depends_on_validity_not_score():
    assert select_candidate({"status": "SUCCEEDED", "tpie": 50.}, {"status": "SUCCEEDED", "tpie": -100.}) == (50., "new_input_heavy_fixed")
    assert select_candidate({"status": "FAILED", "error": "ValueError: QM status: FAILED_UNPHYSICAL_PIEDA"},
                            {"status": "SUCCEEDED", "tpie": -100.}) == (-100., "baseline_retained_after_unphysical_input_pose")


def test_unknown_failure_cannot_be_hidden_by_fallback():
    with pytest.raises(ValueError, match="unregistered failure"):
        select_candidate({"status": "FAILED", "error": "atom mapping mismatch"}, {"status": "SUCCEEDED", "tpie": -100.})
