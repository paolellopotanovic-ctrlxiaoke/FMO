from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from report_endpoint_conformer_tpie import summarize_control


def test_missing_or_failed_results_never_publish_scores():
    plan = {"tasks": [{"directory": "first"}, {"directory": "second"}], "scope": "control"}
    row = {"directory": "first", "result": {"status": "FAILED"}}
    result = summarize_control(plan, [row])
    assert result["status"] == "INCOMPLETE"
    assert result["targets"] == {}
    assert result["missing"] == ["second"]


def test_duplicate_results_rejected():
    plan = {"tasks": [{"directory": "first"}], "scope": "control"}
    row = {"directory": "first", "result": {"status": "SUCCEEDED"}}
    with pytest.raises(ValueError, match="duplicate"):
        summarize_control(plan, [row, row])
