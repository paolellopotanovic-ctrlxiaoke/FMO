import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from report_quanpol_bridge_sampling import adjacent_diagnostic, normalized_log_weights, summarize


def test_log_weights_stable_and_constant_has_full_weight_ess():
    values, log_mean = normalized_log_weights([100000, 100000, 100000])
    np.testing.assert_allclose(values, 1 / 3)
    assert log_mean == 100000
    result = adjacent_diagnostic([-30000.] * 12, [-30000.] * 12, 0, .05, 298.15, -30000)
    assert result["forward_weight_ess_uncorrected"] == pytest.approx(12)
    assert result["difference_forward_reverse_kcal_mol"] == 0
    assert result["diagnostic_only_not_equilibrated_or_autocorrelation_corrected"]


def test_unstarted_and_failed_states_are_not_reported_complete(tmp_path):
    state = tmp_path / "b10"
    state.mkdir()
    (state / "protocol.json").write_text(json.dumps({"status": "ERROR", "error": "loader failed"}))
    result = summarize(tmp_path)
    assert result["status"] == "ERROR"
    assert result["states"]["b0"]["status"] == "PREPARING_OR_QUEUED"
    assert result["full_four_state_free_energy_correction_kcal_mol"] is None
    assert result["new_target_ranking_metrics"] is None
