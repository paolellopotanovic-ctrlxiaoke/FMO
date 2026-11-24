import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from report_bace_extensions import paired_rank_summary


def test_paired_comparison_counts_all_pairs_and_exposes_uniform_improvement():
    result = paired_rank_summary([1, 2, 3, 4], [4, 3, 2, 1], [1, 2, 3, 4], n_bootstrap=40)
    assert result["pair_changes"] == {"corrected": 6, "regressed": 0, "unchanged": 0, "tied": 0}
    assert result["leave_one_out_spearman_delta_range"] == [2., 2.]
    assert result["paired_bootstrap_95_percent_delta"]["spearman"] == [2., 2.]


def test_no_change_has_zero_bootstrap_difference():
    result = paired_rank_summary([1, 1, 3, 4], [2, 1, 4, 3], [2, 1, 4, 3], n_bootstrap=40)
    assert result["pair_changes"]["tied"] == 1
    assert result["paired_bootstrap_95_percent_delta"]["spearman"] == [0., 0.]
