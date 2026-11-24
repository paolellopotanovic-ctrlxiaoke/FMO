import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from analyze_cycle_solvation_replacement import corrected_scores


def test_pair_solvent_is_removed_before_cycle_solvent_added_once():
    scores = corrected_scores(-100., 10., 40., 35., 5.)
    assert scores["electronic_pair_plus_cycle_electrostatic"] == -70.
    assert scores["electronic_pair_plus_cycle_total_solvation"] == -75.
    assert scores["electronic_pair_plus_cycle_total_solvation_and_strain"] == -70.


def test_incomplete_free_minimum_cannot_supply_strain_metric():
    scores = corrected_scores(-100., 10., 40., 35.)
    assert "electronic_pair_plus_cycle_total_solvation_and_strain" not in scores
