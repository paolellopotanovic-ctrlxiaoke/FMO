import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

from report_sidechain_adaptive_controls import pair_errors


def test_pair_errors_distinguish_correct_inverted_and_tied_predictions():
    assert pair_errors([1., 2., 3.], [1., 2., 3.]) == ([0, 0, 0], [2, 2, 2])
    assert pair_errors([1., 2., 3.], [3., 2., 1.]) == ([2, 2, 2], [2, 2, 2])
    assert pair_errors([1., 2.], [1., 1.]) == ([1, 1], [1, 1])


def test_experimental_ties_are_not_comparable():
    assert pair_errors([1., 1., 2.], [3., 2., 1.]) == ([1, 1, 2], [1, 1, 2])
