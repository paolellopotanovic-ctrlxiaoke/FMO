from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / "tools"))
from report_endpoint_energy_comparison import FIXED_ITERATIONS, HARTREE_TO_KCAL_MOL, grouped_summary


def test_matched_difference_cancels_arbitrary_molecule_offsets():
    frames = []
    for leg in ("solution", "bound"):
        for endpoint in (0, 10):
            change = (1 if endpoint == 0 else 4) if leg == "bound" else 0
            for iteration in FIXED_ITERATIONS:
                frames.append({"leg": leg, "endpoint": endpoint, "trajectory": leg,
                    "iteration": iteration, "MM_total_kcal_mol": 1000 * endpoint + change,
                    "DFTB_hartree": -2000 * endpoint + 2 * change / HARTREE_TO_KCAL_MOL})
    result = grouped_summary(frames)
    assert result["B_minus_A_internal_energy_double_difference_kcal_mol"]["MM"] == pytest.approx(3)
    assert result["B_minus_A_internal_energy_double_difference_kcal_mol"]["DFTB"] == pytest.approx(6)
    assert result["DFTB_minus_MM_double_difference_kcal_mol"] == pytest.approx(3)
    with pytest.raises(ValueError, match="incomplete"):
        grouped_summary(frames[:-1])
