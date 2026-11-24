from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / "tools"))
from report_endpoint_repeat_stability import HARTREE_TO_KCAL_MOL, paired_repeats, stability_summary


def frames_for_repeats(bound=(2., 4.), solution=(0., 1., -1.)):
    frames = []
    for leg, values in (("bound", bound), ("solution", solution)):
        for repeat, value in enumerate(values):
            for endpoint in (0, 10):
                for iteration in (400, 800, 1600, 2000):
                    energy = 1000 * endpoint + value * (endpoint == 10)
                    frames.append({"leg": leg, "trajectory": f"{leg}_{repeat}", "endpoint": endpoint,
                        "iteration": iteration, "MM_total_kcal_mol": energy,
                        "DFTB_hartree": 2 * energy / HARTREE_TO_KCAL_MOL})
    return frames


def test_repeat_statistics_do_not_count_cross_products():
    report = stability_summary(frames_for_repeats())
    assert report["repeat_counts"] == {"bound": 2, "solution": 3}
    assert not report["all_planned_repeats_present"]
    assert report["n_independent_bound_solution_combinations"] is None
    mm = report["models"]["MM"]
    assert mm["internal_energy_double_difference_kcal_mol"] == pytest.approx(3)
    assert mm["leg_repeat_std_kcal_mol"]["bound"] == pytest.approx(2 ** .5)
    assert mm["repeat_standard_error_diagnostic_kcal_mol"] == pytest.approx((1 + 1 / 3) ** .5)
    assert len(mm["leave_one_repeat_out"]) == 5
    assert report["models"]["model_difference"]["internal_energy_double_difference_kcal_mol"] == pytest.approx(3)


def test_single_bound_repeat_has_unknown_not_zero_variability():
    result = stability_summary(frames_for_repeats(bound=(2.,)))
    assert result["models"]["MM"]["leg_repeat_std_kcal_mol"]["bound"] is None
    assert result["models"]["MM"]["repeat_standard_error_diagnostic_kcal_mol"] is None


def test_missing_endpoint_fails_closed():
    with pytest.raises(ValueError, match="matching times"):
        paired_repeats(frames_for_repeats()[:-1])
