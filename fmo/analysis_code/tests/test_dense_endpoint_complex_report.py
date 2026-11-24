from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from report_dense_endpoint_complex_tpie import is_water_fragment, repeat_diagnostics


def test_uses_repeat_variation_instead_of_frame_pseudoreplicates():
    values = np.repeat(np.array([[-1.], [2.], [5.]]), 41, axis=1)
    result = repeat_diagnostics(values)
    assert result["mean"] == pytest.approx(2)
    assert result["repeat_std"] == pytest.approx(3)
    assert result["leave_one_repeat_out_means"] == pytest.approx([3.5, 2, .5])
    assert not result["all_three_repeats_same_nonzero_sign"]
    assert not result["convergence_proven"]


def test_fixed_windows_do_not_select_best_time():
    values = np.tile(np.arange(41), (3, 1))
    result = repeat_diagnostics(values)
    assert result["fixed_windows"]["early_1_to_2p9_ns"]["mean"] == pytest.approx(9.5)
    assert result["fixed_windows"]["late_3_to_5_ns"]["mean"] == pytest.approx(30)
    assert result["prefix_sensitivity"]["41"]["mean"] == result["mean"]


def test_incomplete_dense_schedule_rejected():
    with pytest.raises(ValueError, match="complete"):
        repeat_diagnostics(np.ones((3, 40)))


def test_water_fragment_names_include_chain_prefix():
    assert is_water_fragment("AHOH604")
    assert is_water_fragment("HOH1028")
    assert not is_water_fragment("AASP93")
