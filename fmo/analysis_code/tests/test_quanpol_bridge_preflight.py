import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

from preflight_quanpol_bridge import first_step, weights


def test_weights_are_shift_invariant_and_finite():
    assert weights([-10000, 0, 10000], 298.15) == [1, 0, 0]
    assert weights([1, 2, 3], 298.15) == pytest.approx(weights([1001, 1002, 1003], 298.15))
    assert weights([1, 2, 3], 298.15, 0) == pytest.approx([1/3]*3)


def test_first_step_limits_first_reference_ess_only():
    values = [-10, -5, 0, 5, 10]
    coupling = first_step(values, 298.15, .8)
    assert 0 < coupling < .1
    normalized = weights(values, 298.15, coupling)
    assert (1 / sum(value**2 for value in normalized)) / 5 == pytest.approx(.8)
    assert first_step([0, 0, 0], 298.15) == 1


def test_invalid_inputs_are_rejected():
    with pytest.raises(ValueError):
        weights([], 298.15)
    with pytest.raises(ValueError):
        weights([1], 0)
    with pytest.raises(ValueError):
        first_step([1, 2], 298.15, 1.1)
