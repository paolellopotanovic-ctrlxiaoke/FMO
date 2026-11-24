import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from report_intramolecular_qm_mm_correction import component_perturbation_summary, relative_cycle


def test_constant_component_shift_is_recovered():
    result = component_perturbation_summary([4.0] * 20, temperature_kelvin=298.15)
    assert result["delta_f_kcal_mol"] == pytest.approx(4.0)
    assert result["weight_ess"] == pytest.approx(20.0)
    assert result["convergence_proven"] is False


def test_large_component_offset_is_numerically_stable():
    first = component_perturbation_summary([0.0, 1.0, 2.0], temperature_kelvin=298.15)
    second = component_perturbation_summary([1e6, 1e6 + 1.0, 1e6 + 2.0], temperature_kelvin=298.15)
    assert second["delta_f_kcal_mol"] - first["delta_f_kcal_mol"] == pytest.approx(1e6)


def test_relative_cycle_cancels_ligand_specific_reference_shifts():
    base = relative_cycle(-2.0, 1.0, -1.0, 3.0)
    shifted = relative_cycle(-2.0 + 17.0, 1.0 - 23.0, -1.0 + 17.0, 3.0 - 23.0)
    assert base == pytest.approx(shifted)


@pytest.mark.parametrize("values", [[], [1.0], [0.0, math.nan], [[0.0, 1.0]]])
def test_invalid_component_samples_rejected(values):
    with pytest.raises(ValueError):
        component_perturbation_summary(values, temperature_kelvin=298.15)


def test_invalid_temperature_rejected():
    with pytest.raises(ValueError):
        component_perturbation_summary([0.0, 1.0], temperature_kelvin=0.0)
