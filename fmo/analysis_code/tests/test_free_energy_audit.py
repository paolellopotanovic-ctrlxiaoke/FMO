import numpy as np
import pytest

from fmotpiefinal.core.free_energy_audit import (
    perturbation_summary,
    relative_binding_difference,
    relative_qmmm_correction,
)


def test_constant_shift_and_ess():
    result = perturbation_summary([4.0] * 20, temperature_kelvin=298.15)
    assert result["delta_free_energy_kcal_mol"] == pytest.approx(4.0)
    assert result["weight_ess"] == pytest.approx(20.0)
    assert result["convergence_proven"] is False


def test_large_offsets_are_stable():
    values = np.array([0.0, 1.0, 2.0, 3.0])
    first = perturbation_summary(values, temperature_kelvin=300)
    second = perturbation_summary(values + 1e6, temperature_kelvin=300)
    assert second["delta_free_energy_kcal_mol"] - first["delta_free_energy_kcal_mol"] == pytest.approx(1e6)


def test_weight_collapse_is_visible():
    result = perturbation_summary([-100.0] + [0.0] * 99, temperature_kelvin=298.15)
    assert result["weight_ess"] == pytest.approx(1)
    assert result["largest_normalized_weight"] == pytest.approx(1)


@pytest.mark.parametrize("values", [[], [1], [0, float("nan")], [[0, 1]]])
def test_bad_arrays_rejected(values):
    with pytest.raises(ValueError):
        perturbation_summary(values, temperature_kelvin=298.15)


def test_interaction_only_is_not_total_hamiltonian():
    with pytest.raises(ValueError, match="total_potential"):
        perturbation_summary([0, 1], temperature_kelvin=298.15, energy_kind="tpie")


def test_leg_direction():
    assert relative_binding_difference(2, 5) == -3
    with pytest.raises(ValueError):
        relative_binding_difference(2, float("nan"))


def test_relative_qmmm_correction_uses_two_thermodynamic_legs():
    result = relative_qmmm_correction(
        [1.0] * 10, [2.0] * 10, [0.5] * 10, [0.5] * 10,
        temperature_kelvin=298.15, energy_kind="total_potential")
    assert result["delta_binding_correction_kcal_mol"] == pytest.approx(1.0)
    assert result["convergence_proven"] is False
    assert result["ranking_eligible"] is False


def test_relative_qmmm_correction_rejects_interaction_energy():
    with pytest.raises(ValueError, match="total_potential"):
        relative_qmmm_correction(
            [1.0, 1.0], [2.0, 2.0], [0.5, 0.5], [0.5, 0.5],
            temperature_kelvin=298.15, energy_kind="tpie")
