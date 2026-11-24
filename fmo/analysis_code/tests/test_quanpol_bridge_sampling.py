import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from sample_quanpol_bridge import GAS_KJ, compose_gradient, log_acceptance
from sample_quanpol_bridge import read_endpoint_state
from validate_quanpol_embedded_forces import direction_vector
from diagnose_quanpol_periodic_force_gate import evaluate_cutoff_jumps
from diagnose_quanpol_periodic_force_gate import custom_pair_cutoff_energy


def test_metropolis_uses_mm_integration_error_and_exact_target_difference():
    value = log_acceptance(2., -1., 3., .1, 298.15)
    assert value == pytest.approx(-(1. + .1 * 3. * 4.184) / (GAS_KJ * 298.15))
    assert log_acceptance(-2., 1., -3., .1, 298.15) == pytest.approx(-value)
    assert log_acceptance(2., -2., 1e6, 0., 298.15) == 0


@pytest.mark.parametrize("args", [(0, 0, float("nan"), .1, 298), (0, 0, 1, 1.1, 298), (0, 0, 1, .1, 0)])
def test_metropolis_invalid_target_fails_closed(args):
    with pytest.raises(ValueError):
        log_acceptance(*args)


def test_local_full_mapping_scatter_not_saved_particle_order():
    full = np.zeros((5, 3))
    result = compose_gradient(full, np.array([[1, 2, 3], [4, 5, 6]]), [4, 1])
    np.testing.assert_allclose(result[4], [1, 2, 3])
    np.testing.assert_allclose(result[1], [4, 5, 6])
    np.testing.assert_allclose(full, 0)
    with pytest.raises(ValueError):
        compose_gradient(full, np.ones((2, 3)), [1, 1])


def test_collective_checks_include_environment_and_are_reproducible():
    direction = {"kind": "collective_environment", "seed": 47}
    first = direction_vector(direction, (9, 3), 3)
    np.testing.assert_allclose(first[:3], 0)
    assert np.linalg.norm(first) == pytest.approx(1)
    np.testing.assert_array_equal(first, direction_vector(direction, (9, 3), 3))


def test_cutoff_crossing_is_explicit_jump_not_a_force_change():
    prepared = {"first": np.array([0]), "second": np.array([1]),
                "separation": np.array([[.89999, 0, 0]]), "jump": np.array([.002]), "cutoff": .9}
    positions = np.zeros((2, 3))
    vector = np.array([[1., 0, 0], [0, 0, 0]])
    rows = evaluate_cutoff_jumps(prepared, positions, [0, 1], vector, [.001, .00001])
    assert rows[0]["crossing_pairs"] == 1
    assert rows[0]["predicted_gradient_discontinuity"] == pytest.approx(-1)
    assert rows[1]["crossing_pairs"] == 0
    reverse = evaluate_cutoff_jumps(prepared, positions, [0, 1], -vector, [.001])
    assert reverse[0]["predicted_gradient_discontinuity"] == pytest.approx(1)


def test_endpoint_reference_loader_normalizes_numeric_reporter_reference(monkeypatch):
    import openmmtools.utils

    data = {"thermodynamic_states/state10": {"thermodynamic_state": {"_Reporter__compatible_state": "thermodynamic_states/0", "temperature": 298.15}},
            "thermodynamic_states/state0": {"thermodynamic_state": {"standard_system": "validated XML", "temperature": 298.15}}}

    class Reporter:
        def read_dict(self, key):
            return data[key]

    monkeypatch.setattr(openmmtools.utils, "deserialize", lambda value: value)
    state = read_endpoint_state(Reporter(), 10)
    assert state["thermodynamic_state"]["standard_system"] == "validated XML"
    assert "_Reporter__compatible_state" in data["thermodynamic_states/state10"]["thermodynamic_state"]


def test_custom_steric_pair_cutoff_energy_is_not_omitted():
    import openmm

    system = openmm.System()
    system.addParticle(1)
    system.addParticle(1)
    force = openmm.CustomNonbondedForce("scale*4*sqrt(epsilon1*epsilon2)*((0.5*(sigma1+sigma2)/r)^12-(0.5*(sigma1+sigma2)/r)^6)")
    force.addGlobalParameter("scale", 1)
    force.addPerParticleParameter("sigma")
    force.addPerParticleParameter("epsilon")
    force.addParticle([.3, .5])
    force.addParticle([.4, .8])
    force.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffNonPeriodic)
    force.setCutoffDistance(.9)
    system.addForce(force)
    integrator = openmm.VerletIntegrator(.001)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName("Reference"))
    context.setParameter("scale", .7)
    evaluate = custom_pair_cutoff_energy(system, context, .9)
    expected = .7 * 4 * (.5 * .8)**.5 * ((.35 / .9)**12 - (.35 / .9)**6) / 4.184
    assert evaluate(0, 1) == pytest.approx(expected)
    assert evaluate(1, 0) == pytest.approx(expected)
