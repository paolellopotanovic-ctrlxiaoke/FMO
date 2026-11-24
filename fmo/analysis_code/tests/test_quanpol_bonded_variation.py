import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

from analyze_quanpol_bonded_variation import bonded_observables, describe_correlation, physical_constraint_pairs


def test_constraints_use_system_indices_not_saved_coordinate_order():
    import openmm

    system = openmm.System()
    for _ in range(4):
        system.addParticle(1)
    system.addConstraint(0, 3, .1)
    assert physical_constraint_pairs(system, [0, 3, 1, 2], [0, 1]) == {(0, 1)}
    with pytest.raises(ValueError, match="mapping"):
        physical_constraint_pairs(system, [0, 0], [0, 1])


def test_bond_and_angle_units_match_analytical_energy():
    import openmm
    from rdkit import Chem

    molecule = Chem.MolFromSmiles("CO")
    molecule = Chem.AddHs(molecule)
    system = openmm.System()
    bond = openmm.HarmonicBondForce()
    bond.addBond(0, 1, .1, 418.4)
    angle = openmm.HarmonicAngleForce()
    angle.addAngle(0, 1, 2, np.pi / 2, 8.368)
    system.addForce(bond)
    system.addForce(angle)
    terms = bonded_observables(system, molecule, np.array([[0., 0., 0.], [2., 0., 0.], [3., 1., 0.]]), {(0, 1)})
    assert terms[0]["energy_kcal_mol"] == pytest.approx(.5)
    assert terms[0]["constrained_in_source"] is True
    assert terms[0]["has_hydrogen"] is False
    assert terms[1]["has_hydrogen"] is True
    assert terms[1]["observed"] == pytest.approx(135)
    assert terms[1]["energy_kcal_mol"] == pytest.approx((np.pi / 4)**2)


def test_constant_term_does_not_emit_nan_correlation():
    assert describe_correlation([0, 0, 0], [1, 2, 3])["pearson_with_internal_delta_descriptive"] is None
    assert describe_correlation([1, 2, 3], [3, 2, 1])["pearson_with_internal_delta_descriptive"] == pytest.approx(-1)
