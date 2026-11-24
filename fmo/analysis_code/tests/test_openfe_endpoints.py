import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / "tools"))
from analyze_openfe_endpoints import circular_summary, endpoint_indices, fitted_rmsd, minimum_image, torsion_degrees, unwrap_ligand


def test_periodic_minimum_image():
    assert minimum_image(np.array([9.0, 0.0, 0.0]), np.eye(3) * 10).tolist() == pytest.approx([-1, 0, 0])


def test_rmsd_removes_rigid_rotation_not_reflection():
    reference = np.array([[0., 0., 0.], [1., 0., 0.], [0., 2., 0.], [0., 0., 3.]])
    rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    assert fitted_rmsd(reference @ rotation + 4, reference) < 1e-12
    assert fitted_rmsd(reference * [-1, 1, 1], reference) > 0.1


def test_endpoint_excludes_dummy_atoms_and_hydrogens():
    atoms = [{"class": atom_class, "element": element}
             for atom_class, element in [(0, "C"), (.25, "C"), (.5, "N"), (.75, "O"), (.5, "H")]]
    assert endpoint_indices(atoms, 0).tolist() == [1, 2]
    assert endpoint_indices(atoms, 10).tolist() == [2, 3]


def test_torsion_periodic_and_circular_summary():
    coordinates = np.array([[0., 1., 0.], [0., 0., 0.], [1., 0., 0.], [1., 0., 1.]])
    shifted = coordinates.copy()
    shifted[-1] += [10, 0, 0]
    assert abs(torsion_degrees(coordinates, np.eye(3) * 10)) == pytest.approx(90)
    assert torsion_degrees(coordinates, np.eye(3) * 10) == pytest.approx(torsion_degrees(shifted, np.eye(3) * 10))
    result = circular_summary([179, -179])
    assert abs(result["circular_mean_degrees"]) == pytest.approx(180)
    assert result["circular_std_degrees"] < 2


def test_bond_unwrap_does_not_fold_long_ligand():
    coordinates = np.array([[1., 0., 0.], [4., 0., 0.], [7., 0., 0.], [0., 0., 0.]])
    whole = unwrap_ligand(coordinates, [(0, 1), (1, 2), (2, 3)], np.eye(3) * 10)
    assert whole[:, 0].tolist() == pytest.approx([1, 4, 7, 10])
    with pytest.raises(ValueError, match="disconnected"):
        unwrap_ligand(coordinates, [(0, 1)], np.eye(3) * 10)
