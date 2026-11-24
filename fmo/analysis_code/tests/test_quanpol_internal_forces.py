import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

from preflight_quanpol_internal_forces import BOHR_ANGSTROM, HARTREE_KCAL_MOL, parse_gradient, tangent_projection


def test_gradient_parser_units_and_atom_identity():
    text = "FINAL RDFTB ENERGY IS -1.0 AFTER 2 ITERATIONS\nGRADIENT OF THE ENERGY\nUNITS ARE HARTREE/BOHR\n 1 C 1.0 0.0 0.0\nEXECUTION OF GAMESS TERMINATED NORMALLY\n"
    assert parse_gradient(text, ["C"])[0, 0] == pytest.approx(HARTREE_KCAL_MOL / BOHR_ANGSTROM)
    with pytest.raises(ValueError, match="order"):
        parse_gradient(text, ["H"])


def test_constraint_projection_closes_and_removes_bond_stretch():
    coordinates = np.array([[0., 0., 0.], [1., 0., 0.]])
    gradient = np.array([[-1., 2., 0.], [1., 3., 0.]])
    tangent, normal = tangent_projection(gradient, coordinates, [(0, 1)])
    np.testing.assert_allclose(tangent[:, 0], 0, atol=1e-12)
    np.testing.assert_allclose(tangent + normal, gradient)
    assert float(np.sum(tangent * normal)) == pytest.approx(0, abs=1e-12)
