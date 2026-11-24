import sys
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"tools"))
from run_direction_preserving_hydrogen_control import radial_coordinates


def test_direction_preserved_while_only_reference_bond_length_is_transplanted():
    molecule = Chem.AddHs(Chem.MolFromSmiles("[OH]"))
    before = ["O 8 0 0 0", "H 1 1 0 0"]
    reference = ["O 8 0 0 0", "H 1 0 .95 0"]
    result = radial_coordinates(molecule, before, reference)
    assert np.allclose(list(map(float, result[1].split()[2:])), [.95, 0., 0.])
    assert np.allclose(list(map(float, result[0].split()[2:])), [0., 0., 0.])


def test_bad_reference_bond_length_is_not_silently_clamped():
    molecule = Chem.AddHs(Chem.MolFromSmiles("[OH]"))
    with pytest.raises(ValueError, match="unphysical"):
        radial_coordinates(molecule, ["O 8 0 0 0", "H 1 1 0 0"], ["O 8 0 0 0", "H 1 0 10 0"])


def test_chemical_inventory_mismatch_is_rejected():
    with pytest.raises(ValueError, match="inventory mismatch"):
        radial_coordinates(Chem.MolFromSmiles("O"), ["O 8 0 0 0", "H 1 1 0 0"], ["O 8 0 0 0", "H 1 0 .95 0"])
