from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from run_endpoint_local_water_control import select_water_geometry, water_groups


def test_nearest_water_selection_uses_periodic_whole_molecules():
    coordinates = np.array([[9.8, 0, 0], [.76, 0, 0], [9.5, .9, 0], [4., 0, 0], [4.96, 0, 0], [3.7, .9, 0]])
    result = select_water_geometry(coordinates, {"near": [0, 1, 2], "far": [3, 4, 5]}, np.array([[0., 0, 0]]), np.eye(3) * 10, 1)
    assert result[0]["residue"] == "near"
    assert result[0]["minimum_polar_distance_angstrom"] == pytest.approx(.2)
    positions = np.asarray(result[0]["xyz_angstrom"])
    assert np.linalg.norm(positions[1] - positions[0]) == pytest.approx(.96)


def test_invalid_water_bonds_do_not_get_repaired_silently():
    with pytest.raises(ValueError, match="geometry"):
        select_water_geometry(np.array([[0., 0, 0], [2., 0, 0], [0, .9, 0]]), {"bad": [0, 1, 2]},
                              np.array([[3., 0, 0]]), np.eye(3) * 10, 1)


def test_repeated_pdb_water_residue_numbers_do_not_merge_molecules():
    topology = [{"class": 0., "resname": "HOH", "residue": "A 1", "element": element}
                for element in ("O", "H", "H", "O", "H", "H")]
    assert list(water_groups(topology).values()) == [[0, 1, 2], [3, 4, 5]]


def test_partial_water_does_not_pass_inventory():
    topology = [{"class": 0., "resname": "HOH", "residue": "A 1", "element": "H"}]
    with pytest.raises(ValueError, match="incomplete"):
        water_groups(topology)
