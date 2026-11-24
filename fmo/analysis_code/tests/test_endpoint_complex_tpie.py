from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from run_endpoint_complex_tpie import map_receptor, whole_components
from fmotpiefinal.core.structure import Atom


def test_receptor_mapping_tracks_renumbered_water_and_preserves_identity():
    atoms = [Atom(0, "C", (1, 0, 0), "A", "ALA", "1", atom_name="CA"),
             Atom(1, "O", (3, 0, 0), "A", "HOH", "100", "HETATM", "O")]
    topology = [{"name": "CA", "resname": "ALA", "residue": "A   1", "element": "C", "class": 0., "xyz": [1, 0, 0]},
                {"name": "O", "resname": "HOH", "residue": "B 200", "element": "O", "class": 0., "xyz": [3, 0, 0]}]
    indices, modes = map_receptor(atoms, topology)
    assert indices.tolist() == [0, 1]
    assert modes == ["residue_atom_identity", "verified_initial_geometry"]
    with pytest.raises(ValueError, match="unique receptor mapping"):
        map_receptor(atoms, topology[:1])


def test_mapping_rejects_duplicate_atoms_instead_of_guessing():
    atom = Atom(0, "C", (1, 0, 0), "A", "ALA", "1", atom_name="CA")
    entry = {"name": "CA", "resname": "ALA", "residue": "A 1", "element": "C", "class": 0., "xyz": [1, 0, 0]}
    with pytest.raises(ValueError, match="ambiguous"):
        map_receptor([atom], [entry, entry])


def test_whole_components_keeps_bonded_protein_and_water_in_ligand_image():
    xyz = np.array([[9.7, 0, 0], [.7, 0, 0], [4., 0, 0], [4.8, 0, 0]])
    whole = whole_components(xyz, [(0, 1), (2, 3)], np.eye(3) * 10, np.array([0., 0, 0]), [0])
    assert np.linalg.norm(whole[1] - whole[0]) == pytest.approx(1.)
    assert np.linalg.norm(whole[3] - whole[2]) == pytest.approx(.8)
    assert abs(whole[0, 0]) < 1
    assert abs(whole[2, 0]) <= 5


def test_inconsistent_periodic_bond_cycle_fails_closed():
    xyz = np.array([[1., 0, 0], [4., 0, 0], [7., 0, 0]])
    with pytest.raises(ValueError, match="bond cycle"):
        whole_components(xyz, [(0, 1), (1, 2), (2, 0)], np.eye(3) * 10, np.zeros(3), [0])
