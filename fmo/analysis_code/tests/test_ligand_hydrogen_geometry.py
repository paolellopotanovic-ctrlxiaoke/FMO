import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from fmotpiefinal.core.ligand_hydrogen_geometry import audit_hch_angles, rebuild_flagged_carbon_hydrogens


def test_tetrahedral_angles_pass_and_collapsed_hch_is_flagged():
    molecule = Chem.AddHs(Chem.MolFromSmiles("C"))
    AllChem.EmbedMolecule(molecule, randomSeed=7)
    positions = np.asarray(molecule.GetConformer().GetPositions())
    assert audit_hch_angles(molecule, positions)["n_gross_distortions"] == 0
    positions[2] = positions[1]+np.array([.01, .01, .01])
    result = audit_hch_angles(molecule, positions)
    assert result["n_gross_distortions"] >= 1
    assert result["angles"][0]["atoms_zero_based"] == [1, 0, 2]


def test_geometry_screen_rejects_bad_shape_and_thresholds():
    molecule = Chem.AddHs(Chem.MolFromSmiles("C"))
    with pytest.raises(ValueError):
        audit_hch_angles(molecule, np.zeros((1, 3)))
    with pytest.raises(ValueError):
        audit_hch_angles(molecule, np.zeros((5, 3)), minimum_degrees=float("nan"))


def test_non_sp3_angles_are_not_subject_to_tetrahedral_screen():
    molecule = Chem.AddHs(Chem.MolFromSmiles("C=C"))
    AllChem.EmbedMolecule(molecule, randomSeed=7)
    assert audit_hch_angles(molecule, molecule.GetConformer().GetPositions())["angles"] == []


def test_reconstruction_keeps_heavy_coordinates_and_original_atom_count():
    molecule = Chem.AddHs(Chem.MolFromSmiles("CC"))
    AllChem.EmbedMolecule(molecule, randomSeed=7)
    positions = np.asarray(molecule.GetConformer().GetPositions())
    neighbors = [atom.GetIdx() for atom in molecule.GetAtomWithIdx(0).GetNeighbors() if atom.GetAtomicNum() == 1]
    positions[neighbors[1]] = positions[neighbors[0]]+np.array([.01, .01, .01])
    repaired, changed = rebuild_flagged_carbon_hydrogens(molecule, positions)
    assert set(changed) == set(neighbors)
    assert np.array_equal(repaired[:2], positions[:2])
    assert repaired.shape == positions.shape
    assert not audit_hch_angles(molecule, repaired)["n_gross_distortions"]


def test_untriggered_reconstruction_is_exact_no_op():
    molecule = Chem.AddHs(Chem.MolFromSmiles("C"))
    AllChem.EmbedMolecule(molecule, randomSeed=7)
    positions = np.asarray(molecule.GetConformer().GetPositions())
    repaired, changed = rebuild_flagged_carbon_hydrogens(molecule, positions)
    assert not changed
    assert np.array_equal(positions, repaired)
