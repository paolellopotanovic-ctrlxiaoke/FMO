import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from fmotpiefinal.core.conjugated_amine_geometry import conjugated_nh2_groups, group_departure, induced_conjugated_nh2, reconstruct_flagged_nh2


def fixture_geometry():
    molecule = Chem.AddHs(Chem.MolFromSmiles("NC(=[NH+](C))N(C)C"))
    assert AllChem.EmbedMolecule(molecule, randomSeed=20261025) == 0
    points = molecule.GetConformer().GetPositions().copy()
    group = conjugated_nh2_groups(molecule)[0]
    points[group["center"]] = [0., 0., 0.]
    points[group["nitrogen"]] = [1.35, 0., 0.]
    points[group["other_nitrogens"][0]] = [-.675, 1.17, 0.]
    points[group["other_nitrogens"][1]] = [-.675, -1.17, 0.]
    points[group["hydrogens"][0]] = [1.85, np.sqrt(.75), 0.]
    points[group["hydrogens"][1]] = [1.85, -np.sqrt(.75), 0.]
    return molecule, points, group


def test_amide_aniline_and_ammonium_are_not_misclassified_as_guanidino_nh2():
    for smiles in ("CC(=O)N", "Nc1ccccc1", "C[NH3+]", "CO"):
        assert not conjugated_nh2_groups(Chem.AddHs(Chem.MolFromSmiles(smiles)))


def test_induced_deformation_requires_good_input_and_large_new_departure():
    molecule, original, group = fixture_geometry()
    distorted = original.copy()
    distorted[group["hydrogens"][0]] = [1.85, 0., np.sqrt(.75)]
    diagnostics = induced_conjugated_nh2(molecule, original, distorted)
    assert diagnostics[0]["new_deformation_flag"]
    assert not induced_conjugated_nh2(molecule, distorted, distorted)[0]["new_deformation_flag"]


def test_repair_changes_only_flagged_hydrogens_and_preserves_bond_lengths():
    molecule, original, group = fixture_geometry()
    distorted = original.copy()
    distorted[group["hydrogens"][0]] = [1.85, 0., np.sqrt(.75)]
    diagnostics = induced_conjugated_nh2(molecule, original, distorted)
    repaired = reconstruct_flagged_nh2(molecule, distorted, diagnostics)
    assert group_departure(repaired, group) < 1.e-6
    unchanged = [index for index in range(molecule.GetNumAtoms()) if index not in group["hydrogens"]]
    assert np.array_equal(repaired[unchanged], distorted[unchanged])
    for index in group["hydrogens"]:
        assert np.linalg.norm(repaired[index]-repaired[group["nitrogen"]]) == pytest.approx(np.linalg.norm(distorted[index]-distorted[group["nitrogen"]]))


def test_nonfinite_coordinates_fail_closed():
    molecule, points, _ = fixture_geometry()
    points[0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        induced_conjugated_nh2(molecule, points, points)
