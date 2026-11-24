import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from fmotpiefinal.core.sidechain_motion import add_displaced_clash_seeds
from fmotpiefinal.core.sidechain_policy import select_sidechains
from fmotpiefinal.core.structure import Atom, Protein


def inputs():
    molecule = Chem.AddHs(Chem.MolFromSmiles("C"))
    AllChem.EmbedMolecule(molecule, randomSeed=7)
    positions = molecule.GetConformer().GetPositions()
    positions -= positions[0]
    for index, position in enumerate(positions):
        molecule.GetConformer().SetAtomPosition(index, position)
    protein = Protein((Atom(0, "C", (2., 0., 0.), "A", "LYS", "1", "ATOM", "CB"),))
    bound = {"lig": positions + np.array([-2., 0., 0.])}
    molecules = {"lig": molecule}
    return protein, molecules, bound


def test_original_overlap_resolved_by_motion_is_a_new_seed():
    protein, molecules, bound = inputs()
    seed = select_sidechains(protein, molecules, bound)
    assert not seed["selected_residues"]
    result = add_displaced_clash_seeds(protein, molecules, bound, seed)
    assert result["selected_residues"] == [["A", "1", "LYS"]]
    assert result["motion_events"][0]["ligand_displacement_angstrom"] == 2.
    assert not seed["selected_residues"]


def test_motion_seed_preserves_special_network_protection():
    protein, molecules, bound = inputs()
    seed = select_sidechains(protein, molecules, bound)
    seed["nodes"][0]["protection"] = ["METAL_COORDINATION_CANDIDATE"]
    result = add_displaced_clash_seeds(protein, molecules, bound, seed)
    assert not result["selected_residues"]
    assert result["nodes"][0]["action"] == "REVIEW_PROTECTED_NETWORK_KEEP_FIXED"


def test_small_motion_alone_does_not_trigger():
    protein, molecules, bound = inputs()
    seed = select_sidechains(protein, molecules, bound)
    bound["lig"] = molecules["lig"].GetConformer().GetPositions() + np.array([-.1, 0., 0.])
    result = add_displaced_clash_seeds(protein, molecules, bound, seed)
    assert not result["selected_residues"]


def test_invalid_motion_threshold_and_identity_fail_explicitly():
    protein, molecules, bound = inputs()
    seed = select_sidechains(protein, molecules, bound)
    with pytest.raises(ValueError):
        add_displaced_clash_seeds(protein, molecules, bound, seed, float("nan"))
    with pytest.raises(ValueError, match="identities"):
        add_displaced_clash_seeds(protein, molecules, {}, seed)
