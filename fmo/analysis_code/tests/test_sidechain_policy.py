from dataclasses import replace

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

from fmotpiefinal.core.sidechain_policy import SidechainPolicy, amide_changes, select_sidechains
from fmotpiefinal.core.structure import Atom, Protein


def receptor_atom(name="CB", residue="LYS", xyz=(2.0, 0.0, 0.0), number="1", element="C"):
    return Atom(0, element, xyz, "A", residue, number, "ATOM", name)


def methane():
    molecule = Chem.AddHs(Chem.MolFromSmiles("C"))
    AllChem.EmbedMolecule(molecule, randomSeed=2026)
    conformer = molecule.GetConformer()
    offset = np.asarray(conformer.GetAtomPosition(0))
    for index in range(molecule.GetNumAtoms()):
        conformer.SetAtomPosition(index, np.asarray(conformer.GetAtomPosition(index))-offset)
    return molecule


def test_overlap_selects_sidechain_without_using_labels():
    molecule = methane()
    molecule.SetProp("affinity", "100")
    protein = Protein((receptor_atom(),))
    plan = select_sidechains(protein, {"lig": molecule}, {"lig": molecule.GetConformer().GetPositions()})
    molecule.SetProp("affinity", "-100")
    second = select_sidechains(protein, {"lig": molecule}, {"lig": molecule.GetConformer().GetPositions()})
    assert plan == second
    assert plan["selected_residues"] == [["A", "1", "LYS"]]
    assert not plan["experimental_labels_used"]


def test_no_trigger_keeps_baseline_even_inside_contact_radius():
    molecule = methane()
    protein = Protein((receptor_atom(xyz=(3.8, 0, 0)),))
    plan = select_sidechains(protein, {"lig": molecule}, {"lig": molecule.GetConformer().GetPositions()})
    assert plan["selected_residues"] == []
    assert plan["action"] == "KEEP_BASELINE_UNCHANGED"
    assert len(plan["edges"]) == 2


def test_backbone_contact_does_not_select_sidechain():
    molecule = methane()
    protein = Protein((receptor_atom(name="CA"),))
    plan = select_sidechains(protein, {"lig": molecule}, {"lig": molecule.GetConformer().GetPositions()})
    assert not plan["selected_residues"]


def test_coordination_and_disulfides_are_protected():
    molecule = methane()
    bound = {"lig": molecule.GetConformer().GetPositions()}
    nitrogen = receptor_atom(name="NZ", element="N")
    zinc = Atom(1, "ZN", (2.0, 0, 2.0), "A", "ZN", "50", "HETATM", "ZN")
    plan = select_sidechains(Protein((nitrogen, zinc)), {"lig": molecule}, bound)
    assert not plan["selected_residues"]
    assert "METAL_COORDINATION_CANDIDATE" in plan["nodes"][0]["protection"]
    sulfur = receptor_atom(name="SG", residue="CYS", element="S")
    partner = replace(sulfur, index=1, residue_number="2", xyz=(4.1, 0, 0))
    plan = select_sidechains(Protein((sulfur, partner)), {"lig": molecule}, bound)
    assert not plan["selected_residues"]
    assert any(edge["kind"] == "disulfide_distance_candidate" for edge in plan["edges"])


def test_polar_contact_not_called_clash_at_normal_screen_distance():
    molecule = Chem.AddHs(Chem.MolFromSmiles("O"))
    AllChem.EmbedMolecule(molecule, randomSeed=4)
    conformer = molecule.GetConformer()
    origin = np.asarray(conformer.GetAtomPosition(0))
    for index in range(molecule.GetNumAtoms()):
        conformer.SetAtomPosition(index, np.asarray(conformer.GetAtomPosition(index))-origin)
    protein = Protein((receptor_atom(name="NZ", element="N", xyz=(2.4, 0, 0)),))
    plan = select_sidechains(protein, {"lig": molecule}, {"lig": conformer.GetPositions()})
    assert not plan["selected_residues"]


def test_new_amide_deformation_is_detected_but_existing_twist_is_not_new():
    molecule = Chem.AddHs(Chem.MolFromSmiles("CC(=O)NC"))
    AllChem.EmbedMolecule(molecule, randomSeed=4)
    conformer = molecule.GetConformer()
    rdMolTransforms.SetDihedralDeg(conformer, 2, 1, 3, 4, 180.0)
    original = conformer.GetPositions().copy()
    rdMolTransforms.SetDihedralDeg(conformer, 2, 1, 3, 4, 120.0)
    bound = conformer.GetPositions().copy()
    assert any(row["new_deformation_flag"] for row in amide_changes(molecule, original, bound, SidechainPolicy()))
    assert not any(row["new_deformation_flag"] for row in amide_changes(molecule, bound, bound, SidechainPolicy()))


def test_proline_and_nonpeptide_covalent_links_are_never_freed():
    molecule = methane()
    bound = {"lig": molecule.GetConformer().GetPositions()}
    first = receptor_atom()
    second = receptor_atom(number="2", residue="ASN")
    plan = select_sidechains(Protein((first, second)), {"lig": molecule}, bound,
        covalent_links=[(("A", "1", "LYS"), ("A", "2", "ASN"))])
    assert not plan["selected_residues"]
    proline = replace(first, residue_name="PRO")
    assert not select_sidechains(Protein((proline,)), {"lig": molecule}, bound)["selected_residues"]


def test_invalid_inputs_fail_explicitly():
    with pytest.raises(ValueError):
        SidechainPolicy(heavy_overlap_angstrom=float("nan"))
    with pytest.raises(ValueError, match="identities"):
        select_sidechains(Protein(()), {"lig": methane()}, {})
    molecule = methane()
    with pytest.raises(ValueError, match="count/order"):
        select_sidechains(Protein(()), {"lig": molecule}, {"lig": np.zeros((1, 3))})
