import sys
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from audit_cycle_ligand_chemistry import geometry_screen, stereochemistry


def test_stereochemical_inversion_is_not_hidden_by_preserved_atom_order():
    molecule = Chem.AddHs(Chem.MolFromSmiles("F[C@H](Cl)Br"))
    assert AllChem.EmbedMolecule(molecule, randomSeed=20261024) == 0
    points = molecule.GetConformer().GetPositions()
    reference = stereochemistry(molecule, points)
    assert reference
    assert geometry_screen(molecule, points, reference)["passed_screen"]
    mirrored = points*np.asarray([-1., 1., 1.])
    assert geometry_screen(molecule, mirrored, reference)["stereo_changed"]


def test_stretched_hydrogen_bond_is_flagged():
    molecule = Chem.AddHs(Chem.MolFromSmiles("CO"))
    assert AllChem.EmbedMolecule(molecule, randomSeed=20261024) == 0
    points = molecule.GetConformer().GetPositions()
    reference = stereochemistry(molecule, points)
    points[-1] += np.asarray([10., 0., 0.])
    assert geometry_screen(molecule, points, reference)["bond_distance_flags"]
