import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from run_pocket_boundary_control import expansion_keys, physical_content
from fmotpiefinal.core.structure import Atom, Ligand, Protein


def test_expansion_preserves_only_registered_waters():
    atoms = (Atom(0, "C", (1., 0., 0.), "A", "ALA", "1", atom_name="CA"),
             Atom(1, "C", (6., 0., 0.), "A", "ALA", "2", atom_name="CA"),
             Atom(2, "O", (2., 0., 0.), "W", "HOH", "3", "HETATM", "O"),
             Atom(3, "O", (3., 0., 0.), "W", "HOH", "4", "HETATM", "O"))
    ligand = Ligand("test", (Atom(0, "C", (0., 0., 0.)),))
    keys = expansion_keys(Protein(atoms), [ligand], {"A:ALA:1", "W:HOH:3"}, 7.)
    assert keys == {"A:ALA:1", "A:ALA:2", "W:HOH:3"}


def test_expansion_rejects_missing_residue_or_nonexpansion():
    protein = Protein((Atom(0, "C", (0., 0., 0.), "A", "ALA", "1"),))
    with pytest.raises(ValueError, match="must exceed"):
        expansion_keys(protein, [], set(), 5.)
    with pytest.raises(ValueError, match="missing"):
        expansion_keys(protein, [], {"A:ALA:99"}, 7.)


def test_physical_fingerprint_requires_complete_controls():
    with pytest.raises(ValueError, match="expected one"):
        physical_content(" $FMO NFRAG=2 $END")
