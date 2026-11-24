from pathlib import Path

import pytest

from fmotpiefinal.core.structure import Atom, Protein
from fmotpiefinal.pipeline.target_stages import receptor_at_minimized_geometry


def atom(index, name, residue, element):
    return Atom(index=index, atom_name=name, residue_name=residue, element=element,
                residue_number="1", chain="A", xyz=(1.0, 2.0, 3.0), record="ATOM")


def test_verified_openmm_terminal_cap_aliases_preserve_source_identity(monkeypatch):
    original = Protein((atom(0, "HH31", "ACE", "H"), atom(1, "CH3", "NME", "C")))
    minimized = Protein((atom(0, "H1", "ACE", "H"), atom(1, "C", "NME", "C")))
    monkeypatch.setattr("fmotpiefinal.pipeline.target_stages.read_pdb", lambda path: minimized)
    assert receptor_at_minimized_geometry(original, Path("unused.pdb")) == original


def test_noncap_atom_renaming_is_not_accepted(monkeypatch):
    original = Protein((atom(0, "HH31", "LYS", "H"),))
    minimized = Protein((atom(0, "H1", "LYS", "H"),))
    monkeypatch.setattr("fmotpiefinal.pipeline.target_stages.read_pdb", lambda path: minimized)
    with pytest.raises(ValueError, match="identity/order"):
        receptor_at_minimized_geometry(original, Path("unused.pdb"))


def test_cap_element_changes_are_not_accepted(monkeypatch):
    original = Protein((atom(0, "CH3", "NME", "C"),))
    minimized = Protein((atom(0, "C", "NME", "N"),))
    monkeypatch.setattr("fmotpiefinal.pipeline.target_stages.read_pdb", lambda path: minimized)
    with pytest.raises(ValueError, match="identity/order"):
        receptor_at_minimized_geometry(original, Path("unused.pdb"))
