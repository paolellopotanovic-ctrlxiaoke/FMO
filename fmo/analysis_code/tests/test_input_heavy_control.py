import sys
from dataclasses import replace
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from run_input_heavy_control import heavy_displacement
from fmotpiefinal.core.structure import Atom, Ligand


def test_fixed_heavy_geometry_ignores_relaxed_hydrogen():
    ligand = Ligand("test", (Atom(0, "C", (0., 0., 0.)), Atom(1, "H", (1., 0., 0.))))
    optimized = replace(ligand, atoms=(ligand.atoms[0], replace(ligand.atoms[1], xyz=(0., 1., 0.))))
    assert heavy_displacement(ligand, optimized) == 0.


def test_changed_heavy_geometry_is_rejected():
    ligand = Ligand("test", (Atom(0, "C", (0., 0., 0.)),))
    with pytest.raises(ValueError, match="heavy atoms moved"):
        heavy_displacement(ligand, replace(ligand, atoms=(replace(ligand.atoms[0], xyz=(.01, 0., 0.)),)))


def test_changed_elements_are_not_matched_by_position():
    ligand = Ligand("test", (Atom(0, "C", (0., 0., 0.)),))
    with pytest.raises(ValueError, match="identity mismatch"):
        heavy_displacement(ligand, replace(ligand, atoms=(replace(ligand.atoms[0], element="N"),)))
