import importlib.util
from pathlib import Path
import sys

import pytest

from fmotpiefinal.core.structure import Atom, Ligand, Protein


sys.path.insert(0, str(Path(__file__).parents[1] / "tools"))
SPEC = importlib.util.spec_from_file_location("water_reference", Path(__file__).parents[1] / "tools/run_qm_water_reference.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_intact_water_selection_has_no_covalent_cut():
    ligand = Ligand("L", (Atom(0, "C", (0, 0, 0)),))
    atoms = (Atom(1, "O", (3, 0, 0), "A", "HOH", "1"),
             Atom(2, "H", (3.96, 0, 0), "A", "HOH", "1"),
             Atom(3, "H", (2.76, .93, 0), "A", "HOH", "1"))
    key, distance, selected = MODULE.closest_complete_water(Protein(atoms), ligand)
    assert key == "A:HOH:1"
    assert distance == 3
    assert len(selected) == 3
    with pytest.raises(ValueError, match="incomplete water"):
        MODULE.closest_complete_water(Protein(atoms[:2]), ligand)


def test_missing_water_not_synthetically_invented():
    with pytest.raises(ValueError, match="no retained"):
        MODULE.closest_complete_water(Protein(()), Ligand("L", (Atom(0, "C", (0, 0, 0)),)))


def test_optional_correlated_reference_smoke(tmp_path):
    pytest.importorskip("pyscf")
    atoms = [Atom(0, "H", (0, 0, 0)), Atom(1, "H", (0, 0, .74))]
    result = MODULE.reference_energy(atoms, 0, "def2-svpd", set(), tmp_path / "h2", 1)
    assert result["scf_converged"]
    assert result["mp2_hartree"] < result["hf_hartree"]
