import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"tools"))
from analyze_ligand_strain_components import PAIR_FIELDS, strain_components
from fmotpiefinal.core.parser import PairInteraction


def test_component_names_match_the_actual_pieda_parser():
    pair = PairInteraction(1, 2, 0.)
    assert all(hasattr(pair, field+"_kcal_mol") for field in PAIR_FIELDS)


def test_strain_decomposition_keeps_solvent_and_internal_response_separate():
    components = strain_components({"internal_hartree": 0., "total": 3.},
                                   {"internal_hartree": 0., "total": 1.}, 2.)
    assert components["internal_response"] == 0.
    assert components["solvent_response"] == 2.
    assert components["closure_error"] == 0.


def test_inconsistent_energy_inventory_is_rejected():
    with pytest.raises(ValueError, match="do not close"):
        strain_components({"internal_hartree": 0., "total": 3.},
                          {"internal_hartree": 0., "total": 1.}, 5.)
