"""Frozen deck and contact diagnostics reject broken provenance."""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
SPEC = importlib.util.spec_from_file_location("dissect_size_interactions", ROOT / "tools/dissect_size_interactions.py")
DISSECT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DISSECT)


def deck(indices="0 1 -2 0 3 0"):
    return f""" $FMO NFRAG=2
 FRGNAM(1)=ARES1, LIG
 INDAT(1)={indices}
 $END
 $FMOXYZ
 C 6 0.0 0.0 0.0
 H 1 1.0 0.0 0.0
 O 8 3.0 0.0 0.0
 $END
 """


def test_fragment_range_mapping_is_one_based_and_complete():
    names, groups, numbers, coordinates = DISSECT.parse_deck(deck())
    assert names == ["ARES1", "LIG"]
    assert groups == [[0, 1], [2]]
    assert numbers.tolist() == [6, 1, 8]
    assert coordinates.shape == (3, 3)


def test_missing_or_duplicated_fragment_atom_rejected():
    with pytest.raises(ValueError, match="incomplete"):
        DISSECT.parse_deck(deck("0 1 -2 0 2 0"))
    with pytest.raises(ValueError, match="unterminated"):
        DISSECT.parse_deck(deck("0 1 -2 0 3"))


def test_remote_receptor_does_not_bury_ligand_surface():
    free, bound = DISSECT.ligand_surface(np.array([[0., 0., 0.]]), np.array([6]),
        np.array([[100., 0., 0.]]), np.array([6]))
    assert free[0] > 0
    assert np.allclose(free, bound)


def test_near_receptor_buries_some_but_not_all_surface():
    free, bound = DISSECT.ligand_surface(np.array([[0., 0., 0.]]), np.array([6]),
        np.array([[3.5, 0., 0.]]), np.array([6]))
    assert 0 < bound[0] < free[0]


def test_saved_report_keeps_complete_cohorts_and_closing_pair_deltas():
    result = json.loads((ROOT / "docs/SIZE_STRUCTURE_INTERACTIONS.json").read_text())
    assert result["production_changed"] is False
    assert result["new_QM"] is False
    assert sum(group["n"] for group in result["groups"].values()) == 629
    for group in result["groups"].values():
        assert len(group["records"]) == group["n"]
        for pair in group["representative_pairs"]:
            delta = pair["delta_components_energy_second_minus_first"]
            assert delta["tpie"] == pytest.approx(pair["sum_all_fragment_delta"], abs=1e-7)
            assert sum(value for field, value in delta.items() if field != "tpie") == pytest.approx(delta["tpie"], abs=.1)
    assert sum(record["n_crystals"] for record in result["groups"]["OpenBind"]["records"].values()) == 639
