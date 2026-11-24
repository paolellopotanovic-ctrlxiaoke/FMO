from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

import report_weak_tpie_contacts as report
from fmotpiefinal.core.parser import PairInteraction


def test_fragment_charge_inventory_is_explicit_and_complete(tmp_path):
    deck = tmp_path / "input.inp"
    deck.write_text("ICHARG(1)=-1, 0, 1\n FRGNAM(1)=AASP1, AASP2, LIG\n INDAT(1)=0 1 0")
    assert report.fragment_charges(deck) == {"AASP1": -1, "AASP2": 0, "LIG": 1}
    deck.write_text("ICHARG(1)=-1, 0\n FRGNAM(1)=AASP1, AASP2, LIG\n INDAT(1)=0 1 0")
    with pytest.raises(ValueError, match="counts"):
        report.fragment_charges(deck)


def test_named_pairs_reproduces_full_sum_and_handles_both_index_orders(tmp_path, monkeypatch):
    deck = tmp_path / "input.inp"
    deck.write_text("NFRAG=3 FRGNAM(1)=AASP1, AHOH2, LIG\n INDAT(1)=0 1 0")
    (tmp_path / "run").mkdir()
    (tmp_path / "run/gamess.log").write_text("fixture")
    monkeypatch.setattr(report, "parse_pie", lambda path: [
        PairInteraction(1, 3, -12.), PairInteraction(3, 2, 2.)])
    values, fingerprint = report.named_pairs(deck, -10.)
    assert values == {"AASP1": -12., "AHOH2": 2.}
    assert len(fingerprint) == 64
    with pytest.raises(ValueError, match="did not reproduce"):
        report.named_pairs(deck, -11.)


def test_missing_or_duplicate_fragment_inventory_is_not_reported(tmp_path, monkeypatch):
    deck = tmp_path / "input.inp"
    deck.write_text("NFRAG=3 FRGNAM(1)=AASP1, AASP1, LIG\n INDAT(1)=0 1 0")
    with pytest.raises(ValueError, match="inventory"):
        report.named_pairs(deck, -10.)
    deck.write_text("NFRAG=3 FRGNAM(1)=AASP1, AHOH2, LIG\n INDAT(1)=0 1 0")
    monkeypatch.setattr(report, "parse_pie", lambda path: [PairInteraction(1, 3, -10.)])
    with pytest.raises(ValueError, match="did not reproduce"):
        report.named_pairs(deck, -10.)
