import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"tools"))
from run_hydrogen_normalized_strain import evaluate, fixed_heavy_displacement, hydrogen_only_deck, transplant_ligand
from recover_hydrogen_normalized_strain import qa_recovery_deck


DECK = """ $CONTRL RUNTYP=ENERGY SCFTYP=RHF $END
 $FMO NFRAG=1 NBODY=1 NACUT=0 $END
 $FMOXYZ
C 6 0 0 0
H 1 1 0 0
O 8 2 0 0
H 1 3 0 0
 $END
"""


def test_only_actual_hydrogens_are_active_in_cartesian_optimizer():
    deck = hydrogen_only_deck(DECK)
    assert "IACTAT(1)=2,4" in deck and "NZVAR=0" in deck
    assert "MODGRD=32" in deck and "OPTTOL=0.000100" in deck
    assert "METHOD=RFO NSTEP=500" in deck
    assert all(len(line) <= 80 for line in deck.splitlines())


def test_unsupported_or_double_constraints_are_not_silently_replaced():
    with pytest.raises(ValueError):
        hydrogen_only_deck(DECK.replace("RUNTYP=ENERGY", "RUNTYP=OPTIMIZE"))
    with pytest.raises(ValueError):
        hydrogen_only_deck(DECK+" $STATPT IFREEZ(1)=1 $END\n")


def test_fixed_heavy_geometry_and_element_order_are_enforced():
    before = ["C 6 0 0 0", "H 1 1 0 0"]
    assert fixed_heavy_displacement(before, [before[0], "H 1 1.1 0 0"]) == 0.
    with pytest.raises(ValueError, match="heavy atoms moved"):
        fixed_heavy_displacement(before, ["C 6 .1 0 0", before[1]])
    with pytest.raises(ValueError, match="elements changed"):
        fixed_heavy_displacement(before, ["O 8 0 0 0", before[1]])


def test_transplant_cannot_replace_an_unmatched_ligand():
    with pytest.raises(ValueError, match="exact complex suffix"):
        transplant_ligand(DECK.replace("C 6 0 0 0", "C 6 .1 0 0"), DECK,
                          ["C 6 0 0 0", "H 1 1 0 0", "O 8 2 0 0", "H 1 3 0 0"])


def test_failed_case_prevents_complete_cohort_metrics():
    result = evaluate({"one": {"status": "FAILED"}}, {"names": ["one"]}, {}, {"one": 1.})
    assert result["status"] == "INCOMPLETE" and not result["metrics"]


def test_qa_rescue_changes_only_numerical_optimizer_controls():
    original = hydrogen_only_deck(DECK)
    rescued = qa_recovery_deck(original)
    assert "OPTTOL=0.000100" in rescued and "IACTAT(1)=2,4" in rescued
    assert rescued.replace("METHOD=QA", "METHOD=RFO").replace("TRMIN=0.000010", "TRMIN=0.001") == original
    with pytest.raises(ValueError):
        qa_recovery_deck(DECK)
