import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

from decompose_quanpol_local_correction import HARTREE_KCAL_MOL, decompose, isolated_deck, parse_vacuum


def test_isolated_deck_preserves_qm_and_removes_environment():
    text = " $CONTRL RUNTYP=ENERGY $END\n $DFTB NDFTB=3 $END\n $QUANPO NSTEP=0\n SWRB=6000 $END\n $DATA\nQM geometry\n $END\n $FFDATA\nMM parameters\n $END\n"
    vacuum = isolated_deck(text)
    assert "QM geometry" in vacuum
    assert "RUNTYP=ENERGY" in vacuum
    assert "$QUANPO" not in vacuum
    assert "$FFDATA" not in vacuum
    assert "MM parameters" not in vacuum
    with pytest.raises(ValueError):
        isolated_deck(" $DATA\n $END")


def test_vacuum_parser_requires_success_and_unique_energy():
    text = "FINAL RDFTB ENERGY IS -52.0 AFTER 10 ITERATIONS\nEXECUTION OF GAMESS TERMINATED NORMALLY"
    assert parse_vacuum(text) == pytest.approx(-52 * HARTREE_KCAL_MOL)
    with pytest.raises(ValueError):
        parse_vacuum(text + "\nSCF IS UNCONVERGED")
    with pytest.raises(ValueError):
        parse_vacuum(text + "\nFINAL RDFTB ENERGY IS -53.0")


def test_four_state_decomposition_closure_and_missing_point():
    rows = []
    for repeat in (0, 1):
        for leg in ("bound", "solution"):
            for endpoint in (0, 10):
                internal = endpoint * (2 if leg == "bound" else 1) + repeat
                environment = endpoint * (3 if leg == "bound" else 1) - repeat
                rows.append({"repeat": repeat, "iteration": 400, "leg": leg, "endpoint": endpoint,
                             "total_correction_kcal_mol": internal + environment,
                             "internal_correction_kcal_mol": internal,
                             "environment_correction_kcal_mol": environment})
    report = decompose(rows)
    assert report["cycle_statistics"]["total"]["mean"] == 30
    assert report["cycle_statistics"]["internal"]["mean"] == 10
    assert report["cycle_statistics"]["environment"]["mean"] == 20
    assert report["max_abs_closure_kcal_mol"] == 0
    with pytest.raises(ValueError, match="incomplete"):
        decompose(rows[:-1])
    with pytest.raises(ValueError, match="duplicate"):
        decompose(rows + rows[:1])
    rows[0]["total_correction_kcal_mol"] += 1
    with pytest.raises(ValueError, match="does not close"):
        decompose(rows)


def test_rejects_nonfinite_energy_before_statistics():
    row = {"repeat": 0, "iteration": 400, "leg": "bound", "endpoint": 0,
           "total_correction_kcal_mol": float("nan"),
           "internal_correction_kcal_mol": 0, "environment_correction_kcal_mol": 0}
    with pytest.raises(ValueError, match="nonfinite"):
        decompose([row])
