"""Thermodynamic-cycle deck extraction keeps exact source chemistry."""
import importlib.util
from pathlib import Path
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
SPEC = importlib.util.spec_from_file_location("solution_cycle", ROOT / "tools/run_cdk8_solution_cycle.py")
CYCLE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CYCLE)


DECK = """! source
 $CONTRL RUNTYP=ENERGY $END
 $DFTB SCC=.TRUE. DFTB3=.TRUE. $END
 $FMOPRP IPIEDA=1 $END
 $PCM SOLVNT=WATER $END
 $PCMCAV RADII=VANDW $END
 $TESCAV NTSALL=60 $END
 $FMO NFRAG=3 NBODY=2 NACUT=0 MODMOL=1 MOLFRG(1)=3
 ICHARG(1)=1, -1, 0
 FRGNAM(1)=ONE, TWO, LIG
 INDAT(1)=0 1 -2 0 3 -4 0 5 -6 0
 $END
 $FMOHYB
 HOP_C 4 4
 $END
 $FMOBND
 -1 3 HOP_C
 $END
 $FMOXYZ
 C 6 0 0 0
 H 1 1 0 0
 N 7 2 0 0
 H 1 3 0 0
 O 8 4 0 0
 H 1 5 0 0
 $END
 $DATA
 test
 C1
 $END
"""


def test_split_preserves_receptor_and_bound_coordinates(tmp_path):
    decks = CYCLE.split_decks(DECK, tmp_path)
    assert "MODMOL" not in decks["complex"]
    assert "IPIEDA=0" in decks["complex"]
    assert "NFRAG=2 NBODY=2" in decks["receptor"]
    assert "NFRAG=1 NBODY=1" in decks["ligand"]
    assert "O 8 4 0 0" not in decks["receptor"]
    assert "C 6 0 0 0" not in decks["ligand"]
    assert "$FMOBND" not in decks["ligand"]
    assert "-1 3 HOP_C" in decks["receptor"]
    assert "ICHARG(1)=0" in decks["ligand"]
    assert "INDAT(1)=0 1 -2 0" in decks["ligand"]


def test_covalent_ligand_is_not_silently_split(tmp_path):
    with pytest.raises(ValueError, match="covalently"):
        CYCLE.split_decks(DECK.replace("-1 3 HOP_C", "-1 5 HOP_C"), tmp_path)


def test_gas_keeps_qm_model_and_removes_pcm():
    gas = CYCLE.gas_deck(DECK + " $DFT DC=.TRUE. IDCVER=4 $END\n")
    assert "$DFTB " in gas
    assert "$DFT DC=.TRUE. IDCVER=4" in gas
    assert "$PCM " not in gas
    assert "$PCMCAV " not in gas
    assert "$TESCAV " not in gas


def test_wrapped_charge_array_keeps_commas():
    rendered = CYCLE.format_array("ICHARG", [0] * 80)
    lines = rendered.splitlines()
    assert len(lines) > 1
    assert all(line.endswith(",") for line in lines[:-1])
    assert not lines[-1].endswith(",")


def test_energy_direction_and_no_fitted_coefficients():
    result = CYCLE.metrics([-9., -8., -7.], [-30., -20., -10.])
    assert result["spearman"] == pytest.approx(1)
    assert len(CYCLE.ARMS) == 7


def test_free_minimum_strain_has_positive_binding_cost():
    from evaluate_cdk8_solution_cycle import scores
    values = scores(-50., -60., -30., 8.)
    assert values["three_state_pcm_free_min"] == pytest.approx(-22.)
    assert values["tpie_plus_ligand_strain"] == pytest.approx(-42.)
    assert values["tpie_plus_cycle_solvation"] == pytest.approx(-20.)


def test_numerical_recovery_preserves_tolerances_and_physics():
    from recover_cdk8_solution_cycle import recovery_deck
    source = "MAXIT=200 DAMPEX=4.00 $END MCONV(2)=530 MCONV(4)=530 CONV=1.0E-8"
    recovered = recovery_deck(source, False)
    assert "CONV=1.0E-8" in recovered
    assert "MAXIT=200" in recovered
    assert "ITYPMX=-1" in recovered
    assert recovery_deck("NSTEP=200 OPTTOL=0.000100", True) == "NSTEP=1000 OPTTOL=0.000100"
