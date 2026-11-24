"""A matching input keyword does not guarantee matching effective solvation."""
import pytest

from fmotpiefinal.core.solution_cycle import HARTREE_TO_KCAL_MOL, pcm_energy_inventory, require_consistent_pcm_cycle


def log(dispersion=0., repulsion=0., disabled=False):
    total = -10. + 20. + dispersion + repulsion
    return f"""
 EXECUTION OF GAMESS TERMINATED NORMALLY
 ELECTROSTATIC INTERACTION = -10.000 KCAL/MOL
 PIEROTTI CAVITATION ENERGY = 20.000 KCAL/MOL
 DISPERSION FREE ENERGY = {dispersion:.3f} KCAL/MOL
 REPULSION FREE ENERGY = {repulsion:.3f} KCAL/MOL
 TOTAL INTERACTION = {total:.3f} KCAL/MOL
 Free unco+D energy in solvent= {-100. + total / HARTREE_TO_KCAL_MOL:.10f}
 Internal unco+D energy in solvent= -100.0000000000
 """ + ("Setting IDISP=0 for NCAV2= 3" if disabled else "")


def test_pcm_inventory_reconciles_energy_and_components():
    inventory = pcm_energy_inventory(log(-5., 1.))
    assert inventory["total"] == pytest.approx(6.)


def test_automatic_multicavity_disable_rejected_even_with_same_keywords():
    states = ("complex", "receptor", "ligand")
    decks = {state: "$PCM SOLVNT=WATER IEF=-10 ICAV=1 IDISP=1 IFMO=-1 $END" for state in states}
    logs = {state: log(disabled=state != "ligand") for state in states}
    with pytest.raises(ValueError, match="automatically disabled"):
        require_consistent_pcm_cycle(decks, logs)


def test_explicit_idisp0_cycle_is_supported():
    states = ("complex", "receptor", "ligand")
    decks = {state: "$PCM SOLVNT=WATER IEF=-10 ICAV=1 IDISP=0 IFMO=-1 $END" for state in states}
    assert require_consistent_pcm_cycle(decks, {state: log() for state in states})["consistent"]


def test_missing_or_nonclosing_energy_rejected():
    with pytest.raises(ValueError, match="completed"):
        pcm_energy_inventory("")
    with pytest.raises(ValueError, match="components"):
        pcm_energy_inventory(log().replace("TOTAL INTERACTION = 10.000", "TOTAL INTERACTION = 1.000"))


def test_same_pcm_but_different_qm_dispersion_rejected():
    states = ("complex", "receptor", "ligand")
    decks = {state: "$PCM SOLVNT=WATER IEF=-10 ICAV=1 IDISP=0 IFMO=-1 $END" for state in states}
    decks["ligand"] += " $DFT DC=.TRUE. IDCVER=4 $END"
    with pytest.raises(ValueError, match="electronic"):
        require_consistent_pcm_cycle(decks, {state: log() for state in states})
