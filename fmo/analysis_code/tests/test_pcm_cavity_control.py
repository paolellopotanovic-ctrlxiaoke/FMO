import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from run_pcm_cavity_control import documented_cavity
from run_cdk8_solution_cycle import group


def test_documented_recipe_changes_only_cavity_and_nonpolar_switches():
    deck = " $PCM SOLVNT=WATER IEF=-10 ICAV=1 IDISP=1 IFMO=-1 $END\n $PCMCAV RADII=VANDW $END\n $FMOXYZ C 6 1 2 3 $END"
    changed = documented_cavity(deck)
    assert "RADII=DFTB33OB ALPHA(1)=1.0" in changed
    assert "ICAV=0 IDISP=0" in changed
    assert group(changed, "FMOXYZ")[1] == group(deck, "FMOXYZ")[1]
    assert "IFMO=-1" in changed


def test_custom_radii_are_not_overwritten():
    deck = " $PCM ICAV=1 IDISP=1 $END\n $PCMCAV RADII=SUAHF $END"
    with pytest.raises(ValueError, match="custom radii"):
        documented_cavity(deck)
