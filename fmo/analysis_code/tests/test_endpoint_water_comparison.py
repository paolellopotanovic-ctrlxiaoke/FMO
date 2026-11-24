from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from compare_endpoint_water_control import validate_water_only_change


DECK = """ $FMO NFRAG=2
 FRGNAM(1)=AHOH604, LIG
 INDAT(1)=0 1 -3 0 4 -5 0
 $END
 $FMOXYZ
 O 8 1 0 0
 H 1 2 0 0
 H 1 1 1 0
 N 7 5 0 0
 H 1 6 0 0
 $END
"""


def test_only_water_coordinates_can_change():
    updated = DECK.replace("O 8 1 0 0", "O 8 2 0 0")
    assert validate_water_only_change(DECK, updated) == 3


def test_ligand_change_rejected():
    with pytest.raises(ValueError, match="protein, ligand or caps"):
        validate_water_only_change(DECK, DECK.replace("N 7 5 0 0", "N 7 7 0 0"))


def test_qm_setting_change_rejected():
    with pytest.raises(ValueError, match="controls"):
        validate_water_only_change(DECK, DECK.replace("NFRAG=2", "NFRAG=3"))
