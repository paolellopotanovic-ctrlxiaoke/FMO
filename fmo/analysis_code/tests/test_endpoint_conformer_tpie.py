from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from run_endpoint_conformer_tpie import conformer_deck
from run_cdk8_solution_cycle import group


SOURCE = """ $FMO NFRAG=2 NBODY=2 MODMOL=1 $END
 $PCM IDISP=1 $END
 $FMOXYZ
 O 8 9 9 9
 C 6 0 0 0
 N 7 1 0 0
 O 8 0 1 0
 H 1 0 0 1
 $END
"""


def test_rigid_motion_removed_and_receptor_model_preserved():
    original = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    moved = original @ np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]) + 7
    deck, audit = conformer_deck(SOURCE, ["C", "N", "O", "H"], moved)
    assert deck[:deck.index("$FMOXYZ")] == SOURCE[:SOURCE.index("$FMOXYZ")]
    rows = group(deck, "FMOXYZ")[1].strip().splitlines()
    assert rows[0] == "O 8 9 9 9"
    assert np.allclose(np.array([row.split()[2:] for row in rows[1:]], float), original)
    assert audit["ligand_internal_rmsd_angstrom"] < 1e-10


def test_atom_order_and_nonfinite_coordinates_rejected():
    with pytest.raises(ValueError, match="inventory"):
        conformer_deck(SOURCE, ["N", "C", "O", "H"], np.zeros((4, 3)))
    with pytest.raises(ValueError, match="finite"):
        conformer_deck(SOURCE, ["C", "N", "O", "H"], np.full((4, 3), np.nan))
