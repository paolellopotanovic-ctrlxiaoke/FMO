from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

from run_ligand_pose_control import minimization_arguments


def test_controls_change_only_one_geometry_constraint():
    fixed = minimization_arguments("input_heavy_fixed")
    unrestrained = minimization_arguments("hydrogens_unrestrained")
    assert fixed == {"restraint_k_kj_mol_a2": 10.0, "move_protein": False, "fix_ligand_heavy": True}
    assert unrestrained == {"restraint_k_kj_mol_a2": 10.0, "move_protein": False, "restrain_ligand_hydrogens": False}
    with pytest.raises(ValueError):
        minimization_arguments("pick_best")
