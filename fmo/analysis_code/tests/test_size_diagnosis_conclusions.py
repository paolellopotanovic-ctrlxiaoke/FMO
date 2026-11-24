"""Local geometry conclusions retain provenance and are not corrected scores."""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
SPEC = importlib.util.spec_from_file_location("conclude_size_diagnosis", ROOT / "tools/conclude_size_diagnosis.py")
CONCLUSIONS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CONCLUSIONS)


def test_kabsch_corrects_rotation_without_reflection():
    first = np.array([[0., 0., 0.], [1., 0., 0.], [0., 2., 0.], [0., 0., 3.]])
    rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    second = first @ rotation + np.array([3., 4., 5.])
    fitted_rotation, translation = CONCLUSIONS.kabsch(first, second)
    assert np.allclose(first @ fitted_rotation + translation, second)
    assert np.linalg.det(fitted_rotation) == pytest.approx(1)


def test_donor_hydrogen_acceptor_angle():
    assert CONCLUSIONS.angle(np.array([0., 0., 0.]), np.array([1., 0., 0.]), np.array([3., 0., 0.])) == pytest.approx(180)


def test_conclusions_do_not_claim_subset_as_fixed_performance():
    result = json.loads((ROOT / "docs/SIZE_DIAGNOSIS_CONCLUSIONS.json").read_text())
    assert not result["new_QM"]
    assert not result["production_changed"]
    assert len(result["OpenBind_charge_audit"]) == 6
    assert sum(bool(row["charged_amide_N_matches"]) for row in result["OpenBind_charge_audit"]) == 2
    assert result["OpenBind_308_318_receptor_alignment"]["n_CA"] == 139
