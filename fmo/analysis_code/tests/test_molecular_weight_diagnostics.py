"""Descriptive size diagnostics must not become fitted production scores."""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
SPEC = importlib.util.spec_from_file_location("diagnose_molecular_weight", ROOT / "tools/diagnose_molecular_weight.py")
DIAGNOSTICS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DIAGNOSTICS)


def test_linear_mass_residualization_removes_shared_trend():
    mass = np.arange(1.0, 7.0)
    pattern = np.array([1., -2., 1., 1., -2., 1.])
    observed = mass * 3 + pattern
    predicted = mass * 7 - pattern
    assert DIAGNOSTICS.partial(observed, predicted, mass) == pytest.approx(-1)


def test_near_mass_pairs_report_ties_and_dependent_pair_counts():
    labels = np.array([1., 2., 3.])
    mass = np.array([100., 100., 130.])
    result = DIAGNOSTICS.pair_order(labels, mass, mass, 10)
    assert result["pairs"] == 1
    assert result["ties"] == 1
    assert result["accuracy_excluding_ties"] is None


def test_scaffold_centering_removes_between_group_offsets():
    centered = DIAGNOSTICS.center([1., 3., 10., 14.], ["a", "a", "b", "b"])
    assert centered.tolist() == [-1., 1., -2., 2.]


def test_bootstrap_reproducible():
    values = np.arange(10.)
    assert DIAGNOSTICS.bootstrap(values, values, 100) == pytest.approx([1., 1.])
