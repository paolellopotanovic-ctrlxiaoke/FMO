import importlib.util
from pathlib import Path
import sys

import pytest


sys.path.insert(0, str(Path(__file__).parents[1] / "tools"))
SPEC = importlib.util.spec_from_file_location("partial_openfe", Path(__file__).parents[1] / "tools/report_partial_openfe.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_partial_diagnostic_is_not_complete_free_energy(tmp_path):
    path = tmp_path / "analysis.yaml"
    path.write_text("- iteration: 400\n  percent_complete: 20.0\n  mbar_analysis:\n"
        "    free_energy_in_kT: 4.0\n    standard_error_in_kT: 0.5\n"
        "    number_of_uncorrelated_samples: 100\n    n_equilibrium_iterations: 10\n")
    row = MODULE.leg_snapshot(path, .6, 1.0)
    assert row["binding_B_minus_A_kcal_mol"] == pytest.approx(1.4)
    assert row["bound_internal_mbar_se_kcal_mol"] == pytest.approx(.3)
    assert row["production_complete"] is False
    assert row["full_uncertainty_estimated"] is False
