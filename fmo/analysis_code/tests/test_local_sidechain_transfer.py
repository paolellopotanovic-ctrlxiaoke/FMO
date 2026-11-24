from pathlib import Path
import sys

import pytest


sys.path.insert(0, str(Path(__file__).parents[1] / "tools"))
from run_local_sidechain_transfer import directory_key, evaluate, topology_signature


def test_output_keys_preserve_unique_registered_names():
    assert directory_key("17h flip") == "17h_flip"
    assert directory_key("21 (Xtal)") == "21__Xtal_"


def test_topology_comparison_ignores_energy_controls_but_not_charge():
    reference = "$FMO NFRAG=1 ICHARG(1)=0 FRGNAM(1)=LIG INDAT(1)=0 1 -3 0 $END"
    assert topology_signature(reference) == topology_signature(reference + "\n $PCM IDISP=0 $END")
    assert topology_signature(reference) != topology_signature(reference.replace("ICHARG(1)=0", "ICHARG(1)=1"))


def test_incomplete_cohort_never_gets_subset_metrics(tmp_path):
    workspace = tmp_path / "workspace"
    (workspace / "data").mkdir(parents=True)
    (workspace / "data/fmo87_cohort.csv").write_text("target,ligand_name,experimental_dg_kcal_mol\n")
    output = tmp_path / "output"
    output.mkdir()
    report = evaluate(workspace, tmp_path, {"CHK1": {"example": {"status": "SUCCEEDED", "tpie": -10}}}, output)
    cohort = report["groups"]["CHK1"]
    assert cohort["n_expected"] == 13
    assert cohort["n_succeeded"] == 1
    assert "local_sidechains" not in cohort
    assert "baseline" not in cohort


def test_failure_remains_explicit(tmp_path):
    (tmp_path / "data").mkdir()
    (tmp_path / "data/fmo87_cohort.csv").write_text("target,ligand_name,experimental_dg_kcal_mol\n")
    report = evaluate(tmp_path, tmp_path, {"FXa": {"5f": {"status": "FAILED", "error": "SCF"}}}, tmp_path)
    assert report["groups"]["FXa"]["failures"]["5f"]["error"] == "SCF"
    assert report["groups"]["FXa"]["n_succeeded"] == 0
