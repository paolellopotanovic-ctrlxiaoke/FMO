import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from run_dense_endpoint_complex_tpie import reuse_result
from run_matched_openfe import digest


def make_task(tmp_path, directory):
    path = tmp_path / directory
    path.mkdir()
    (path / "input.inp").write_text("fixed deck")
    source = tmp_path / "source.inp"
    source.write_text("fixed source")
    task = {"target": "T", "name": "A", "trajectory": "trajectory", "iteration": 400,
            "replica": 2, "endpoint": 0, "directory": str(path), "source": str(source),
            "input_sha256": digest(path / "input.inp"), "source_sha256": digest(source),
            "runtime": {"version": "same"}, "parameters": "same", "ligand_fragment": 2}
    return task


def test_exact_reuse_preserves_result_and_provenance(tmp_path):
    cached = {**make_task(tmp_path, "old"), "result": {"status": "SUCCEEDED", "tpie_kcal_mol": -3}}
    Path(cached["directory"], "summary.json").write_text(json.dumps(cached))
    task = make_task(tmp_path, "new")
    row = reuse_result(task, cached)
    assert row["result"] == cached["result"]
    assert row["directory"] == task["directory"]
    assert row["reused_summary_sha256"]


def test_changed_input_cannot_be_reused(tmp_path):
    cached = {**make_task(tmp_path, "old"), "result": {"status": "SUCCEEDED"}}
    task = make_task(tmp_path, "new")
    Path(cached["directory"], "input.inp").write_text("changed")
    with pytest.raises(ValueError, match="digest"):
        reuse_result(task, cached)


def test_changed_chemistry_cannot_be_reused(tmp_path):
    cached = {**make_task(tmp_path, "old"), "result": {"status": "SUCCEEDED"}}
    task = {**make_task(tmp_path, "new"), "ligand_fragment": 3}
    with pytest.raises(ValueError, match="chemistry"):
        reuse_result(task, cached)
