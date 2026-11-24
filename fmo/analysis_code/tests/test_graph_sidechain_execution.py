from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).parents[1] / "tools"))
from run_graph_sidechain_transfer import execute


def test_no_trigger_reuses_baseline_explicitly_without_running_qm(tmp_path, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("unchanged baseline must not start minimization or QM")

    monkeypatch.setattr("run_graph_sidechain_transfer.prepare_job", forbidden)
    monkeypatch.setattr("run_graph_sidechain_transfer.score_job", forbidden)
    cohort = {"target": "example", "local_residues": [], "decision_plan_sha256": "plan",
        "baseline": {"lig": {"tpie": -42.0, "input_sha256": "input", "log_sha256": "log", "n_pairs": 2}}}
    target, name, result = execute((cohort, "lig", str(tmp_path)))
    assert (target, name) == ("example", "lig")
    assert result["tpie"] == result["baseline_tpie"] == -42.0
    assert result["execution_path"] == "UNCHANGED_VERIFIED_BASELINE_NO_STRUCTURAL_TRIGGER"
    assert result["new_qm_executed"] is False
    assert (tmp_path / "example/lig/summary.json").is_file()


def test_failed_preparation_is_kept_and_never_scored(tmp_path, monkeypatch):
    def prepare(payload):
        return "example", "lig", {"status": "FAILED", "error": "missing model"}

    def forbidden(*args, **kwargs):
        raise AssertionError("failed preparation must not reach QM")

    monkeypatch.setattr("run_graph_sidechain_transfer.prepare_job", prepare)
    monkeypatch.setattr("run_graph_sidechain_transfer.score_job", forbidden)
    cohort = {"target": "example", "local_residues": [["A", "1", "LYS"]], "decision_plan_sha256": "plan"}
    _, _, result = execute((cohort, "lig", str(tmp_path)))
    assert result["status"] == "FAILED"
    assert result["error"] == "missing model"
    assert result["execution_path"] == "GRAPH_SELECTED_PREPARATION_FAILED"


def test_complete_evidence_distinguishes_new_qm_from_unchanged_reuse():
    import json
    evidence = Path(__file__).parents[1] / "docs/GRAPH_SIDECHAIN_RESULTS.zh-CN.json"
    report = json.loads(evidence.read_text())
    assert report["final_status"]["status"] == "COMPLETE"
    assert report["final_status"]["n_succeeded"] == report["final_status"]["n_requested"] == 270
    assert report["total_new_qm"] == 10
    assert report["total_unchanged_reused"] == 260
    assert not report["production_default_changed"]
    for target, group in report["groups"].items():
        assert not group["failures"]
        if target != "CDK8":
            assert group["baseline"] == group["graph_guided"]
            assert group["n_new_qm"] == 0
        else:
            assert group["n_new_qm"] == 10
            assert len(group["independent_geometry_audit"]) == 10
