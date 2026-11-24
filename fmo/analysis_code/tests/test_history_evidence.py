"""Historical claims stay auditable without a local research workspace."""
import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("audit_history_evidence", ROOT / "tools/audit_history_evidence.py")
HISTORY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(HISTORY)


def test_portable_evidence_snapshot():
    result = HISTORY.verify(ROOT)
    assert result["status"] == "PASS"
    assert result["source_count"] == 40
    assert result["workspace_sources_verified"] is False


def test_evidence_excerpts_use_original_line_numbers():
    result = HISTORY.excerpts("first\n## matched\nthird\n", ["## matched"])
    assert result == [{"anchor": "## matched", "start_line": 2,
                       "end_line": 3, "text": "## matched\nthird"}]


def test_missing_evidence_anchor_is_not_silently_skipped():
    with pytest.raises(ValueError, match="anchor missing"):
        HISTORY.excerpts("different report", ["## required"])


def test_changed_metric_snapshot_is_rejected(tmp_path):
    destination = tmp_path / "docs"
    destination.mkdir()
    (destination / "HISTORY_EVIDENCE.json").write_bytes((ROOT / "docs/HISTORY_EVIDENCE.json").read_bytes())
    (destination / "HISTORY_METRICS.json").write_text("{}")
    with pytest.raises(ValueError, match="metric snapshot changed"):
        HISTORY.verify(tmp_path)


def test_changed_original_source_is_rejected(tmp_path):
    root = tmp_path / "release"
    (root / "docs").mkdir(parents=True)
    (root / "docs/metrics.json").write_text("{}")
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    source = workspace / "source.md"
    source.write_text("original")
    evidence = {"metrics_file": "docs/metrics.json", "metrics_sha256": HISTORY.digest(root / "docs/metrics.json"),
                "sources": {"H01": {"workspace_relative_path": "source.md", "sha256": HISTORY.digest(source)}}}
    (root / "docs/HISTORY_EVIDENCE.json").write_text(json.dumps(evidence))
    assert HISTORY.verify(root, workspace)["workspace_sources_verified"] is True
    source.write_text("modified")
    with pytest.raises(ValueError, match="source fingerprint changed: H01"):
        HISTORY.verify(root, workspace)


def test_latest_report_table_matches_frozen_metrics():
    metrics = json.loads((ROOT / "docs/HISTORY_METRICS.json").read_text())
    report = (ROOT / "docs/IMPROVEMENT_HISTORY.zh-CN.md").read_text()
    groups = metrics["latest_2004_groups"]
    assert len(groups) == 14
    assert sum(group.get("n_crystals", group.get("n")) for group in groups.values()) == 909
    assert sum(group.get("n_compounds", group.get("n")) for group in groups.values()) == 760
    for name, group in groups.items():
        count = group.get("n_compounds", group.get("n"))
        values = (group["raw"]["pearson"], group["softcap"]["pearson"],
                  group["raw"]["spearman"], group["raw"]["kendall"])
        formatted = " | ".join(f"{value:.3f}".replace("-", "−") for value in values)
        assert f"| {name} | {count} | {formatted} |" in report
        assert group["raw"]["spearman"] == group["softcap"]["spearman"]
        assert group["raw"]["kendall"] == group["softcap"]["kendall"]
    status = metrics["latest_2004_final_status"]
    assert status["status"] == "COMPLETE_INDEPENDENTLY_VERIFIED"
    assert status["n_requested"] == status["n_verified"] == 909


def test_k10_report_table_matches_paired_evidence():
    metrics = json.loads((ROOT / "docs/HISTORY_METRICS.json").read_text())
    report = (ROOT / "docs/IMPROVEMENT_HISTORY.zh-CN.md").read_text()
    for setting in ("k100", "k25", "k10"):
        result = metrics["restraint_dose"]["metrics"][f"tpie_{setting}"]
        assert f"| wet、{setting} | {result['pearson']:.6f} | {result['kendall']:.6f} |" in report
    for score, values in metrics["latest_original87_weighted"].items():
        assert f"{values['pearson']:.6f}" in report
        for metric in ("spearman", "kendall"):
            assert f"{values[metric]:.6f}" in report
