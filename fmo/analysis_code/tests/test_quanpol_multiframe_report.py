import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

import pytest

from report_quanpol_multiframe import render_markdown, single_state_weight_ess, summarize


def run_report(tmp_path, *, parity=0.01, omit=None, labels_used=False):
    results = []
    for repeat in (0, 1):
        for iteration in (400, 1200, 2000):
            for leg in ("bound", "solution"):
                for endpoint in (0, 10):
                    key = (repeat, iteration, leg, endpoint)
                    if key == omit:
                        continue
                    correction = endpoint + (repeat * 2) + (iteration / 400)
                    if leg == "solution":
                        correction = 0.7 * endpoint + (repeat * 2) + (iteration / 400) - 3
                    log_path = tmp_path / f"r{repeat}_i{iteration}_{leg}_{endpoint}.out"
                    log_path.write_text("QM CENTER                     ENERGY = 0.00000000 KCAL/MOL\n")
                    results.append({
                        "repeat": repeat,
                        "iteration": iteration,
                        "leg": leg,
                        "endpoint": endpoint,
                        "delta_u_kcal_mol": correction,
                        "mm_parity_difference_kcal_mol": parity,
                        "qm": {"status": "SUCCEEDED", "scf_unconverged": False,
                               "log": str(log_path),
                               "potential_kcal_mol": correction},
                        "mm": {"status": "SUCCEEDED", "scf_unconverged": False,
                               "potential_kcal_mol": 0.0},
                    })
    path = tmp_path / "run.json"
    import json
    path.write_text(json.dumps({
        "status": "COMPLETE",
        "target": "CDK8",
        "labels_used": labels_used,
        "times_iterations": [400, 1200, 2000],
        "driver_sha256": "frozen",
        "results": results,
    }))
    return path


def test_summarizes_complete_frames_and_checks_mm_parity(tmp_path):
    report = summarize([run_report(tmp_path)])

    assert report["paired_single_points"] == 24
    assert report["fixed_frame_cycles"] == 6
    assert report["cycle_correction_kcal_mol"]["mean"] == pytest.approx(3)
    assert report["all_matched_mm_parity_passed_0p1_kcal_mol"]
    assert report["max_abs_breakdown_closure_kcal_mol"] == pytest.approx(0)
    ess = report["single_state_zwanzig_weight_ess_diagnostic"]["bound_0"]["effective_sample_size"]
    assert 1 <= ess <= 6
    markdown = render_markdown(report)
    assert "### 单向权重诊断（不是自由能结果）" in markdown
    assert "每状态 6 个" in markdown
    assert "周期 PME" in markdown
    assert markdown.index("| iteration 400 |") < markdown.index("### 单向权重诊断")
    assert "ESS 只是稀疏帧" not in markdown


def test_single_state_weight_ess_is_shift_invariant():
    ess, _, free_energy = single_state_weight_ess([-3.0, -2.0, -1.0])
    shifted, _, shifted_free_energy = single_state_weight_ess([97.0, 98.0, 99.0])
    assert shifted == pytest.approx(ess)
    assert shifted_free_energy == pytest.approx(free_energy + 100)


def test_rejects_incomplete_or_label_selected_frames(tmp_path):
    with pytest.raises(ValueError, match="incomplete frame matrix"):
        summarize([run_report(tmp_path, omit=(0, 400, "bound", 0))])

    other = tmp_path / "labelled"
    other.mkdir()
    with pytest.raises(ValueError, match="label-dependent"):
        summarize([run_report(other, labels_used=True)])


def test_rejects_failed_qm_or_bad_mm_parity(tmp_path):
    with pytest.raises(ValueError, match="MM parity"):
        summarize([run_report(tmp_path, parity=0.2)])
