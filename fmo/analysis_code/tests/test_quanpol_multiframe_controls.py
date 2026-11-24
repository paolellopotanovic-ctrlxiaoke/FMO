import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

from compare_quanpol_multiframe_controls import build_comparison
from quanpol_constants import DEFAULT_SWRB_ANGSTROM, validate_swrb


def test_safe_default_disables_qm_center_restraint():
    assert DEFAULT_SWRB_ANGSTROM == 6000.0
    validate_swrb(DEFAULT_SWRB_ANGSTROM)
    with pytest.raises(ValueError, match="QM-center restraint"):
        validate_swrb(1001.0)
    validate_swrb(1001.0, allow_qm_center_restraint=True)


def make_run(tmp_path, name, swrb, corrections, center_values):
    results = []
    endpoints = (("bound", 0), ("bound", 10), ("solution", 0), ("solution", 10))
    for (leg, endpoint), delta_u, center in zip(endpoints, corrections, center_values, strict=True):
        frame_dir = tmp_path / name / "repeat_0" / "iteration_0400" / leg / f"endpoint_{endpoint:02d}"
        frame_dir.mkdir(parents=True, exist_ok=True)
        log = frame_dir / "qm_matched.stdout"
        log.write_text(f" QM CENTER                     ENERGY = {center:.8f}  KCAL/MOL\n")
        (frame_dir / "mm_matched.stdout").write_text("matched MM\n")
        (frame_dir / "physical_coordinates.npz").write_bytes(b"same coordinates")
        (frame_dir / "environment.pdb").write_text("same environment\n")
        stem = f"r0_i400_{leg}_{endpoint}"
        (frame_dir / f"qm_{stem}.inp").write_text(
            f"SWRB={swrb}\nCOORDINATES\nC 6.0 0.0 0.0 0.0\nSTOP\nQM input\n")
        (frame_dir / f"mm_{stem}.inp").write_text("matched MM deck\n")
        results.append({
            "repeat": 0,
            "iteration": 400,
            "leg": leg,
            "endpoint": endpoint,
            "trajectory": "fixed trajectory",
            "replica": 3 if leg == "bound" else 10,
            "delta_u_kcal_mol": delta_u,
            "qm": {"log": str(log), "potential_kcal_mol": delta_u, "status": "SUCCEEDED",
                   "scf_unconverged": False},
            "mm": {"log": str(frame_dir / "mm_matched.stdout"), "potential_kcal_mol": 0.0,
                   "status": "SUCCEEDED", "scf_unconverged": False},
            "openmm_ffdata": {"potential_kcal_mol": 0.0},
        })
    path = tmp_path / f"{name}.json"
    path.write_text(json.dumps({"status": "COMPLETE", "target": "CDK8", "labels_used": False,
                                "swrb_angstrom": swrb, "times_iterations": [400], "results": results}))
    return path


def test_center_penalty_change_closes_against_energy_component(tmp_path):
    with_center = make_run(tmp_path, "with", 1001, [1.0, 3.0, 2.0, 4.0], [0.0, 0.0, 10.0, 0.0])
    without_center = make_run(tmp_path, "without", 6000, [1.0, 3.0, -8.0, 4.0], [0.0, 0.0, 0.0, 0.0])

    report = build_comparison(with_center, without_center)

    assert report["qm_center_cycle_mean_kcal_mol"] == pytest.approx(10.0)
    assert report["max_abs_closure_error_kcal_mol"] == pytest.approx(0.0)
    assert report["mm_max_abs_change_kcal_mol"] == pytest.approx(0.0)


def test_control_comparison_requires_paired_frames(tmp_path):
    with_center = make_run(tmp_path, "with", 1001, [1.0, 3.0, 2.0, 4.0], [0.0] * 4)
    without_center = make_run(tmp_path, "without", 6000, [1.0, 3.0, 2.0, 4.0], [0.0] * 4)
    data = json.loads(without_center.read_text())
    data["results"].pop()
    without_center.write_text(json.dumps(data))

    with pytest.raises(ValueError, match="incomplete"):
        build_comparison(with_center, without_center)


def test_control_comparison_rejects_non_swrb_input_changes(tmp_path):
    with_center = make_run(tmp_path, "with", 1001, [1.0, 3.0, 2.0, 4.0], [0.0] * 4)
    without_center = make_run(tmp_path, "without", 6000, [1.0, 3.0, 2.0, 4.0], [0.0] * 4)
    qm_deck = (tmp_path / "without" / "repeat_0" / "iteration_0400" / "bound"
               / "endpoint_00" / "qm_r0_i400_bound_0.inp")
    qm_deck.write_text(qm_deck.read_text().replace("QM input", "changed QM model"))

    with pytest.raises(ValueError, match="non-SWRB input"):
        build_comparison(with_center, without_center)
