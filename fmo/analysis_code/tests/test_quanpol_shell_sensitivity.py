import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

from compare_quanpol_shell_sensitivity import summarize


def make_run(root, shell, corrections, *, changed_ligand=False):
    rows = []
    for leg in ("bound", "solution"):
        for endpoint in (0, 10):
            frame = root / f"{leg}_{endpoint}"
            frame.mkdir(parents=True)
            coordinates = np.asarray([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
            if changed_ligand and endpoint == 10:
                coordinates[1, 0] += 0.01
            np.savez(frame / "physical_coordinates.npz", ligand_xyz=coordinates,
                     environment_xyz=np.asarray([[2.0, 0.0, 0.0], [3.0, 0.0, 0.0]]),
                     environment_original_indices=np.asarray([5, 6]))
            (frame / "input_manifest.json").write_text(json.dumps({
                "trajectory": "same.nc", "replica": 2, "qm_atoms": 2, "mm_atoms": 10,
                "residues": [
                    {"name": "HOH", "indices": [1, 2, 3]},
                    {"name": "NA", "indices": [4]},
                    {"name": "ALA", "indices": [5, 6]},
                ],
            }))
            delta = corrections[(leg, endpoint)]
            rows.append({
                "repeat": 0, "iteration": 100, "leg": leg, "endpoint": endpoint,
                "trajectory": "same.nc", "replica": 2,
                "qm": {"status": "SUCCEEDED", "scf_unconverged": False,
                       "log": str(frame / "qm.stdout")},
                "mm": {"status": "SUCCEEDED", "scf_unconverged": False},
                "mm_parity_difference_kcal_mol": 0.01,
                "delta_u_kcal_mol": delta,
            })
    status = root / "status.json"
    status.write_text(json.dumps({"status": "COMPLETE", "labels_used": False,
                                 "target": "CDK8", "shell_angstrom": shell,
                                 "swrb_angstrom": 6000.0,
                                 "qm_center_restraint_allowed": False, "results": rows}))
    return status


def corrections(offset):
    return {("bound", 0): offset, ("bound", 10): offset + 3,
            ("solution", 0): offset + 1, ("solution", 10): offset + 2}


def test_shell_compare_holds_frame_fixed_and_reports_cycle(tmp_path):
    small = make_run(tmp_path / "small", 8.0, corrections(0.0))
    large = make_run(tmp_path / "large", 16.0, corrections(10.0))

    report = summarize([(8.0, small), (16.0, large)])

    assert report["fixed_four_state_cycles"] == 1
    assert [row["correction_kcal_mol"] for row in report["cycles"][0]["shell_corrections"]] == pytest.approx([2, 2])
    assert report["cycle_summary_by_shell_kcal_mol"]["8.0"]["mean_kcal_mol"] == pytest.approx(2)
    assert report["frame_states"][0]["shell_results"][0]["water_residues"] == 1
    assert report["frame_states"][0]["shell_results"][0]["ion_residues"] == 1


def test_shell_compare_rejects_changed_ligand_geometry(tmp_path):
    small = make_run(tmp_path / "small", 8.0, corrections(0.0))
    large = make_run(tmp_path / "large", 16.0, corrections(0.0), changed_ligand=True)

    with pytest.raises(ValueError, match="changed the sampled frame or QM ligand geometry"):
        summarize([(8.0, small), (16.0, large)])
