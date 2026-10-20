import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from fmotpiefinal.scoring import REFERENCE_SCHEMA, score_run


def create_run(tmp_path, *, rescue=False):
    root = tmp_path / "run"
    (root / "geometry").mkdir(parents=True)
    (root / "geometry/geometry_gate.json").write_text(json.dumps({
        "schema": "fmo-new-target-geometry-gate-v1", "target": "Example", "n_ligands": 3}))
    if rescue:
        source = root / "rescues/registered/blind_predictions.json"
    else:
        source = root / "qm/full_target_fmo/blind_predictions.json"
    source.parent.mkdir(parents=True)
    source.write_text(json.dumps({"n_requested": 3, "n_succeeded": 3,
                                  "predictions": [{"ligand": name, "status": "SUCCEEDED", "tpie": value}
                                                  for name, value in (("a", -10), ("b", -20), ("c", -30))]}))
    status = {"ranking_eligible": True, "target": "Example", "n_ligands": 3,
              "status": "RESCUED_ALTERNATE_PROTOCOL" if rescue else
                        "BLIND_FMO_COMPLETE_REQUIRES_EXTERNAL_VALIDATION",
              "blind_predictions_sha256": hashlib.sha256(source.read_bytes()).hexdigest()}
    if rescue:
        status.update(selected_arm="registered", arms=[{"arm": "registered", "complete": True}])
    (root / "adaptive_status.json").write_text(json.dumps(status))
    return root, source


@pytest.mark.parametrize("rescue", [False, True])
def test_completed_run_preserves_rank_and_is_repeatable(tmp_path, rescue):
    root, _ = create_run(tmp_path, rescue=rescue)
    result = score_run(root)
    assert result["scoring"]["mode"] == "cohort_local_descriptive_only"
    assert [row["ligand"] for row in result["ranked"]] == ["c", "b", "a"]
    assert [row["softcap_score"] for row in result["ranked"]] == sorted(
        [row["softcap_score"] for row in result["ranked"]], reverse=True)
    assert (root / "ranking/local_reference.json").is_file()
    assert (root / "ranking/pocket_sensitivity.json").is_file()
    assert (root / "ranking/preparation_sensitivity.json").is_file()
    preparation = json.loads((root / "ranking/preparation_sensitivity.json").read_text())
    assert preparation["ranking_modified"] is False
    assert preparation["experimental_labels_used"] is False
    assert score_run(root) == result


def test_frozen_reference_does_not_use_local_median(tmp_path):
    root, _ = create_run(tmp_path)
    reference = tmp_path / "reference.json"
    reference.write_text(json.dumps({"schema": REFERENCE_SCHEMA, "factor": 3.0,
                                     "mad_conversion": 1.4826,
                                     "reference": {"Example": {"median": 100.0,
                                                               "scaled_mad": 20.0}}}))
    result = score_run(root, reference)
    assert result["scoring"]["mode"] == "previously_frozen_target_reference"
    assert result["ranked"][0]["softcap_score"] > 0
    from fmotpiefinal.preparation_sensitivity import audit_selected_preparation
    before = (root / "ranking/ranking.json").read_bytes()
    assert audit_selected_preparation(root, result)["ranking_modified"] is False
    assert (root / "ranking/ranking.json").read_bytes() == before
    with pytest.raises(ValueError, match="frozen ranking differs"):
        score_run(root)


def test_single_crystal_registered_rescue_fingerprint(tmp_path):
    root, source = create_run(tmp_path, rescue=True)
    status_path = root / "adaptive_status.json"
    status = json.loads(status_path.read_text())
    status.pop("blind_predictions_sha256")
    status["arms"] = [{"arm": "registered", "tpie": -20,
                       "prediction_sha256": hashlib.sha256(source.read_bytes()).hexdigest()}]
    status_path.write_text(json.dumps(status))
    assert score_run(root)["n"] == 3


def test_builtin_reference_is_used_for_developed_target(tmp_path):
    root, _ = create_run(tmp_path)
    gate_path = root / "geometry/geometry_gate.json"
    gate = json.loads(gate_path.read_text())
    gate["target"] = "Galectin3"
    gate_path.write_text(json.dumps(gate))
    status_path = root / "adaptive_status.json"
    status = json.loads(status_path.read_text())
    status["target"] = "Galectin3"
    status_path.write_text(json.dumps(status))
    result = score_run(root)
    assert result["scoring"]["mode"] == "bundled_retrospective_target_reference"
    assert result["scoring"]["reference_development_set_retrospective"] is True
    assert result["scoring"]["reference_sha256"]


@pytest.mark.parametrize("failure", ["fingerprint", "coverage", "duplicate", "review", "bad_arm"])
def test_incomplete_or_inconsistent_run_fails_closed(tmp_path, failure):
    root, source = create_run(tmp_path, rescue=failure == "bad_arm")
    status_path = root / "adaptive_status.json"
    status = json.loads(status_path.read_text())
    if failure == "fingerprint":
        source.write_text(source.read_text() + "\n")
    elif failure in ("coverage", "duplicate"):
        data = json.loads(source.read_text())
        if failure == "coverage":
            data["n_succeeded"] = 2
        else:
            data["predictions"][1]["ligand"] = "a"
        source.write_text(json.dumps(data))
        status["blind_predictions_sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
    elif failure == "review":
        status["ranking_eligible"] = False
    else:
        status["selected_arm"] = "../other"
    status_path.write_text(json.dumps(status))
    with pytest.raises(ValueError):
        score_run(root)
    assert not (root / "ranking/ranking.json").exists()


def test_unrecognized_target_reference_fails_closed(tmp_path):
    root, _ = create_run(tmp_path)
    reference = tmp_path / "reference.json"
    reference.write_text(json.dumps({"schema": REFERENCE_SCHEMA, "factor": 3.0,
                                     "mad_conversion": 1.4826, "reference": {}}))
    with pytest.raises(ValueError, match="no compatible frozen reference"):
        score_run(root, reference)


def test_follow_on_cohort_uses_prior_unlabeled_reference(tmp_path):
    earlier, _ = create_run(tmp_path / "earlier")
    score_run(earlier)
    later, _ = create_run(tmp_path / "later")
    result = score_run(later, earlier / "ranking/local_reference.json")
    assert result["scoring"]["mode"] == "previously_frozen_target_reference"
    assert not (later / "ranking/local_reference.json").exists()


def test_zero_mad_keeps_physical_score_without_inventing_scale(tmp_path):
    root, source = create_run(tmp_path)
    data = json.loads(source.read_text())
    for entry in data["predictions"]:
        entry["tpie"] = -10
    source.write_text(json.dumps(data))
    status_path = root / "adaptive_status.json"
    status = json.loads(status_path.read_text())
    status["blind_predictions_sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
    status_path.write_text(json.dumps(status))
    report = score_run(root)
    assert not report["softcap_defined"]
    assert all(row["softcap_score"] is None for row in report["ranked"])


def test_installed_module_cli_scores_validated_run(tmp_path):
    root, _ = create_run(tmp_path)
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(Path(sys.modules[score_run.__module__].__file__).resolve().parents[1])
    completed = subprocess.run([sys.executable, "-m", "fmotpiefinal.pipeline.cli",
                                "score", "--run-dir", str(root)], text=True, capture_output=True,
                               env=environment)
    assert completed.returncode == 0, completed.stderr
    assert json.loads(completed.stdout)["n"] == 3
    assert (root / "ranking/ranking.csv").is_file()


def test_label_free_cohort_diagnostics_do_not_modify_rank(tmp_path):
    root, _ = create_run(tmp_path)
    (root / "intake").mkdir()
    (root / "intake/intake.json").write_text(json.dumps({"ligands": [
        {"ligand_id": name, "canonical_smiles": smiles, "formal_charge": 0}
        for name, smiles in (("a", "C"), ("b", "CC"), ("c", "CCC"))]}))
    result = score_run(root)
    diagnostic = json.loads((root / "ranking/diagnostics.json").read_text())
    assert diagnostic["experimental_labels_used"] is False
    assert diagnostic["ranking_modified"] is False
    assert diagnostic["heavy_atom_count_min_max"] == [1, 3]
    assert diagnostic["score_heavy_atom_pearson"] == 1.0
    assert [row["ligand"] for row in result["ranked"]] == ["c", "b", "a"]
    assert score_run(root) == result


def test_diagnostics_fail_closed_for_changed_intake_identities(tmp_path):
    root, _ = create_run(tmp_path)
    (root / "intake").mkdir()
    (root / "intake/intake.json").write_text(json.dumps({"ligands": [
        {"ligand_id": name, "canonical_smiles": "C", "formal_charge": 0}
        for name in ("a", "b", "wrong")]}))
    with pytest.raises(ValueError, match="intake ligand identities"):
        score_run(root)


def test_pocket_sensitivity_audit_is_separate_from_frozen_scores(tmp_path):
    root, _ = create_run(tmp_path)
    pocket = root / "qm/pocket_selection_audit.json"
    pocket.parent.mkdir(exist_ok=True)
    pocket.write_text(json.dumps({"n_ligands": 3, "ligand_names": ["b", "a", "c"],
                                  "n_residues_all_atom": 12, "n_residues_heavy_atom": 10,
                                  "membership_differences": [{"residue": "A:HOH:12"}],
                                  "selection_production": "all_atom"}))
    original = score_run(root)
    report = json.loads((root / "ranking/pocket_sensitivity.json").read_text())
    assert report["hydrogen_sensitive_residues"] == 1
    assert report["experimental_labels_used"] is False
    assert score_run(root) == original
    pocket.write_text(pocket.read_text().replace('"b", "a", "c"', '"a", "a", "c"'))
    with pytest.raises(ValueError, match="pocket membership audit"):
        score_run(root)
