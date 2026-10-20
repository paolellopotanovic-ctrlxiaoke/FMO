"""Preflight diagnoses missing dependencies without starting a calculation."""
import pytest

from fmotpiefinal import doctor


def test_dependency_failure_is_a_structured_action(monkeypatch):
    def missing(name):
        raise ImportError("dependency missing")
    monkeypatch.setattr(doctor.importlib, "import_module", missing)
    monkeypatch.setattr(doctor.shutil, "which", lambda name: None)
    report = doctor.inspect_environment()
    assert report["status"] == "ACTION_REQUIRED"
    assert not report["qm_paths_checked"]
    assert all(not check["ok"] for check in report["checks"])


def test_partial_qm_configuration_cannot_pass(tmp_path, monkeypatch):
    monkeypatch.setattr(doctor.importlib, "import_module", lambda name: object())
    monkeypatch.setattr(doctor.shutil, "which", lambda name: "/bin/pdb2pqr")
    report = doctor.inspect_environment(gamess_dir=tmp_path)
    assert report["status"] == "ACTION_REQUIRED"
    assert report["qm_paths_checked"]
    assert not next(check["ok"] for check in report["checks"] if check["name"] == "all_QM_paths_supplied")


def test_existing_report_is_not_overwritten(tmp_path, monkeypatch):
    output = tmp_path / "doctor.json"
    output.write_text("original")
    monkeypatch.setattr(doctor.sys, "argv", ["fmotpiefinal-doctor", "--output", str(output)])
    monkeypatch.setattr(doctor, "inspect_environment", lambda *args: {"status": "PASS"})
    with pytest.raises(FileExistsError):
        doctor.main()
    assert output.read_text() == "original"
