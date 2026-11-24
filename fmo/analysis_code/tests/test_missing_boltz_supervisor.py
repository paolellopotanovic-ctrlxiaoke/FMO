"""Background orchestration guards without contacting MSA or launching GPU jobs."""
import importlib.util
import json
from pathlib import Path
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("run_missing_boltz2", ROOT / "tools/run_missing_boltz2.py")
SUPERVISOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SUPERVISOR)


def test_exact_missing_cohort_inventory():
    assert SUPERVISOR.EXPECTED == {"CHK1": 13, "CDK8": 10, "BACE": 36, "HSP90": 11, "PTP1B": 23, "FAAH": 24}
    assert sum(SUPERVISOR.EXPECTED.values()) == 117


def test_frozen_files_cannot_be_replaced_with_different_inputs(tmp_path):
    path = tmp_path / "inputs/frozen.json"
    SUPERVISOR.write_json(path, {"labels_used": False}, frozen=True)
    SUPERVISOR.write_json(path, {"labels_used": False}, frozen=True)
    with pytest.raises(ValueError, match="frozen input changed"):
        SUPERVISOR.write_json(path, {"labels_used": True}, frozen=True)
    assert json.loads(path.read_text()) == {"labels_used": False}


def test_subprocess_failure_is_logged_not_hidden(tmp_path):
    outcome = SUPERVISOR.logged_process([sys.executable, "-c", "print('failed case'); raise SystemExit(3)"],
                                       tmp_path / "failed.log", 10)
    assert outcome["returncode"] == 3
    assert outcome["wall_seconds"] > 0
    assert "failed case" in (tmp_path / "failed.log").read_text()


def test_subprocess_timeout_kills_its_own_process_group(tmp_path):
    outcome = SUPERVISOR.logged_process([sys.executable, "-c", "import time; time.sleep(20)"],
                                       tmp_path / "timeout.log", 0.05)
    assert outcome["returncode"] == 124
    assert outcome["wall_seconds"] < 10
