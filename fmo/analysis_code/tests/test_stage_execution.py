import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from fmotpiefinal.core.stage_execution import execute_stage


def command():
    return [sys.executable, "-m", "fmotpiefinal.pipeline.target_stages",
            "--stage", "prepare"]


def test_reuses_arguments_and_restores_threads(monkeypatch):
    from fmotpiefinal.pipeline import target_stages as runner

    observed = []

    def stage(argv):
        observed.append(argv)
        os.environ["OPENMM_CPU_THREADS"] = "8"

    monkeypatch.setattr(runner, "main", stage)
    monkeypatch.delenv("FMOTPIE_STAGE_EXECUTION", raising=False)
    monkeypatch.setenv("OPENMM_CPU_THREADS", "1")
    assert execute_stage(command(), Path.cwd()) == "inprocess"
    assert observed == [["--stage", "prepare"]]
    assert os.environ["OPENMM_CPU_THREADS"] == "1"


@pytest.mark.parametrize("failure", [ValueError("bad structure"), SystemExit(2)])
def test_failure_preserves_subprocess_failure_contract(monkeypatch, failure):
    from fmotpiefinal.pipeline import target_stages as runner

    def stage(argv):
        os.environ["OPENMM_CPU_THREADS"] = "8"
        raise failure

    monkeypatch.setattr(runner, "main", stage)
    monkeypatch.delenv("FMOTPIE_STAGE_EXECUTION", raising=False)
    monkeypatch.delenv("OPENMM_CPU_THREADS", raising=False)
    with pytest.raises(subprocess.CalledProcessError) as caught:
        execute_stage(command(), Path.cwd())
    assert caught.value.__cause__ is failure
    assert "OPENMM_CPU_THREADS" not in os.environ


def test_legacy_subprocess_mode(monkeypatch, tmp_path):
    observed = []
    monkeypatch.setenv("FMOTPIE_STAGE_EXECUTION", "subprocess")
    monkeypatch.setattr(subprocess, "run", lambda *args, **kwargs: observed.append((args, kwargs)))
    assert execute_stage(command(), tmp_path) == "subprocess"
    assert observed == [((command(),), {"cwd": tmp_path, "check": True})]


def test_rejects_unregistered_entrypoint(monkeypatch):
    monkeypatch.delenv("FMOTPIE_STAGE_EXECUTION", raising=False)
    with pytest.raises(ValueError, match="registered"):
        execute_stage([sys.executable, "-m", "other"], Path.cwd())


def test_rejects_unknown_execution_mode(monkeypatch):
    monkeypatch.setenv("FMOTPIE_STAGE_EXECUTION", "other")
    with pytest.raises(ValueError, match="must be"):
        execute_stage(command(), Path.cwd())


def test_holdout_import_does_not_load_minimizer():
    result = subprocess.run([sys.executable, "-c",
        "import sys,json; import fmotpiefinal.pipeline.target_stages; "
        "print(json.dumps('fmotpiefinal.core.minimization' in sys.modules))"],
        capture_output=True, text=True, check=True)
    assert json.loads(result.stdout) is False
