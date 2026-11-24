import json
import sys

import pytest

from fmotpiefinal.pipeline import cli


def test_cli_reports_known_input_failures_without_traceback(tmp_path, monkeypatch, capsys):
    output = tmp_path / "run"
    output.mkdir()
    monkeypatch.setattr(sys, "argv", ["fmotpiefinal", "blind-run", "--output", str(output)])

    def invalid_input():
        raise ValueError("unsupported component: preserve source and resolve model")

    monkeypatch.setattr(cli, "_main", invalid_input)
    with pytest.raises(SystemExit) as error:
        cli.main()
    assert error.value.code == 2
    report = json.loads(capsys.readouterr().err)
    assert report["status"] == "REVIEW_REQUIRED"
    assert not report["ranking_eligible"]
    assert "unsupported component" in report["message"]
    assert json.loads((output / "tool_error.json").read_text())["error_type"] == "ValueError"


def test_cli_does_not_overwrite_existing_failure_evidence(tmp_path, monkeypatch):
    output = tmp_path / "run"
    output.mkdir()
    prior = output / "tool_error.json"
    prior.write_text("pinned previous diagnostic\n")
    monkeypatch.setattr(sys, "argv", ["fmotpiefinal", "--output", str(output)])
    monkeypatch.setattr(cli, "_main", lambda: (_ for _ in ()).throw(FileNotFoundError("missing input")))
    with pytest.raises(SystemExit):
        cli.main()
    assert prior.read_text() == "pinned previous diagnostic\n"
