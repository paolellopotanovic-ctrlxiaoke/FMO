import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / "tools"))
from finalize_openfe_endpoint_audit import analysis_commands, validate_completed_geometry


def test_finalizer_schedules_analysis_not_new_md(tmp_path):
    commands = analysis_commands(tmp_path, tmp_path, tmp_path, tmp_path, tmp_path)
    assert len(commands) == 6
    assert all("run_matched_openfe.py" not in command[0] for command in commands)
    assert "--reuse" in commands[3] and "--stride" in commands[3]


def test_finalizer_requires_all_independent_repeats(tmp_path):
    rows = [{"endpoint": endpoint, "frame_count": 41} for endpoint in (0, 10) for repeat in range(3)]
    data = {"targets": {target: {leg: {"complete_repeats": list(rows)} for leg in ("solution", "bound")}
                        for target in ("BACE", "CDK8")}}
    path = tmp_path / "geometry.json"
    path.write_text(json.dumps(data))
    validate_completed_geometry(path)
    data["targets"]["CDK8"]["bound"]["complete_repeats"].pop()
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="three completed"):
        validate_completed_geometry(path)
