import importlib.util
import json
from pathlib import Path


SPEC = importlib.util.spec_from_file_location("gpu_hold", Path(__file__).parents[1] / "tools/hold_gpu_after_run.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def write_leg(root, target, leg, status):
    path = root / "production" / target / "edge_000" / leg / "status.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"status": status}))


def test_missing_or_running_legs_do_not_reserve(tmp_path):
    assert MODULE.completion_state(tmp_path)[0] == "WAITING"
    write_leg(tmp_path, "CDK8", "solution", "SUCCEEDED")
    write_leg(tmp_path, "CDK8", "bound", "RUNNING")
    assert MODULE.completion_state(tmp_path)[0] == "WAITING"


def test_all_four_legs_required(tmp_path):
    for target in ("CDK8", "BACE"):
        for leg in ("solution", "bound"):
            write_leg(tmp_path, target, leg, "SUCCEEDED")
    assert MODULE.completion_state(tmp_path)[0] == "READY"


def test_failure_does_not_hide_as_completion(tmp_path):
    write_leg(tmp_path, "BACE", "bound", "FAILED")
    assert MODULE.completion_state(tmp_path)[0] == "FAILED"
