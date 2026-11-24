"""Finalize only existing first-edge endpoint diagnostics after production succeeds."""
from __future__ import annotations

import argparse
import fcntl
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

from hold_gpu_after_run import completion_state
from run_matched_openfe import digest, write_json


def analysis_commands(tools, run, workspace, output, reuse):
    geometry = output / "endpoint_structure.json"
    mm = output / "endpoint_internal_energy.json"
    qm = output / "qm"
    comparison = output / "energy_comparison.json"
    return [
        [str(tools / "report_partial_openfe.py"), "--workspace", str(workspace), "--run", str(run), "--output", str(output / "free_energy_snapshot.json")],
        [str(tools / "analyze_openfe_endpoints.py"), "--run", str(run), "--output", str(geometry)],
        [str(tools / "analyze_endpoint_internal_energy.py"), "--run", str(run), "--diagnostics", str(geometry), "--output", str(mm)],
        [str(tools / "run_endpoint_dftb_internal.py"), "--run", str(run), "--diagnostics", str(geometry),
         "--workspace", str(workspace), "--output", str(qm), "--stride", "40", "--workers", "8", "--reuse", str(reuse)],
        [str(tools / "report_endpoint_energy_comparison.py"), "--mm", str(mm), "--qm", str(qm / "status.json"), "--output", str(comparison)],
        [str(tools / "report_endpoint_repeat_stability.py"), "--comparison", str(comparison), "--output", str(output / "repeat_stability.json")]]


def validate_completed_geometry(path):
    data = json.loads(path.read_text())
    for target in ("BACE", "CDK8"):
        for leg in ("bound", "solution"):
            rows = data["targets"][target][leg]["complete_repeats"]
            counts = {endpoint: sum(row["endpoint"] == endpoint for row in rows) for endpoint in (0, 10)}
            if counts != {0: 3, 10: 3} or any(row["frame_count"] != 41 for row in rows):
                raise ValueError("exactly three completed 41-frame repeats per physical endpoint required")


def run_stage(command, log_path, environment):
    with log_path.open("a") as log:
        child = subprocess.Popen([sys.executable, *command], stdout=log, stderr=subprocess.STDOUT,
                                 env=environment, start_new_session=True)
        try:
            returncode = child.wait()
            if returncode:
                raise RuntimeError(f"analysis stage failed with code {returncode}; inspect {log_path}")
        finally:
            if child.poll() is None:
                os.killpg(child.pid, signal.SIGTERM)
                try:
                    child.wait(timeout=20)
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid, signal.SIGKILL)
                    child.wait()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reuse", type=Path, required=True)
    parser.add_argument("--interval", type=int, default=60)
    args = parser.parse_args()
    if args.interval < 10:
        parser.error("polling interval must be at least ten seconds")
    tools = Path(__file__).resolve().parent
    run, workspace, output, reuse = (path.resolve() for path in (args.run, args.workspace, args.output, args.reuse))
    output.mkdir(parents=True, exist_ok=True)
    commands = analysis_commands(tools, run, workspace, output, reuse)
    scripts = {command[0]: digest(command[0]) for command in commands}
    scripts[__file__] = digest(__file__)
    frozen = {"schema": "first-edge-post-production-audit-v1", "commands": commands,
              "scripts_sha256": scripts, "source_qm_sha256": digest(reuse / "status.json"),
              "production_protocol_sha256": digest(run / "frozen_protocol.json"),
              "maximum_endpoint_singlepoints": 984, "expected_reusable_singlepoints": 738,
              "new_MD_scheduled": False, "new_network_edges_scheduled": False,
              "GPU_required": False, "main_pipeline_changed": False}
    def interrupted(signum, frame):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, interrupted)
    with (output / "watcher.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        frozen_path = output / "frozen_schedule.json"
        if frozen_path.exists() and json.loads(frozen_path.read_text()) != frozen:
            raise ValueError("scheduled audit changed; use new output")
        write_json(frozen_path, frozen)
        def status(state, **details):
            write_json(output / "status.json", {"status": state, "pid": os.getpid(),
                       "updated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), **details})
        try:
            while True:
                state, records = completion_state(run)
                if state == "FAILED":
                    status("FAILED_PRODUCTION", legs=records)
                    return 1
                if state == "READY":
                    break
                status("WAITING_FOR_EXISTING_PRODUCTION", legs=records)
                time.sleep(args.interval)
            for path, expected in scripts.items():
                if digest(path) != expected:
                    raise ValueError(f"scheduled analysis script changed: {path}")
            if digest(reuse / "status.json") != frozen["source_qm_sha256"] or digest(run / "frozen_protocol.json") != frozen["production_protocol_sha256"]:
                raise ValueError("frozen production or reuse source changed")
            environment = os.environ.copy()
            environment.update(CUDA_VISIBLE_DEVICES="", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
            for index, command in enumerate(commands):
                if index == 2:
                    validate_completed_geometry(output / "endpoint_structure.json")
                status("RUNNING_AUDIT", stage=index + 1, command=command)
                run_stage(command, output / f"stage_{index + 1}.stdout", environment)
            status("COMPLETE", scope="First predefined edge per target, three-repeat structure/energy stability diagnostics only; not cohort ranking or a QM free-energy correction.")
            return 0
        except KeyboardInterrupt:
            status("INTERRUPTED")
            return 130
        except Exception as error:
            status("FAILED_AUDIT", error=f"{type(error).__name__}: {error}")
            raise


if __name__ == "__main__":
    raise SystemExit(main())
