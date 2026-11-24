"""Resume a frozen transfer once its original driver exits, then write audit report."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from run_cdk8_solution_cycle import digest, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wait-pid", required=True, type=int)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--workspace", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()
    root = args.run.resolve()
    status = {"status": "WAITING_FOR_ORIGINAL_DRIVER", "wait_pid": args.wait_pid,
              "pid": os.getpid(), "openbind_requested": False}
    write_json(root / "supervisor_status.json", status)
    while True:
        try:
            os.kill(args.wait_pid, 0)
        except ProcessLookupError:
            break
        time.sleep(60)
    plan = json.loads((root / "frozen_protocol.json").read_text())
    tools = Path(__file__).resolve().parent
    status.update(status="REVALIDATING_WITH_TERMINAL_CAP_ALIAS_SUPPORT",
        cap_identity_implementation_sha256=digest(tools.parent / "src/fmotpiefinal/pipeline/target_stages.py"))
    write_json(root / "supervisor_status.json", status)
    command = [sys.executable, str(tools / "run_local_sidechain_transfer.py"),
        "--workspace", str(args.workspace.resolve()), "--baseline", plan["baseline"],
        "--output", str(root), "--targets", *plan["targets"],
        "--prepare-workers", "4", "--qm-workers", "8", "--ncores", "2"]
    with (root / "cap_alias_resume.stdout").open("w") as handle:
        returncode = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT).returncode
    status.update(status="REPORTING" if returncode == 0 else "FAILED_RESUME", resume_returncode=returncode)
    write_json(root / "supervisor_status.json", status)
    if returncode != 0:
        return
    command = [sys.executable, str(tools / "report_local_sidechain_transfer.py"),
        "--workspace", str(args.workspace.resolve()), "--run", str(root), "--output", str(args.report.resolve())]
    with (root / "report.stdout").open("w") as handle:
        returncode = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT).returncode
    status.update(status="FINISHED" if returncode == 0 else "FAILED_REPORT", report_returncode=returncode)
    if returncode == 0:
        status["calculation_status"] = json.loads((root / "final_status.json").read_text())["status"]
        status["report_sha256"] = digest(args.report.resolve())
    write_json(root / "supervisor_status.json", status)


if __name__ == "__main__":
    main()
