"""Reuse a target process without changing stage arguments or failure gates."""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys


def execute_stage(command, cwd):
    mode = os.environ.get("FMOTPIE_STAGE_EXECUTION", "inprocess")
    if mode not in {"inprocess", "subprocess"}:
        raise ValueError("FMOTPIE_STAGE_EXECUTION must be inprocess or subprocess")
    prefix = [sys.executable, "-m", "fmotpiefinal.pipeline.target_stages"]
    if list(command[:3]) != prefix:
        raise ValueError("stage execution requires the registered holdout entry point")
    if mode == "subprocess":
        subprocess.run(command, cwd=cwd, check=True)
        return mode
    if Path.cwd().resolve() != Path(cwd).resolve():
        raise ValueError("in-process stages require the workflow working directory")
    previous_threads = os.environ.get("OPENMM_CPU_THREADS")
    try:
        from fmotpiefinal.pipeline.target_stages import main

        main(list(command[3:]))
    except SystemExit as error:
        if error.code not in (None, 0):
            raise subprocess.CalledProcessError(error.code, command) from error
    except Exception as error:
        raise subprocess.CalledProcessError(1, command) from error
    finally:
        if previous_threads is None:
            os.environ.pop("OPENMM_CPU_THREADS", None)
        else:
            os.environ["OPENMM_CPU_THREADS"] = previous_threads
    return mode
