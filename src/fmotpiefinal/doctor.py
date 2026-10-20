"""Offline preflight for preparation dependencies and the explicit GAMESS runtime."""
from __future__ import annotations

import argparse
import importlib
import json
import os
from pathlib import Path
import shutil
import sys


MODULES = ("numpy", "scipy", "rdkit", "gemmi", "openmm", "pdbfixer",
           "openmmforcefields", "openff.toolkit", "openff.nagl")


def inspect_environment(gamess_dir=None, library_dir=None, parameter_dir=None) -> dict:
    checks = []
    for name in MODULES:
        try:
            importlib.import_module(name)
            checks.append({"name": name, "ok": True})
        except Exception as error:
            checks.append({"name": name, "ok": False,
                           "detail": f"{type(error).__name__}: {error}"})
    executable = shutil.which("pdb2pqr")
    checks.append({"name": "pdb2pqr", "ok": executable is not None,
                   "detail": executable or "Install PDB2PQR and add pdb2pqr to PATH"})
    supplied = (gamess_dir, library_dir, parameter_dir)
    qm_requested = any(value is not None for value in supplied)
    if qm_requested:
        checks.append({"name": "all_QM_paths_supplied", "ok": all(value is not None for value in supplied)})
        for name, directory, child, executable_required in (
                ("GAMESS launcher", gamess_dir, "rungms", False),
                ("GAMESS binary", gamess_dir, "gamess.00.x", True)):
            path = Path(directory) / child if directory is not None else None
            valid = path is not None and path.is_file()
            checks.append({"name": name, "ok": bool(valid and
                           (not executable_required or os.access(path, os.X_OK))),
                           "detail": str(path) if path is not None else "missing --gamess-dir"})
        checks.append({"name": "runtime library directory", "ok": library_dir is not None and Path(library_dir).is_dir()})
        parameters = Path(parameter_dir) if parameter_dir is not None else None
        checks.append({"name": "DFTB parameter directory", "ok": parameters is not None and
                       parameters.is_dir() and (parameters / "C-C.skf").is_file()})
        checks.append({"name": "GAMESS tcsh shell", "ok": Path("/bin/tcsh").is_file() and
                       os.access("/bin/tcsh", os.X_OK)})
    return {"schema": "fmotpiefinal-environment-preflight-v1", "python": sys.version.split()[0],
            "status": "PASS" if all(check["ok"] for check in checks) else "ACTION_REQUIRED",
            "qm_paths_checked": qm_requested, "checks": checks,
            "scope": "Offline dependency/path inspection only; does not prove parameter coverage, model asset availability or successful QM",
            "next_step": "Run a small complete blind-run target before batch calculation"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gamess-dir", type=Path)
    parser.add_argument("--library-dir", type=Path)
    parser.add_argument("--dftb-parameter-dir", type=Path)
    parser.add_argument("--output", type=Path, help="optional JSON report; existing files are not overwritten")
    args = parser.parse_args()
    report = inspect_environment(args.gamess_dir, args.library_dir, args.dftb_parameter_dir)
    content = json.dumps(report, indent=2) + "\n"
    if args.output:
        with args.output.open("x") as handle:
            handle.write(content)
    print(content, end="")
    raise SystemExit(0 if report["status"] == "PASS" else 2)


if __name__ == "__main__":
    main()
