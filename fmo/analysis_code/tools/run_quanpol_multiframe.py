"""Run frozen multi-repeat, multi-time endpoint diagnostics without affinity selection."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
from pathlib import Path
import statistics
import subprocess
import sys
import time

import numpy as np

from audit_quanpol_mm_parity import evaluate_export
from build_quanpol_endpoint_pilot import build
from build_quanpol_mm_only import make
from export_quanpol_endpoint_forcefield import export
from quanpol_constants import DEFAULT_SWRB_ANGSTROM, validate_swrb
from report_quanpol_endpoint_pilot import parse
from run_matched_openfe import digest, write_json

def execute(deck, log, args):
    environment = os.environ.copy()
    environment.update({"LD_LIBRARY_PATH": str(Path(sys.executable).resolve().parents[1] / "lib") + ":" + environment.get("LD_LIBRARY_PATH", ""),
                        "OMPI_ALLOW_RUN_AS_ROOT": "1", "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM": "1",
                        "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"})
    started = time.monotonic()
    with log.open("w") as handle:
        process = subprocess.run(["timeout", "180", "tcsh", "-f", str(args.gamess / "rungms"), str(deck), "00", "2"],
                                 cwd=deck.parent, env=environment, stdout=handle, stderr=subprocess.STDOUT)
    result = parse(log)
    if process.returncode or result["status"] != "SUCCEEDED" or result["scf_unconverged"]:
        raise ValueError(f"GAMESS failed: {log}")
    return {**result, "log": str(log), "deck_sha256": digest(deck), "wall_seconds": time.monotonic()-started}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--gamess", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--target", default="CDK8")
    parser.add_argument("--repeat-indices", type=int, nargs="+", default=[0, 1, 2])
    parser.add_argument("--iterations", type=int, nargs="+", default=[400, 800, 1200, 1600, 2000])
    parser.add_argument("--shell", type=float, default=12.0,
                        help="water/ion environment radius around the ligand in angstrom")
    parser.add_argument("--swrb", type=float, default=DEFAULT_SWRB_ANGSTROM,
                        help="QuanPol outer switching radius in angstrom; safe default disables the GAMESS QM-center restraint branch")
    parser.add_argument("--allow-qm-center-restraint", action="store_true",
                        help="allow a SWRB that activates GAMESS' QM-center restraint for explicit diagnostic controls")
    parser.add_argument("--reuse-status", type=Path,
                        help="reuse completed paired frames from a prior compatible status file")
    args = parser.parse_args()
    if not args.repeat_indices or len(set(args.repeat_indices)) != len(args.repeat_indices) or min(args.repeat_indices) < 0:
        parser.error("repeat indices must be unique nonnegative integers")
    if not args.iterations or len(set(args.iterations)) != len(args.iterations) or min(args.iterations) < 0:
        parser.error("frame iterations must be unique nonnegative integers")
    if not np.isfinite(args.shell) or args.shell <= 0:
        parser.error("shell radius must be finite and positive")
    try:
        validate_swrb(args.swrb, args.allow_qm_center_restraint)
    except ValueError as error:
        parser.error(str(error) + "; use the safe default or explicitly opt in for a control")
    args.run, args.gamess, args.output = (path.resolve() for path in (args.run, args.gamess, args.output))
    args.output.mkdir(parents=True, exist_ok=True)
    report = {"schema": "quanpol-multiframe-v2", "status": "RUNNING", "labels_used": False,
              "target": args.target, "times_iterations": args.iterations, "repeats": args.repeat_indices,
              "scope": "one alchemical edge, finite shell DFTB3/MM diagnostic; not affinity ranking or FEP",
              "shell_angstrom": args.shell,
              "swrb_angstrom": args.swrb,
              "qm_center_restraint_allowed": args.allow_qm_center_restraint,
              "driver_sha256": digest(__file__), "results": []}
    if args.reuse_status:
        source_path = args.reuse_status.resolve()
        source = json.loads(source_path.read_text())
        if source.get("status") != "COMPLETE" or source.get("target") != args.target or source.get("labels_used"):
            parser.error("reuse-status must be a complete label-independent run for this target")
        if float(source.get("swrb_angstrom", 1001.0)) != args.swrb:
            parser.error("reuse-status must use the same SWRB setting")
        if float(source.get("shell_angstrom", 12.0)) != args.shell:
            parser.error("reuse-status must use the same water/ion shell setting")
        report["reused_status_sha256"] = hashlib.sha256(source_path.read_bytes()).hexdigest()
        report["results"] = [row for row in source["results"]
                             if int(row["repeat"]) in args.repeat_indices
                             and int(row["iteration"]) in args.iterations]
    write_json(args.output / "frozen_design.json", {key: value for key, value in report.items() if key not in ("results", "status")})
    write_json(args.output / "status.json", report)
    try:
        complete = {(int(row["repeat"]), int(row["iteration"]), row["leg"], int(row["endpoint"]))
                    for row in report["results"]}
        for repeat in args.repeat_indices:
            for iteration in args.iterations:
                directory = args.output / f"repeat_{repeat}" / f"iteration_{iteration:04d}"
                build(argparse.Namespace(run=args.run, gamess=args.gamess, output=directory,
                                         target=args.target, iteration=iteration, repeat=repeat, shell=args.shell))
                tasks = json.loads((directory / "plan.json").read_text())["tasks"]
                for task in tasks:
                    task_key = (repeat, iteration, task["leg"], int(task["endpoint"]))
                    if task_key in complete:
                        print(f"REUSED {len(complete)}/({len(args.repeat_indices)*len(args.iterations)*4}) {task_key}", flush=True)
                        continue
                    export(task)
                    frame_dir = Path(task["deck"]).parent
                    parameter_export = json.loads((frame_dir / "parameter_export.json").read_text())
                    qm_source = Path(parameter_export["deck"])
                    stem = f"r{repeat}_i{iteration}_{task['leg']}_{task['endpoint']}"
                    qm_deck = frame_dir / f"qm_{stem}.inp"
                    mm_deck = frame_dir / f"mm_{stem}.inp"
                    qm_text, replacements = re.subn(
                        r"\bSWRB\s*=\s*[-+0-9.Ee]+", f"SWRB={args.swrb:g}", qm_source.read_text())
                    if replacements != 1:
                        raise ValueError("expected exactly one QuanPol SWRB field")
                    qm_deck.write_text(qm_text)
                    make(qm_source, mm_deck, parameter_export["atom_count"])
                    mm = execute(mm_deck, frame_dir / "mm_matched.stdout", args)
                    reference = evaluate_export(mm_deck)
                    parity = mm["potential_kcal_mol"] - reference["potential_kcal_mol"]
                    if abs(parity) > .1:
                        raise ValueError(f"MM parity failed by {parity} kcal/mol: {frame_dir}")
                    qm = execute(qm_deck, frame_dir / "qm_matched.stdout", args)
                    if "qm_center_kcal_mol" not in qm:
                        raise ValueError(f"GAMESS did not report QM-center energy: {frame_dir}")
                    if not args.allow_qm_center_restraint and abs(qm["qm_center_kcal_mol"]) > 1e-8:
                        raise ValueError(f"unexpected nonzero QM-center restraint energy: {frame_dir}")
                    result = {"repeat": repeat, "iteration": iteration, "leg": task["leg"], "endpoint": task["endpoint"],
                              "trajectory": task["trajectory"], "replica": task["replica"],
                              "qm": qm, "mm": mm, "openmm_ffdata": reference,
                              "mm_parity_difference_kcal_mol": parity,
                              "delta_u_kcal_mol": qm["potential_kcal_mol"]-mm["potential_kcal_mol"]}
                    report["results"].append(result)
                    complete.add(task_key)
                    write_json(frame_dir / "paired_result.json", result)
                    write_json(args.output / "status.json", report)
                    print(f"COMPLETE {len(complete)}/{len(args.repeat_indices)*len(args.iterations)*4} {stem} deltaU={result['delta_u_kcal_mol']:.6f} parity={parity:.6f}", flush=True)
        cycles = []
        for repeat in args.repeat_indices:
            for iteration in args.iterations:
                subset = {(row["leg"], row["endpoint"]): row["delta_u_kcal_mol"] for row in report["results"]
                          if row["repeat"] == repeat and row["iteration"] == iteration}
                cycles.append({"repeat": repeat, "iteration": iteration,
                               "correction_kcal_mol": subset["bound", 10]-subset["bound", 0]-subset["solution", 10]+subset["solution", 0]})
        report["cycles"] = cycles
        report["repeat_means_kcal_mol"] = [float(np.mean([row["correction_kcal_mol"] for row in cycles if row["repeat"] == repeat])) for repeat in args.repeat_indices]
        report["mean_potential_correction_kcal_mol"] = float(np.mean(report["repeat_means_kcal_mol"]))
        report["between_repeat_sd_kcal_mol"] = (
            statistics.stdev(report["repeat_means_kcal_mol"])
            if len(report["repeat_means_kcal_mol"]) > 1 else None
        )
        report["status"] = "COMPLETE"
    except Exception as error:
        report["status"] = "FAILED"
        report["error"] = repr(error)
        raise
    finally:
        write_json(args.output / "status.json", report)


if __name__ == "__main__":
    main()
