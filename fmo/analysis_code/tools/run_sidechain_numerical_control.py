"""Whole-cohort strict numerical control without geometry or model changes."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import io
import json
import multiprocessing
from pathlib import Path
import re
import time

from fmotpiefinal.core.fmo_convergence import final_fmo_convergence
from fmotpiefinal.core.parser import parse_pie
from fmotpiefinal.core.runner import run_gamess
from run_cdk8_solution_cycle import digest, write_json
from run_local_sidechain_transfer import directory_key


def strict_damping_control(deck):
    pattern = r"\bCNVDMP\s*=\s*\d+\b"
    if len(re.findall(pattern, deck, re.I)) != 1:
        raise ValueError("exactly one registered CNVDMP control required")
    return re.sub(pattern, "CNVDMP=0", deck, flags=re.I)


def strict_diis_control(deck):
    candidate = strict_damping_control(deck)
    for pattern, replacement in ((r"\bDIIS=\.FALSE\.", "DIIS=.TRUE."),
                                 (r"\bSOSCF=\.TRUE\.", "SOSCF=.FALSE."),
                                 (r"\bMCONV\(2\)=576\b", "MCONV(2)=530"),
                                 (r"\bMCONV\(4\)=576\b", "MCONV(4)=530")):
        candidate, count = re.subn(pattern, replacement, candidate, flags=re.I)
        if count != 1:
            raise ValueError("one registered SOSCF numerical setting required")
    return candidate


def tight_soscf_control(deck):
    candidate = strict_damping_control(deck)
    section = re.search(r"\$SCF\s+(.*?)\$END", candidate, re.S | re.I)
    if section is None:
        raise ValueError("one SCF section required")
    replacement, count = re.subn(r"\bCONV=1\.0E-8\b", "CONV=1.0E-10", section[0], flags=re.I)
    if count != 1:
        raise ValueError("registered SCF 1e-8 tolerance required")
    return candidate[:section.start()]+replacement+candidate[section.end():]


def execute(payload):
    name, cohort, source_root, output, control = payload
    parent = Path(output) / directory_key(name)
    source = Path(source_root) / directory_key(name) / "input.inp"
    try:
        render = {"strict_soscf": strict_damping_control, "strict_diis": strict_diis_control,
            "tight_soscf": tight_soscf_control}[control]
        deck = render(source.read_text())
        parent.mkdir(parents=True, exist_ok=True)
        path = parent / "input.inp"
        if path.exists() and path.read_text() != deck:
            raise ValueError("refusing changed numerical-control input")
        path.write_text(deck)
        count = int(re.search(r"NFRAG=(\d+)", deck)[1])
        runtime = cohort["runtime"][name]
        record_path = parent / "run/result.json"
        if record_path.exists():
            record = json.loads(record_path.read_text())
            if record["input_sha256"] != digest(path):
                raise ValueError("numerical-control cached input mismatch")
        else:
            record = run_gamess(path, parent / "run", gamess_dir=runtime["gamess_dir"],
                library_dir=runtime["library_dir"], ligand_fragment=count, ncores=2,
                timeout_seconds=1800, parameter_dir=cohort["parameters"])
        if record["status"] != "SUCCEEDED":
            raise ValueError(f"QM status: {record['status']}")
        log = parent / "run/gamess.log"
        convergence = final_fmo_convergence(deck, io.StringIO(log.read_text()))
        pairs = [pair for pair in parse_pie(log) if count in (pair.fragment_i, pair.fragment_j)]
        if len(pairs) != count-1:
            raise ValueError("incomplete ligand pair inventory")
        result = {"status": "SUCCEEDED", "tpie": sum(pair.value_kcal_mol for pair in pairs),
            "pair_count": len(pairs), "convergence": convergence, "input_sha256": digest(path),
            "source_input_sha256": digest(source), "log_sha256": digest(log),
            "qm_wall_seconds": record.get("wall_seconds")}
    except Exception as error:
        result = {"status": "FAILED", "error": f"{type(error).__name__}: {error}"}
    write_json(parent / "summary.json", result)
    return name, result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--target", default="FXa")
    parser.add_argument("--workers", default=8, type=int)
    parser.add_argument("--register-only", action="store_true")
    parser.add_argument("--control", choices=("strict_soscf", "strict_diis", "tight_soscf"), default="strict_soscf")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("workers must be positive")
    source = args.source.resolve()
    output = args.output.resolve()
    source_plan = json.loads((source / "frozen_protocol.json").read_text())
    cohort = source_plan["cohorts"][args.target]
    source_root = source / args.target
    inputs = {str(source_root / directory_key(name) / "input.inp"):
        digest(source_root / directory_key(name) / "input.inp") for name in cohort["names"]}
    plan = {"schema": "sidechain-numerical-control-v1", "target": args.target,
        "cohort": cohort, "source_root": str(source_root), "input_sha256": inputs,
        "source_progress_sha256": digest(source / "progress.json"),
        "source_protocol_sha256": digest(source / "frozen_protocol.json"),
        "code_sha256": {str(Path(__file__).resolve()): digest(__file__)},
        "control": args.control,
        "settings": "CNVDMP=0; strict_diis additionally uses DIIS plus damping (MCONV=530) instead of SOSCF (576)",
        "tight_soscf_settings": "SCF CONV tightened from 1e-8 to 1e-10; outer FMO CONV remains 1e-7",
        "invariants": "geometry, atoms, fragment charges, QM/PCM model unchanged; no tolerance relaxed",
        "same_geometry_tpie_agreement_tolerance_kcal_mol": 0.1,
        "full_cohort_required": True, "labels_used_for_selection": False,
        "production_default_changed": False}
    frozen = output / "frozen_protocol.json"
    if frozen.exists() and json.loads(frozen.read_text()) != plan:
        raise ValueError("refusing changed frozen numerical protocol")
    write_json(frozen, plan)
    if args.register_only:
        return
    started = time.monotonic()
    rows = {}
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context("spawn")) as pool:
        futures = [pool.submit(execute, (name, cohort, str(source_root), str(output), args.control)) for name in cohort["names"]]
        for future in as_completed(futures):
            name, result = future.result()
            rows[name] = result
            write_json(output / "progress.json", rows)
            print(name, result["status"], result.get("tpie", result.get("error")), flush=True)
    previous = json.loads((source / "progress.json").read_text())[args.target]
    compared = {name: rows[name]["tpie"]-previous[name]["tpie"] for name in cohort["names"]
        if rows[name]["status"] == "SUCCEEDED" and previous[name]["status"] == "SUCCEEDED"}
    failures = {name: row for name, row in rows.items() if row["status"] != "SUCCEEDED"}
    result = {"status": "COMPLETE" if not failures else "INCOMPLETE_WITH_FAILURES",
        "n_requested": len(cohort["names"]), "n_succeeded": len(rows)-len(failures),
        "failures": failures, "same_geometry_tpie_differences_kcal_mol": compared,
        "same_geometry_agreement_passed": bool(compared) and max(map(abs, compared.values())) <= 0.1,
        "wall_seconds": time.monotonic()-started, "production_default_changed": False}
    write_json(output / "final_status.json", result)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
