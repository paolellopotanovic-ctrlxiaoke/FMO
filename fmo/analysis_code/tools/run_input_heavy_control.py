"""Preserve every input ligand heavy atom and relax ligand H before matched tPIE."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import multiprocessing
from pathlib import Path
import time

import numpy as np

from fmotpiefinal.core.bound_geometry import load_bound_ligand
from fmotpiefinal.core.structure import read_pdb, read_sdf
from fmotpiefinal.pipeline.target_stages import receptor_at_minimized_geometry
from run_cdk8_solution_cycle import digest, metrics, write_json
from run_local_sidechain_transfer import directory_key, score_job


def heavy_displacement(original, optimized):
    if len(original.atoms) != len(optimized.atoms) or any(
            left.element != right.element for left, right in zip(original.atoms, optimized.atoms, strict=True)):
        raise ValueError("input and optimized ligand identity mismatch")
    shifts = [float(np.linalg.norm(np.asarray(left.xyz)-right.xyz))
              for left, right in zip(original.atoms, optimized.atoms, strict=True) if left.element.upper() != "H"]
    if not shifts:
        raise ValueError("ligand contains no heavy atoms")
    maximum = max(shifts)
    if maximum > .002:
        raise ValueError("fixed input heavy atoms moved")
    return maximum


def execute(payload):
    cohort, name, output = payload
    parent = Path(output) / cohort["target"] / directory_key(name)
    try:
        status = parent / "preparation_status.json"
        if not status.exists():
            from fmotpiefinal.core.minimization import minimize_complex
            result = minimize_complex(cohort["protein"], cohort["sdf"], name, parent / "minimization",
                restraint_k_kj_mol_a2=10., move_protein=False, fix_ligand_heavy=True,
                restrain_ligand_hydrogens=False)
            write_json(status, {"status": "SUCCEEDED", "complex": str(result.complex_pdb)})
        prepared = json.loads(status.read_text())
        if prepared["status"] != "SUCCEEDED":
            raise ValueError("cached preparation did not succeed")
        original = next(ligand for ligand in read_sdf(cohort["sdf"]) if ligand.name == name)
        optimized = load_bound_ligand(cohort["sdf"], name, prepared["complex"])
        maximum = heavy_displacement(original, optimized)
        protein = read_pdb(cohort["protein"])
        receptor = receptor_at_minimized_geometry(protein, Path(prepared["complex"]))
        if len(protein.atoms) != len(receptor.atoms):
            raise ValueError("fixed receptor atom inventory changed")
        receptor_shift = max(float(np.linalg.norm(np.asarray(left.xyz)-right.xyz))
                             for left, right in zip(protein.atoms, receptor.atoms, strict=True))
        if receptor_shift > .002:
            raise ValueError("fixed receptor moved")
        row = score_job(cohort, name, Path(output), 1)[2]
        row.update(max_input_heavy_displacement_A=maximum, max_fixed_receptor_displacement_A=receptor_shift,
                   optimized_complex_sha256=digest(prepared["complex"]))
    except Exception as error:
        row = {"status": "FAILED", "error": f"{type(error).__name__}: {error}"}
    write_json(parent / "summary.json", row)
    return name, row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--target", default="BACE")
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("positive workers required")
    start = time.time()
    source, output = args.source.resolve(), args.output.resolve()
    cohort = json.loads(source.read_text())["cohorts"][args.target]
    for path, expected in cohort["input_sha256"].items():
        if digest(path) != expected:
            raise ValueError("registered baseline evidence changed")
    cohort = dict(cohort, local_residues=[])
    frozen = {"schema": "input-heavy-fixed-ligand-hydrogen-control-v1", "cohort": cohort,
        "source_sha256": digest(source), "driver_sha256": digest(__file__),
        "minimizer_sha256": digest(Path(__file__).parents[1] / "src/fmotpiefinal/core/minimization.py"),
        "labels_used_for_generation": False, "production_default_changed": False,
        "treatment": "Input ligand heavy coordinates fixed; only ligand H minimized without position restraints; whole receptor/waters fixed; original 5A pocket and numerical controls",
        "interpretation": "Tests k10 heavy-geometry changes, not crystal-pose accuracy or free energy"}
    path = output / "frozen_protocol.json"
    if path.exists() and json.loads(path.read_text()) != frozen:
        raise ValueError("refusing changed frozen experiment")
    write_json(path, frozen)
    rows = {}
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context("spawn")) as pool:
        futures = [pool.submit(execute, (cohort, name, str(output))) for name in cohort["names"]]
        for future in as_completed(futures):
            name, row = future.result()
            rows[name] = row
            write_json(output / "progress.json", {"n_expected": len(futures), "n_finished": len(rows),
                "n_succeeded": sum(item["status"] == "SUCCEEDED" for item in rows.values()), "rows": rows})
            print(name, row["status"], row.get("tpie", row.get("error")), flush=True)
    complete = len(rows) == len(cohort["names"]) and all(row["status"] == "SUCCEEDED" for row in rows.values())
    report = {"status": "COMPLETE_VERIFIED" if complete else "INCOMPLETE", "rows": rows,
              "wall_seconds": time.time()-start, "production_default_changed": False}
    if complete:
        from report_local_sidechain_transfer import load_truth
        truth = load_truth(args.workspace.resolve(), args.target)
        names = sorted(rows)
        if set(truth) != set(names):
            raise ValueError("full cohort evaluation identity mismatch")
        observed = [truth[name] for name in names]
        report["baseline_metrics"] = metrics(observed, [-rows[name]["baseline_tpie"] for name in names])
        report["input_heavy_metrics"] = metrics(observed, [-rows[name]["tpie"] for name in names])
    write_json(output / "evaluation.json", report)
    write_json(output / "status.json", {key: value for key, value in report.items() if key != "rows"})


if __name__ == "__main__":
    main()
