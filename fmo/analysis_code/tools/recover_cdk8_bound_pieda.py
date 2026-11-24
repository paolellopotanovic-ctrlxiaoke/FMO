"""Add missing PIEDA output without repeating a converged QM optimization."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import io
import json
import math
from pathlib import Path
import re
import sys
import time

from run_cdk8_bound_qm import converged_coordinates
from run_cdk8_solution_cycle import digest, group, write_json


def pieda_deck(deck):
    if "RUNTYP=ENERGY" not in deck or re.search(r"\bMODMOL\s*=", deck):
        raise ValueError("PIEDA recovery requires an independent full-FMO energy deck")
    if deck.count("IPIEDA=0") != 1:
        raise ValueError("expected exactly one disabled PIEDA control")
    return deck.replace("IPIEDA=0", "IPIEDA=1", 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--wait-seconds", type=int, default=15000)
    args = parser.parse_args()
    workspace, reference, run = args.workspace.resolve(), args.reference.resolve(), args.run.resolve()
    sys.path.insert(0, str(workspace / "src"))
    from sophosqm_baseline.relative_total_energy import execute_deck
    from sophosqm_baseline.parser import parse_pie
    from fmotpiefinal.core.fmo_convergence import final_fmo_convergence
    protocol = json.loads((run / "frozen_protocol.json").read_text())
    names = protocol["names"]
    plan = {"driver_sha256": digest(__file__), "names": names,
        "change": "IPIEDA=0 to IPIEDA=1 only; reuse final ENERGY geometry",
        "changes_hamiltonian": False, "labels_used": False}
    plan_path = run / "pieda_recovery_plan.json"
    if plan_path.exists() and json.loads(plan_path.read_text()) != plan:
        raise ValueError("recovery plan changed")
    write_json(plan_path, plan)

    def recover(name):
        parent = run / name
        summary = parent / "summary.json"
        deadline = time.monotonic() + args.wait_seconds
        while not summary.exists():
            if time.monotonic() >= deadline:
                return name, {"status": "FAILED_WAIT_TIMEOUT"}
            time.sleep(15)
        previous = json.loads(summary.read_text())
        if previous.get("pieda_recovery"):
            return name, previous
        try:
            reference_deck = (reference / f"pcm/{name}/complex/input.inp").read_text()
            ligand_size = len(group((reference / f"pcm/{name}/ligand/input.inp").read_text(), "FMOXYZ")[1].strip().splitlines())
            coordinates = converged_coordinates((parent / "optimization/gamess.log").read_text(), reference_deck, ligand_size)
            complex_result = json.loads((parent / "complex/result.json").read_text())
            ligand_result = json.loads((parent / "ligand/result.json").read_text())
            if any(result["status"] != "SUCCEEDED" for result in (complex_result, ligand_result)):
                raise ValueError("original final energy checks are incomplete")
            source = parent / "complex/input.inp"
            deck = pieda_deck(source.read_text())
            directory = parent / "pieda"
            directory.mkdir(exist_ok=True)
            (directory / "input.inp").write_text(deck)
            result = execute_deck(deck, directory, gamess_dir=complex_result["gamess_dir"],
                library_dir=complex_result["library_dir"], ncores=2, timeout=1800)
            if result["status"] != "SUCCEEDED":
                raise ValueError("PIEDA recovery did not complete")
            audit = final_fmo_convergence(deck, io.StringIO((directory / "gamess.log").read_text()))
            if not math.isclose(result["energy_hartree"], complex_result["energy_hartree"], abs_tol=1e-7, rel_tol=0):
                raise ValueError("PIEDA enabled total energy differs from reference")
            count = int(re.search(r"NFRAG=(\d+)", deck)[1])
            pairs = [pair for pair in parse_pie(directory / "gamess.log")
                     if count in (pair.fragment_i, pair.fragment_j)]
            if len(pairs) != count - 1:
                raise ValueError("incomplete ligand pair inventory")
            if not (parent / "summary_before_pieda_recovery.json").exists():
                write_json(parent / "summary_before_pieda_recovery.json", previous)
            result = {"status": "SUCCEEDED", "energies_hartree": {
                "complex": complex_result["energy_hartree"], "ligand": ligand_result["energy_hartree"]},
                "tpie_kcal_mol": sum(pair.value_kcal_mol for pair in pairs),
                "optimized_coordinates": coordinates, "pieda_recovery": True,
                "strict_convergence": audit, "reference_energy_input_sha256": digest(source),
                "pair_count": len(pairs)}
            write_json(summary, result)
        except (ValueError, IndexError, FileNotFoundError) as error:
            result = {"status": "FAILED_PIEDA_RECOVERY", "error": str(error)}
        print(name, result["status"], result.get("error", ""), flush=True)
        write_json(parent / "pieda_recovery.json", result)
        return name, result

    results = {}
    with ThreadPoolExecutor(max_workers=10) as executor:
        for name, result in executor.map(recover, names):
            results[name] = result
            write_json(run / "pieda_recovery_execution.json", results)


if __name__ == "__main__":
    main()
