"""Recover failed matched-cycle jobs without changing geometry or energy model."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import re
import sys

from run_cdk8_solution_cycle import write_json, digest


def recovery_deck(deck, optimization, profile="density_soscf"):
    if optimization:
        if profile == "charge_diis_rfo":
            return deck.replace("NSTEP=200 OPTTOL=0.000100",
                "NSTEP=500 OPTTOL=0.000100 METHOD=RFO TRMIN=0.001", 1)
        return deck.replace("NSTEP=200 ", "NSTEP=1000 ", 1)
    if profile == "charge_diis_rfo":
        return deck.replace("DAMPEX=4.00 $END", "DAMPEX=4.00 ITYPMX=2 $END")
    deck = deck.replace("DAMPEX=4.00 $END", "DAMPEX=4.00 ITYPMX=-1 $END")
    if profile == "diis_level_shift":
        deck = deck.replace("MCONV(2)=530 MCONV(4)=530", "MCONV(2)=532 MCONV(4)=532")
        return deck.replace("DIIS=.FALSE. SOSCF=.TRUE. $END",
                            "DIIS=.TRUE. SOSCF=.FALSE.\n MAXDII=20 ETHRSH=10.0 FSHIFT=0.5 $END")
    if profile == "density_diis":
        deck = deck.replace("MCONV(2)=530 MCONV(4)=530", "MCONV(2)=528 MCONV(4)=528")
        return deck.replace("DIIS=.FALSE. SOSCF=.TRUE. $END",
                            "DIIS=.TRUE. SOSCF=.FALSE.\n MAXDII=20 ETHRSH=10.0 $END")
    deck = deck.replace("MCONV(2)=530 MCONV(4)=530", "MCONV(2)=576 MCONV(4)=576")
    if profile == "soscf_early":
        deck = deck.replace("CONV=1.0E-8 DIIS=.FALSE. SOSCF=.TRUE. $END",
                            "CONV=1.0E-8\n DIIS=.FALSE. SOSCF=.TRUE. SOGTOL=1.0 $END")
    return deck


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--recovery-name", default="recovery_v2")
    parser.add_argument("--gas-only", action="store_true")
    parser.add_argument("--profile", choices=("density_soscf", "charge_diis_rfo", "soscf_early", "density_diis", "diis_level_shift"), default="density_soscf")
    args = parser.parse_args()
    workspace, root = args.workspace.resolve(), args.run.resolve()
    sys.path.insert(0, str(workspace / "src"))
    from sophosqm_baseline.relative_total_energy import execute_deck, parse_converged_ligand_optimization
    frozen = json.loads((root / "frozen_protocol.json").read_text())
    tasks = {}
    for key in frozen["deck_sha256"]:
        if args.gas_only and not key.startswith("gas/"):
            continue
        if any(json.loads(path.read_text())["status"] == "SUCCEEDED"
               for path in root.glob("recovery*/" + key + "/result.json")):
            continue
        result = json.loads((root / key / "result.json").read_text())
        if result["status"] != "SUCCEEDED":
            deck = (root / key / "input.inp").read_text()
            if digest(root / key / "input.inp") != frozen["deck_sha256"][key]:
                raise ValueError("original frozen deck changed")
            tasks[key] = (recovery_deck(deck, key.startswith("optimized/"), args.profile), result)
    plan = {"schema": "cdk8-numerical-recovery-v1", "jobs": list(tasks),
        "driver_sha256": digest(__file__), "labels_used": False,
        "profile": args.profile,
        "policy": "All failed gas full-energy jobs: density SOSCF, strict unchanged tolerances, MAXIT200 (GAMESS limit). All failed ligand optimizations: NSTEP1000 from original bound pose; unchanged tolerance and model. Preserve original attempts."}
    recovery = root / args.recovery_name
    write_json(recovery / "plan.json", plan)

    def execute(key, item):
        deck, previous = item
        directory = recovery / key
        directory.mkdir(parents=True, exist_ok=True)
        parameter_source = (root / key).parent / "params"
        link = directory.parent / "params"
        if not link.exists():
            link.symlink_to(parameter_source.resolve(), target_is_directory=True)
        (directory / "input.inp").write_text(deck)
        result = execute_deck(deck, directory, gamess_dir=previous["gamess_dir"],
            library_dir=previous["library_dir"], timeout=1800, ncores=2)
        if key.startswith("optimized/"):
            try:
                bound = (root / key.replace("optimized/", "pcm/") / "input.inp").read_text()
                coordinates, energy = parse_converged_ligand_optimization(
                    (directory / "gamess.log").read_text(errors="replace"), ligand_deck=bound)
                if result.get("returncode") != 0:
                    raise ValueError("nonzero process exit")
                result.update(status="SUCCEEDED", energy_hartree=energy, optimized_coordinates=coordinates)
                result.pop("error", None)
            except Exception as error:
                result.update(status="FAILED_OPTIMIZATION", error=repr(error))
        result["original_job"] = key
        write_json(directory / "result.json", result)
        print(key, result["status"], round(result["wall_seconds"], 2), flush=True)
        return key, result

    recovered = {}
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(execute, key, item) for key, item in tasks.items()]
        for future in as_completed(futures):
            key, result = future.result()
            recovered[key] = result
            write_json(recovery / "progress.json", {"requested": len(tasks),
                "completed": len(recovered), "succeeded": sum(value["status"] == "SUCCEEDED" for value in recovered.values())})
    write_json(recovery / "execution.json", recovered)


if __name__ == "__main__":
    main()
