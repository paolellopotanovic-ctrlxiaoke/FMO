"""Matched input-pose single points diagnose k10 effects without choosing poses."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import io
import json
from pathlib import Path
import re
import sys

import numpy as np
from rdkit import Chem

from run_cdk8_bound_qm import replace_coordinates
from run_cdk8_solution_cycle import digest, group, metrics, split_decks, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    workspace, reference, output = args.workspace.resolve(), args.reference.resolve(), args.output.resolve()
    sys.path.insert(0, str(workspace / "src"))
    from sophosqm_baseline.relative_total_energy import execute_deck
    from sophosqm_baseline.parser import parse_pie
    from fmotpiefinal.core.fmo_convergence import final_fmo_convergence
    from fmotpiefinal.core.solution_cycle import require_consistent_pcm_cycle
    names = json.loads((reference / "frozen_protocol.json").read_text())["names"]
    sdf = workspace / "vendor/public_binding_free_energy_benchmark/fep_benchmark_inputs/structure_inputs/misc/cdk8_koehler_ligands.sdf"
    molecules = {mol.GetProp("_Name").strip(): mol for mol in Chem.SDMolSupplier(str(sdf), removeHs=False)}
    parameters = (reference / f"pcm/{names[0]}/params").resolve(strict=True)
    plan = {"driver_sha256": digest(__file__), "source_sdf_sha256": digest(sdf), "names": names,
        "labels_used_for_generation": False, "arms": ["input_pose_tpie", "input_pose_three_state"],
        "change": "Original SDF ligand coordinates only; same prepared receptor/water/model/charge/H",
        "warning": "Diagnostic bypass of k10, not a recommended default or all-independent-crystal claim"}
    output.mkdir(parents=True, exist_ok=True)
    plan_path = output / "frozen_protocol.json"
    if plan_path.exists() and json.loads(plan_path.read_text()) != plan:
        raise ValueError("input-pose diagnostic plan changed")
    write_json(plan_path, plan)

    def execute(name):
        parent = output / name
        parent.mkdir(exist_ok=True)
        (parent / "params").symlink_to(parameters, target_is_directory=True) if not (parent / "params").exists() else None
        source = (reference / f"pcm/{name}/complex/input.inp").read_text()
        rows = group(source, "FMOXYZ")[1].strip().splitlines()
        molecule = molecules[name]
        size = molecule.GetNumAtoms()
        if [row.split()[0].upper() for row in rows[-size:]] != [atom.GetSymbol().upper() for atom in molecule.GetAtoms()]:
            raise ValueError("input ligand atom identity/order mismatch")
        xyz = molecule.GetConformer().GetPositions()
        positions = [f"{atom.GetSymbol().upper()} {atom.GetAtomicNum()} " + " ".join(f"{value:.10f}" for value in xyz[index])
                     for index, atom in enumerate(molecule.GetAtoms())]
        states = split_decks(replace_coordinates(source, rows[:-size] + positions), parameters)
        states["complex"] = states["complex"].replace("IPIEDA=0", "IPIEDA=1")
        runtime = json.loads((reference / f"pcm/{name}/complex/result.json").read_text())
        results = {}
        try:
            for state in ("complex", "ligand"):
                directory = parent / state
                directory.mkdir(exist_ok=True)
                (directory / "input.inp").write_text(states[state])
                results[state] = execute_deck(states[state], directory, gamess_dir=runtime["gamess_dir"],
                    library_dir=runtime["library_dir"], ncores=2, timeout=1800)
                if results[state]["status"] != "SUCCEEDED":
                    raise ValueError(f"input-pose {state} failed")
                if state == "complex":
                    final_fmo_convergence(states[state], io.StringIO((directory / "gamess.log").read_text()))
            paths = {state: parent / state for state in ("complex", "ligand")}
            paths["receptor"] = reference / "pcm/shared/receptor"
            audit = require_consistent_pcm_cycle(
                {state: (path / "input.inp").read_text() for state, path in paths.items()},
                {state: (path / "gamess.log").read_text() for state, path in paths.items()})
            count = int(re.search(r"NFRAG=(\d+)", source)[1])
            pairs = [pair for pair in parse_pie(parent / "complex/gamess.log")
                     if count in (pair.fragment_i, pair.fragment_j)]
            if len(pairs) != count - 1:
                raise ValueError("input-pose pair inventory incomplete")
            receptor_energy = json.loads((paths["receptor"] / "result.json").read_text())["energy_hartree"]
            result = {"status": "SUCCEEDED", "input_pose_tpie": sum(pair.value_kcal_mol for pair in pairs),
                "input_pose_three_state": (results["complex"]["energy_hartree"] - receptor_energy - results["ligand"]["energy_hartree"]) * 627.5094740631,
                "pcm_audit": audit}
        except ValueError as error:
            result = {"status": "FAILED", "error": str(error)}
        write_json(parent / "summary.json", result)
        print(name, result["status"], result.get("error", ""), flush=True)
        return name, result

    results = {}
    with ThreadPoolExecutor(max_workers=5) as executor:
        for future in as_completed([executor.submit(execute, name) for name in names]):
            name, result = future.result()
            results[name] = result
            write_json(output / "execution.json", results)
    complete = all(result["status"] == "SUCCEEDED" for result in results.values())
    rows = [{"ligand": name, "experimental_dg": float(molecules[name].GetProp("r_exp_dg")), **results[name]} for name in names]
    write_json(output / "evaluation.json", {"status": "COMPLETE" if complete else "INCOMPLETE", "n": len(rows), "rows": rows,
        "metrics": {arm: metrics([row["experimental_dg"] for row in rows], [row[arm] for row in rows]) for arm in plan["arms"]} if complete else {}})


if __name__ == "__main__":
    main()
