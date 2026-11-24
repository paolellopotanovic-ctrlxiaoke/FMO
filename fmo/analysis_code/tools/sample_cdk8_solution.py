"""Fixed-seed solution conformer sensitivity, without selecting by affinity."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import sys

from rdkit import Chem
from rdkit.Chem import AllChem

from run_cdk8_solution_cycle import group, write_json, digest, metrics, HARTREE_TO_KCAL


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--conformers", type=int, default=4)
    parser.add_argument("--recover", action="store_true")
    parser.add_argument("--recovery-profile", choices=("qa", "rfo_perturbed"), default="qa")
    args = parser.parse_args()
    workspace, root = args.workspace.resolve(), args.run.resolve()
    sys.path.insert(0, str(workspace / "src"))
    from sophosqm_baseline.relative_total_energy import execute_deck, parse_converged_ligand_optimization
    frozen = json.loads((root / "frozen_protocol.json").read_text())
    sdf = workspace / "vendor/public_binding_free_energy_benchmark/fep_benchmark_inputs/structure_inputs/misc/cdk8_koehler_ligands.sdf"
    molecules = {molecule.GetProp("_Name").strip(): molecule for molecule in Chem.SDMolSupplier(str(sdf), removeHs=False)}
    sampling = root / "solution_sampling4"
    if args.recover:
        recover(workspace, root, sampling, execute_deck, parse_converged_ligand_optimization, args.recovery_profile)
        return
    tasks = {}
    inputs = {}
    for name in frozen["names"]:
        molecule = Chem.Mol(molecules[name])
        molecule.RemoveAllConformers()
        parameters = AllChem.ETKDGv3()
        parameters.randomSeed = 20261021
        parameters.pruneRmsThresh = -1.0
        parameters.numThreads = 1
        identifiers = list(AllChem.EmbedMultipleConfs(molecule, numConfs=args.conformers, params=parameters))
        if len(identifiers) != args.conformers:
            raise ValueError("complete deterministic embedding required")
        properties = AllChem.MMFFGetMoleculeProperties(molecule)
        if properties is None:
            raise ValueError("MMFF parameters missing")
        bound = (root / f"pcm/{name}/ligand/input.inp").read_text()
        optimized = (root / f"optimized/{name}/ligand/input.inp").read_text().replace("NSTEP=200", "NSTEP=500")
        xyz_group = group(optimized, "FMOXYZ")
        atom_lines = group(bound, "FMOXYZ")[1].strip().splitlines()
        for identifier in identifiers:
            forcefield = AllChem.MMFFGetMoleculeForceField(molecule, properties, confId=identifier)
            forcefield.Minimize(maxIts=500)
            coordinates = molecule.GetConformer(identifier).GetPositions()
            xyz = [" ".join(line.split()[:2]) + " " + " ".join(f"{value:.10f}" for value in position)
                   for line, position in zip(atom_lines, coordinates, strict=True)]
            deck = optimized[:xyz_group.start()] + "$FMOXYZ\n" + "\n".join(xyz) + "\n $END" + optimized[xyz_group.end():]
            key = f"{name}/c{identifier}"
            directory = sampling / key
            directory.mkdir(parents=True, exist_ok=True)
            link = directory.parent / "params"
            if not link.exists():
                link.symlink_to((root / f"pcm/{name}/params").resolve(), target_is_directory=True)
            (directory / "input.inp").write_text(deck)
            inputs[key] = digest(directory / "input.inp")
            tasks[key] = (deck, bound)
    write_json(sampling / "frozen_protocol.json", {"schema": "cdk8-solution-multistart-v1",
        "driver_sha256": digest(__file__), "n_ligands": 10, "n_starts": args.conformers,
        "seed": 20261021, "input_sha256": inputs, "labels_used_for_generation": False,
        "policy": "Four ETKDGv3 starts per ligand, MMFF preconditioning, identical QM/PCM IDISP0 RFO minima; take minimum energy over all starts including bound-start reference; no label selection; not conformational free energy"})
    runtime = json.loads((root / "pcm/12/ligand/result.json").read_text())

    def execute(key, values):
        deck, bound = values
        directory = sampling / key
        result = execute_deck(deck, directory, gamess_dir=runtime["gamess_dir"],
            library_dir=runtime["library_dir"], timeout=1200, ncores=2)
        try:
            coordinates, energy = parse_converged_ligand_optimization(
                (directory / "gamess.log").read_text(errors="replace"), ligand_deck=bound)
            if result.get("returncode") != 0:
                raise ValueError("nonzero process exit")
            result.update(status="SUCCEEDED", energy_hartree=energy, optimized_coordinates=coordinates)
            result.pop("error", None)
        except Exception as error:
            result.update(status="FAILED_OPTIMIZATION", error=repr(error))
        write_json(directory / "result.json", result)
        print(key, result["status"], round(result["wall_seconds"], 2), flush=True)
        return key, result

    results = {}
    with ThreadPoolExecutor(max_workers=8) as executor:
        for future in as_completed([executor.submit(execute, key, values) for key, values in tasks.items()]):
            key, result = future.result()
            results[key] = result
            write_json(sampling / "progress.json", {"completed": len(results), "requested": len(tasks),
                "succeeded": sum(value["status"] == "SUCCEEDED" for value in results.values())})
    incomplete = [key for key, value in results.items() if value["status"] != "SUCCEEDED"]
    rows = []
    baseline = json.loads((root / "evaluation.json").read_text())
    for row in baseline["rows"]:
        name = row["ligand"]
        candidates = [json.loads((root / f"optimized/{name}/ligand/result.json").read_text())["energy_hartree"]]
        candidates.extend(value["energy_hartree"] for key, value in results.items()
                          if key.split("/")[0] == name and value["status"] == "SUCCEEDED")
        bound = json.loads((root / f"pcm/{name}/ligand/result.json").read_text())["energy_hartree"]
        strain = (bound - min(candidates)) * HARTREE_TO_KCAL
        rows.append({"ligand": name, "experimental_dg": float(molecules[name].GetProp("r_exp_dg")),
            "solution_min_strain": strain, "three_state_pcm_multistart": row["three_state_pcm_bound"] + strain,
            "tpie_plus_multistart_strain": row["tpie_pcm"] + strain,
            "n_successful_minima_including_bound_start": len(candidates)})
    report = {"n": 10, "n_jobs": len(tasks), "failed_starts": incomplete, "rows": rows,
              "metrics": {}, "warning": "Multistart minimum-energy proxy, not entropy or free energy; not blind validation"}
    if not incomplete:
        report["metrics"] = {arm: metrics([row["experimental_dg"] for row in rows], [row[arm] for row in rows])
            for arm in ("three_state_pcm_multistart", "tpie_plus_multistart_strain")}
    write_json(sampling / "evaluation.json", report)
    print(json.dumps({"failed": incomplete, "metrics": report["metrics"]}, indent=2), flush=True)


def recover(workspace, root, sampling, execute_deck, parse_optimization, profile):
    import math
    results = {}
    pending = []
    for path in sorted(sampling.glob("*/c*/result.json")):
        result = json.loads(path.read_text())
        key = str(path.parent.relative_to(sampling))
        for recovered_path in sorted(path.parent.glob("recovery*/result.json")):
            recovered_result = json.loads(recovered_path.read_text())
            if recovered_result["status"] == "SUCCEEDED":
                result = recovered_result
        if result["status"] != "SUCCEEDED":
            pending.append((key, result))
        else:
            results[key] = result
    write_json(sampling / f"recovery_{profile}_plan.json", {"failed_keys": [key for key, _ in pending],
        "labels_used": False, "profile": profile,
        "policy": "Same original ETKDG/MMFF geometry; QA trust radius0.1 or deterministic 0.01A perturbation with RFO; unchanged energy model and gradient tolerance, NSTEP500"})
    for key, previous in pending:
        directory = sampling / key / ("recovery" if profile == "qa" else "recovery_perturbed")
        directory.mkdir(parents=True, exist_ok=True)
        name = key.split("/")[0]
        link = directory.parent / "params"
        if not link.exists():
            link.symlink_to((root / f"pcm/{name}/params").resolve(), target_is_directory=True)
        deck = (sampling / key / "input.inp").read_text()
        if profile == "qa":
            deck = deck.replace("METHOD=RFO TRMIN=0.001", "METHOD=QA TRMIN=0.0001 DXMAX=0.1")
        else:
            xyz_group = group(deck, "FMOXYZ")
            coordinates = []
            for index, line in enumerate(xyz_group[1].strip().splitlines()):
                fields = line.split()
                position = [float(value) + 0.01 * math.sin((index + 1) * (axis + 1) * 1.7)
                            for axis, value in enumerate(fields[2:])]
                coordinates.append(" ".join(fields[:2]) + " " + " ".join(f"{value:.10f}" for value in position))
            deck = deck[:xyz_group.start()] + "$FMOXYZ\n" + "\n".join(coordinates) + "\n $END" + deck[xyz_group.end():]
        (directory / "input.inp").write_text(deck)
        result = execute_deck(deck, directory, gamess_dir=previous["gamess_dir"],
            library_dir=previous["library_dir"], timeout=1200, ncores=2)
        try:
            bound = (root / f"pcm/{name}/ligand/input.inp").read_text()
            coordinates, energy = parse_optimization((directory / "gamess.log").read_text(errors="replace"), ligand_deck=bound)
            if result.get("returncode") != 0:
                raise ValueError("nonzero process exit")
            result.update(status="SUCCEEDED", energy_hartree=energy, optimized_coordinates=coordinates)
            result.pop("error", None)
        except Exception as error:
            result.update(status="FAILED_OPTIMIZATION", error=repr(error))
        write_json(directory / "result.json", result)
        results[key] = result
        print(key, result["status"], flush=True)
    original = json.loads((sampling / "evaluation.json").read_text())
    for row in original["rows"]:
        name = row["ligand"]
        candidates = [json.loads((root / f"optimized/{name}/ligand/result.json").read_text())["energy_hartree"]]
        candidates.extend(value["energy_hartree"] for key, value in results.items()
            if key.split("/")[0] == name and value["status"] == "SUCCEEDED")
        bound = json.loads((root / f"pcm/{name}/ligand/result.json").read_text())["energy_hartree"]
        strain = (bound - min(candidates)) * HARTREE_TO_KCAL
        previous_strain = row["solution_min_strain"]
        for arm in ("three_state_pcm_multistart", "tpie_plus_multistart_strain"):
            row[arm] += strain - previous_strain
        row.update(solution_min_strain=strain, n_successful_minima_including_bound_start=len(candidates))
    original["failed_starts"] = [key for key, result in results.items() if result["status"] != "SUCCEEDED"]
    if not original["failed_starts"]:
        original["metrics"] = {arm: metrics([row["experimental_dg"] for row in original["rows"]],
            [row[arm] for row in original["rows"]]) for arm in ("three_state_pcm_multistart", "tpie_plus_multistart_strain")}
    write_json(sampling / "final_evaluation.json", original)
    print(json.dumps({"failed": original["failed_starts"], "metrics": original["metrics"]}, indent=2))


if __name__ == "__main__":
    main()
