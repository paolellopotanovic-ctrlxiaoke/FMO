"""Frozen, unit-coefficient local ligand-strain transfer on all saved cohorts."""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import io
import json
import multiprocessing
from pathlib import Path
import sys
import time

import numpy as np
from rdkit import Chem

from audit_cycle_ligand_chemistry import geometry_screen, stereochemistry
from fmotpiefinal.core.fmo_convergence import final_fmo_convergence
from fmotpiefinal.core.parser import parse_pie
from fmotpiefinal.core.runner import charge_aware_pair_outliers
from recover_cdk8_solution_cycle import recovery_deck
from run_cdk8_solution_cycle import digest, group, write_json
from run_cohort_solution_cycle import cycle_decks, execute
from run_local_sidechain_transfer import directory_key


REFERENCE_RUNS = {
    "BACE": "BACE_matched_solution_cycle_20261024/ligand_minimum_rfo_recovery",
    "Galectin3": "Galectin3_matched_solution_cycle_20261024",
    "CDK8": "CDK8_current_main_solution_cycle_20261024",
}
RECOVERY_POLICY = {
    "optimization": ["original_200", "rfo_500", "qa_500"],
    "single_point": ["original", "density_soscf", "density_diis"],
    "selection": "First strict converged, chemically valid result; never minimum energy or label agreement",
    "thresholds": "Unchanged OPTTOL=0.000100, SCF CONV and FMO CONV; all restarts from original coordinates",
}


def optimization_recovery(deck, profile):
    if profile == "original_200":
        return deck
    if profile == "rfo_500":
        return recovery_deck(deck, True, "charge_diis_rfo")
    if profile == "qa_500":
        return deck.replace("NSTEP=200 OPTTOL=0.000100", "NSTEP=500 OPTTOL=0.000100 METHOD=QA TRMIN=0.000010", 1)
    raise ValueError("unregistered optimizer recovery")


def read_molecule(path, name):
    molecules = [molecule for molecule in Chem.SDMolSupplier(str(path), removeHs=False)
                 if molecule is not None and molecule.GetProp("_Name").strip() == name]
    if len(molecules) != 1:
        raise ValueError(f"unique chemical identity required: {name}")
    return molecules[0]


def make_case(target, name, source, sdf, runtime, baseline, parameters):
    source = Path(source)
    if digest(source) != baseline["input_sha256"] or digest(source.parent / "run/gamess.log") != baseline["log_sha256"]:
        raise ValueError(f"baseline fingerprint changed: {target}/{name}")
    molecule = read_molecule(sdf, name)
    return {"target": target, "name": name, "source": str(source), "source_sha256": digest(source),
            "source_log_sha256": baseline["log_sha256"], "baseline_tpie": baseline["tpie"],
            "molblock": Chem.MolToMolBlock(molecule), "parameters": str(parameters),
            "runtime": {key: runtime[key] for key in ("gamess_dir", "library_dir", "ligand_fragment")}}


def prepare(args):
    workspace = args.workspace.resolve()
    registry = json.loads(args.registry.read_text())
    baseline = workspace / "FMOtpiefinal_HIV1lab/runs/all_accelerated_robust_20261020"
    sys.path.insert(0, str(workspace / "src"))
    from sophosqm_baseline.relative_total_energy import build_ligand_optimization_deck
    cases = {}
    for target, cohort in registry["cohorts"].items():
        for name in cohort["names"]:
            case = make_case(target, name, cohort["sources"][name], cohort["sdf"],
                             cohort["runtime"][name], cohort["baseline"][name], cohort["parameters"])
            if target in REFERENCE_RUNS:
                root = workspace / "FMOtpiefinal_HIV1lab/runs" / REFERENCE_RUNS[target]
                protocol = root / ("recovered_protocol.json" if target == "BACE" else "frozen_protocol.json")
                prior = json.loads(protocol.read_text())
                evaluation_path = root / "evaluation.json"
                evaluation = json.loads(evaluation_path.read_text())
                if evaluation["status"] != "COMPLETE_VERIFIED":
                    raise ValueError("reference experiment is incomplete")
                case["reuse"] = {state: {"task": prior["tasks"][prior["case_tasks"][name][state]],
                                         "row": evaluation["rows"][prior["case_tasks"][name][state]]}
                                 for state in ("complex", "ligand", "free_ligand_minimum")}
                case["reuse_evidence"] = {str(protocol): digest(protocol), str(evaluation_path): digest(evaluation_path)}
            cases[f"{target}/{name}"] = case
    registered = json.loads((baseline / "frozen_protocol.json").read_text())
    for key, job in registered["jobs"].items():
        if not key.startswith("OpenBind/"):
            continue
        root = baseline / key
        verification = json.loads((root / "independent_verification.json").read_text())
        if verification["status"] != "COMPLETE_INDEPENDENTLY_VERIFIED" or verification["n_verified"] != 1:
            raise ValueError("OpenBind complete reference required")
        record = verification["records"][0]
        name = record["ligand"]
        matches = [path for path in (root / "run").rglob("input.inp")
                   if path.parent.name == directory_key(name) and path.parent.parent.name == "qm"
                   and digest(path) == record["input_sha256"]]
        if not matches:
            raise ValueError(f"verified OpenBind input missing: {name}")
        source = min(matches, key=lambda path: len(str(path)))
        sdf = source.parents[2] / "ligands_no_properties.sdf"
        if not sdf.exists():
            sdf = root / "run/qm/full_target_fmo/ligands_no_properties.sdf"
        runtime = json.loads((source.parent / "run/result.json").read_text())
        parameters = (source.parent / "params").resolve(strict=True)
        cases[f"OpenBind/{name}"] = make_case("OpenBind", name, source, sdf, runtime, record, parameters)
    for key, case in cases.items():
        states = cycle_decks(Path(case["source"]).read_text(), Path(case["parameters"]))
        states.pop("receptor")
        states["free_ligand_minimum"] = build_ligand_optimization_deck(states["ligand"], nstep=200, opttol=.0001)
        molecule = Chem.MolFromMolBlock(case["molblock"], removeHs=False)
        rows = group(states["ligand"], "FMOXYZ")[1].strip().splitlines()
        if [int(float(row.split()[1])) for row in rows] != [atom.GetAtomicNum() for atom in molecule.GetAtoms()]:
            raise ValueError(f"molecular atom inventory differs from ligand deck: {key}")
        case["tasks"] = {}
        for state, deck in states.items():
            parent = args.output.resolve() / case["target"] / directory_key(case["name"]) / state / "initial"
            if "reuse" in case:
                task = case["reuse"][state]["task"]
                old = Path(task["input"]).read_text()
                allowed = ([optimization_recovery(deck, profile) for profile in RECOVERY_POLICY["optimization"]]
                           if state == "free_ligand_minimum" else [deck])
                if old not in allowed:
                    raise ValueError(f"reference task does not match frozen transfer model: {key}/{state}")
            else:
                parent.mkdir(parents=True, exist_ok=True)
                link = parent.parent / "params"
                if not link.exists():
                    link.symlink_to(Path(case["parameters"]), target_is_directory=True)
                path = parent / "input.inp"
                if path.exists() and path.read_text() != deck:
                    raise ValueError("refusing changed frozen deck")
                path.write_text(deck)
                task = {"key": key + "/" + state, "state": state, "input": str(path),
                        "input_sha256": digest(path), "runtime": case["runtime"]}
                if state == "free_ligand_minimum":
                    task["bound_ligand_deck"] = states["ligand"]
            case["tasks"][state] = task
    plan = {"schema": "all-cohort-local-ligand-strain-transfer-v1", "created_utc": datetime.now(timezone.utc).isoformat(),
            "cases": cases, "n": len(cases), "targets": dict(Counter(case["target"] for case in cases.values())),
            "source_registry_sha256": digest(args.registry), "driver_sha256": digest(__file__),
            "reuse_n": sum("reuse" in case for case in cases.values()), "recovery": RECOVERY_POLICY,
            "weights": {"tpie": 1., "local_ligand_strain": 1.}, "primary": "matched_tpie_plus_ligand_strain",
            "model": "DFTB3/D3/PCM, IDISP=0, full AB FMO2 PIEDA; L-bound and one local L-free FMO1 at the same chemical state",
            "scope": "Reuse verified prepared/k10 coordinates, no new H reconstruction or structural model selection; not a raw-PDB preparation rerun or full free energy",
            "labels_used_for_generation": False, "production_default_changed": False,
            "openbind_aggregation": "Fixed arithmetic mean over all registered crystals of each compound; never choose the best crystal"}
    if plan["n"] != 909 or plan["targets"].get("OpenBind") != 639 or plan["reuse_n"] != 72:
        raise ValueError("expected all 909 structures and 72 reusable references")
    write_json(args.output / "frozen_protocol.json", plan)
    return plan


def verify_state(task, row, case, workspace):
    sys.path.insert(0, str(Path(workspace) / "src"))
    from sophosqm_baseline.relative_total_energy import parse_converged_ligand_optimization, parse_full_energy
    path = Path(task["input"])
    log = path.parent / "gamess.log"
    if (row["status"] != "SUCCEEDED" or digest(path) != row["input_sha256"] or digest(log) != row["log_sha256"]
            or digest(Path(case["runtime"]["gamess_dir"]) / "gamess.00.x") != row["executable_sha256"]):
        raise ValueError("input/log/executable changed or QM incomplete")
    deck, text = path.read_text(), log.read_text()
    if task["state"] == "free_ligand_minimum":
        coordinates, energy = parse_converged_ligand_optimization(text, ligand_deck=task["bound_ligand_deck"])
        molecule = Chem.MolFromMolBlock(case["molblock"], removeHs=False)
        bound = np.asarray([line.split()[2:] for line in group(task["bound_ligand_deck"], "FMOXYZ")[1].strip().splitlines()], float)
        screen = geometry_screen(molecule, np.asarray([line.split()[2:] for line in coordinates], float), stereochemistry(molecule, bound))
        if not screen["passed_screen"]:
            raise ValueError("optimized ligand changed stereochemistry, bond integrity or H parent")
    else:
        energy = parse_full_energy(text, input_text=deck)
    if abs(energy - row["energy_hartree"]) > 1.e-8:
        raise ValueError("independent energy extraction differs")
    if task["state"] == "complex":
        final_fmo_convergence(deck, io.StringIO(text))
        ligand = case["runtime"]["ligand_fragment"]
        pairs = [pair for pair in parse_pie(log) if ligand in (pair.fragment_i, pair.fragment_j)]
        if len(pairs) != ligand - 1 or abs(sum(pair.value_kcal_mol for pair in pairs) - row["tpie"]) > 1.e-8:
            raise ValueError("ligand PIEDA inventory or total differs")
        if charge_aware_pair_outliers(pairs, deck, 100.)[0]:
            raise ValueError("PIEDA pair quality screen failed")
    return {"task": task, "row": row}


def worker(payload):
    key, case, workspace, output = payload
    parent = Path(output) / case["target"] / directory_key(case["name"])
    attempts, accepted = [], {}
    try:
        if digest(case["source"]) != case["source_sha256"]:
            raise ValueError("registered source changed")
        for path, expected in case.get("reuse_evidence", {}).items():
            if digest(path) != expected:
                raise ValueError("reference protocol or evaluation changed")
        for state, initial in case["tasks"].items():
            if "reuse" in case:
                accepted[state] = verify_state(initial, case["reuse"][state]["row"], case, workspace)
                continue
            deck = Path(initial["input"]).read_text()
            profiles = RECOVERY_POLICY["optimization" if state == "free_ligand_minimum" else "single_point"]
            for index, profile in enumerate(profiles):
                task = dict(initial)
                if index:
                    rendered = (optimization_recovery(deck, profile) if state == "free_ligand_minimum"
                                else recovery_deck(deck, False, profile))
                    directory = Path(initial["input"]).parent.parent / profile
                    directory.mkdir(parents=True, exist_ok=True)
                    path = directory / "input.inp"
                    if path.exists() and path.read_text() != rendered:
                        raise ValueError("frozen numerical recovery changed")
                    path.write_text(rendered)
                    task.update(input=str(path), input_sha256=digest(path))
                try:
                    _, row = execute((task, workspace))
                    accepted[state] = verify_state(task, row, case, workspace)
                    attempts.append({"state": state, "profile": profile, "status": "SUCCEEDED"})
                    break
                except Exception as error:
                    attempts.append({"state": state, "profile": profile, "status": "FAILED", "error": str(error)})
            if state not in accepted:
                raise ValueError(f"all registered strict attempts failed for {state}")
        strain = (accepted["ligand"]["row"]["energy_hartree"] - accepted["free_ligand_minimum"]["row"]["energy_hartree"]) * 627.5094740631
        tpie = accepted["complex"]["row"]["tpie"]
        result = {"status": "SUCCEEDED_INDEPENDENTLY_VERIFIED", "target": case["target"], "name": case["name"],
                  "baseline_tpie": case["baseline_tpie"], "matched_tpie": tpie, "local_ligand_strain": strain,
                  "matched_tpie_plus_ligand_strain": tpie + strain, "states": accepted, "attempts": attempts,
                  "reused_reference": "reuse" in case}
    except Exception as error:
        result = {"status": "FAILED", "target": case["target"], "name": case["name"], "error": str(error), "attempts": attempts}
    write_json(parent / "case_result.json", result)
    return key, {field: value for field, value in result.items() if field != "states"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=20)
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("positive worker count required")
    path = args.output / "frozen_protocol.json"
    frozen = json.loads(path.read_text()) if path.exists() else prepare(args)
    if digest(__file__) != frozen["driver_sha256"] or digest(args.registry) != frozen["source_registry_sha256"]:
        raise ValueError("frozen implementation or registry changed")
    if args.prepare_only:
        print(json.dumps({key: frozen[key] for key in ("n", "targets", "reuse_n")}), flush=True)
        return
    started, rows = time.time(), {}
    ordered = sorted(frozen["cases"], key=lambda key: (key.startswith("OpenBind/"), "reuse" not in frozen["cases"][key], key))
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context("spawn")) as pool:
        futures = [pool.submit(worker, (key, frozen["cases"][key], str(args.workspace.resolve()), str(args.output.resolve()))) for key in ordered]
        for future in as_completed(futures):
            key, row = future.result()
            rows[key] = row
            counts = {target: {"expected": expected, "finished": sum(value["target"] == target for value in rows.values()),
                              "succeeded": sum(value["target"] == target and value["status"] == "SUCCEEDED_INDEPENDENTLY_VERIFIED" for value in rows.values())}
                      for target, expected in frozen["targets"].items()}
            write_json(args.output / "progress.json", {"status": "RUNNING", "n_expected": frozen["n"], "n_finished": len(rows),
                "n_succeeded": sum(value["status"] == "SUCCEEDED_INDEPENDENTLY_VERIFIED" for value in rows.values()), "targets": counts, "wall_seconds": time.time() - started})
            write_json(args.output / "case_summaries.json", rows)
            print(key, row["status"], row.get("error", ""), flush=True)
    success = all(row["status"] == "SUCCEEDED_INDEPENDENTLY_VERIFIED" for row in rows.values())
    write_json(args.output / "final_status.json", {"status": "COMPLETE_INDEPENDENTLY_VERIFIED" if success else "INCOMPLETE", "n": len(rows),
               "wall_seconds": time.time() - started, "production_default_changed": False})


if __name__ == "__main__":
    main()
