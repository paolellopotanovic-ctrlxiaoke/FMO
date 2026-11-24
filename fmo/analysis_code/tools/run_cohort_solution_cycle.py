"""Matched same-geometry PCM three-state and free-ligand controls, without fitting."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import io
import json
import multiprocessing
from pathlib import Path
import re
import sys
import time

from fmotpiefinal.core.fmo_convergence import final_fmo_convergence
from fmotpiefinal.core.parser import parse_pie
from fmotpiefinal.core.solution_cycle import require_consistent_pcm_cycle
from run_cdk8_solution_cycle import digest, group, metrics, split_decks, write_json
from run_local_sidechain_transfer import directory_key


def cycle_decks(source, parameters):
    source = re.sub(r"\bIDISP=1\b", "IDISP=0", source)
    states = split_decks(source, parameters)
    states["complex"] = states["complex"].replace("IPIEDA=0", "IPIEDA=1")
    return states


def canonical_receptor(deck):
    deck = re.sub(r"^!.*\n", "! Matched common receptor\n", deck, count=1)
    data = group(deck, "DATA")
    lines = data[1].strip().splitlines()
    lines[0] = "Matched common receptor"
    return deck[:data.start(1)]+"\n".join(lines)+"\n "+deck[data.end(1):]


def execute(payload):
    task, workspace = payload
    sys.path.insert(0, str(Path(workspace) / "src"))
    from sophosqm_baseline.relative_total_energy import execute_deck, parse_converged_ligand_optimization
    path = Path(task["input"])
    try:
        if digest(path) != task["input_sha256"]:
            raise ValueError("frozen deck changed")
        binary = Path(task["runtime"]["gamess_dir"]) / "gamess.00.x"
        summary = path.parent / "verified_summary.json"
        if summary.exists():
            cached = json.loads(summary.read_text())
            if (cached["status"] == "SUCCEEDED" and cached["input_sha256"] == digest(path)
                    and cached["executable_sha256"] == digest(binary)
                    and cached["log_sha256"] == digest(path.parent / "gamess.log")):
                return task["key"], cached
        result = execute_deck(path.read_text(), path.parent, gamess_dir=task["runtime"]["gamess_dir"],
            library_dir=task["runtime"]["library_dir"], timeout=2400, ncores=1)
        text = (path.parent / "gamess.log").read_text()
        optimized = None
        if task["state"] == "free_ligand_minimum":
            if result.get("returncode") != 0:
                raise ValueError("ligand optimizer process failed")
            optimized = parse_converged_ligand_optimization(text, ligand_deck=task["bound_ligand_deck"])
            result.update(status="SUCCEEDED", energy_hartree=optimized[1])
        elif result["status"] != "SUCCEEDED":
            raise ValueError(f"QM did not pass: {result['status']}: {result.get('error')}")
        row = {"status": "SUCCEEDED", "energy_hartree": result["energy_hartree"],
            "input_sha256": digest(path), "log_sha256": digest(path.parent / "gamess.log"),
            "executable_sha256": digest(binary),
            "wall_seconds": result["wall_seconds"]}
        if task["state"] in {"complex", "receptor"}:
            row["convergence"] = final_fmo_convergence(path.read_text(), io.StringIO(text))
        if task["state"] == "complex":
            count = int(re.search(r"NFRAG=(\d+)", path.read_text())[1])
            pairs = [pair for pair in parse_pie(path.parent / "gamess.log")
                     if count in (pair.fragment_i, pair.fragment_j)]
            if len(pairs) != count-1:
                raise ValueError("full AB ligand-pair count mismatch")
            row["tpie"] = sum(pair.value_kcal_mol for pair in pairs)
        if optimized is not None:
            row.update(energy_hartree=optimized[1], optimized_coordinates=optimized[0])
    except Exception as error:
        row = {"status": "FAILED", "error": f"{type(error).__name__}: {error}"}
    write_json(path.parent / "verified_summary.json", row)
    return task["key"], row


def evaluate(rows, frozen, workspace):
    from report_local_sidechain_transfer import load_truth
    names = sorted(frozen["names"])
    result = {"rows": rows, "n_expected": len(frozen["tasks"]),
        "n_succeeded": sum(row["status"] == "SUCCEEDED" for row in rows.values()),
        "status": "COMPLETE_VERIFIED" if all(row["status"] == "SUCCEEDED" for row in rows.values()) else "INCOMPLETE",
        "production_default_changed": False, "scores": {}, "metrics": {}}
    truth = load_truth(workspace, frozen["target"])
    if set(truth) != set(names):
        raise ValueError("full cohort evaluation identity mismatch")
    for name in names:
        keys = {state: frozen["case_tasks"][name][state] for state in ("complex", "receptor", "ligand")}
        if any(rows[key]["status"] != "SUCCEEDED" for key in keys.values()):
            continue
        decks = {state: Path(frozen["tasks"][key]["input"]).read_text() for state, key in keys.items()}
        logs = {state: (Path(frozen["tasks"][key]["input"]).parent / "gamess.log").read_text()
                for state, key in keys.items()}
        audit = require_consistent_pcm_cycle(decks, logs)
        energies = {state: rows[key]["energy_hartree"] for state, key in keys.items()}
        binding = (energies["complex"]-energies["receptor"]-energies["ligand"])*627.5094740631
        cavity = audit["states"]["complex"]["cavitation"]-audit["states"]["receptor"]["cavitation"]-audit["states"]["ligand"]["cavitation"]
        tpie = rows[keys["complex"]]["tpie"]
        scores = {"matched_tpie": tpie, "same_geometry_three_state": binding,
                  "tpie_plus_cavity": tpie+cavity, "cavity_delta": cavity,
                  "baseline_tpie": frozen["baseline"][name], "pcm_audit": audit}
        free = rows[frozen["case_tasks"][name]["free_ligand_minimum"]]
        if free["status"] == "SUCCEEDED":
            strain = (energies["ligand"]-free["energy_hartree"])*627.5094740631
            scores.update(ligand_strain=strain, free_minimum_three_state=binding+strain,
                          tpie_plus_ligand_strain=tpie+strain)
        result["scores"][name] = scores
    for score in ("baseline_tpie", "matched_tpie", "same_geometry_three_state", "tpie_plus_cavity",
                  "free_minimum_three_state", "tpie_plus_ligand_strain"):
        if len(result["scores"]) == len(names) and all(score in result["scores"][name] for name in names):
            result["metrics"][score] = metrics([truth[name] for name in names],
                                               [-result["scores"][name][score] for name in names])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--target", default="BACE")
    parser.add_argument("--workers", type=int, default=18)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("positive workers required")
    workspace, output, source = args.workspace.resolve(), args.output.resolve(), args.source.resolve()
    sys.path.insert(0, str(workspace / "src"))
    from sophosqm_baseline.relative_total_energy import build_ligand_optimization_deck
    start = time.time()
    cohort = json.loads(source.read_text())["cohorts"][args.target]
    for path, expected in cohort["input_sha256"].items():
        if digest(path) != expected:
            raise ValueError("source evidence changed")
    parameters = Path(cohort["parameters"])
    tasks, case_tasks = {}, {}
    receptor_decks = {}
    for name in cohort["names"]:
        states = cycle_decks(Path(cohort["sources"][name]).read_text(), parameters)
        states["receptor"] = canonical_receptor(states["receptor"])
        states["free_ligand_minimum"] = build_ligand_optimization_deck(states["ligand"], nstep=200, opttol=.0001)
        case_tasks[name] = {}
        for state, deck in states.items():
            if state == "receptor" and deck in receptor_decks:
                case_tasks[name][state] = receptor_decks[deck]
                continue
            key = f"{name}/{state}"
            case_tasks[name][state] = key
            if state == "receptor":
                receptor_decks[deck] = key
            parent = output / args.target / directory_key(name) / state
            parent.mkdir(parents=True, exist_ok=True)
            link = parent.parent / "params"
            if link.is_symlink() and link.resolve() != parameters.resolve():
                raise ValueError("parameter link mismatch")
            if not link.exists():
                link.symlink_to(parameters.resolve(), target_is_directory=True)
            path = parent / "input.inp"
            if path.exists() and path.read_text() != deck:
                raise ValueError("refusing changed frozen deck")
            path.write_text(deck)
            tasks[key] = {"key": key, "name": name, "state": state, "input": str(path),
                "input_sha256": digest(path), "runtime": cohort["runtime"][name]}
            if state == "free_ligand_minimum":
                tasks[key]["bound_ligand_deck"] = states["ligand"]
    frozen = {"schema": "cohort-matched-pcm-three-state-v1", "target": args.target,
        "names": cohort["names"], "source_sha256": digest(source), "tasks": tasks, "case_tasks": case_tasks,
        "baseline": {name: cohort["baseline"][name]["tpie"] for name in cohort["names"]},
        "driver_sha256": digest(__file__), "parameters": str(parameters),
        "labels_used_for_generation": False, "production_default_changed": False,
        "matching": "Exact source atoms/caps/charges/HOPs in A, B and AB; IDISP=0 all states, same VANDW cavity and DFTB3/D3; full AB all-pair FMO2 not selective total energy",
        "free_reference": "Single PCM FMO1 local ligand minimum at unchanged chemical state; not a conformational ensemble or binding free energy",
        "arms": ["matched_tpie", "same_geometry_three_state", "tpie_plus_cavity", "free_minimum_three_state", "tpie_plus_ligand_strain"],
        "hybrid_warning": "tPIE plus ligand strain is a diagnostic unit-coefficient hybrid, not a closed free-energy cycle"}
    path = output / "frozen_protocol.json"
    if path.exists() and json.loads(path.read_text()) != frozen:
        raise ValueError("changed frozen experiment")
    write_json(path, frozen)
    write_json(output / "status.json", {"status": "RUNNING", "n": len(tasks), "n_shared_receptor_tasks": len(receptor_decks)})
    rows = {}
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context("spawn")) as pool:
        futures = [pool.submit(execute, (task, str(workspace))) for task in tasks.values()]
        for future in as_completed(futures):
            key, row = future.result()
            rows[key] = row
            write_json(output / "progress.json", {"n_expected": len(tasks), "n_finished": len(rows),
                "n_succeeded": sum(item["status"] == "SUCCEEDED" for item in rows.values()), "rows": rows})
            print(key, row["status"], row.get("energy_hartree", row.get("error")), flush=True)
    report = evaluate(rows, frozen, workspace)
    report["wall_seconds"] = time.time()-start
    write_json(output / "evaluation.json", report)
    write_json(output / "status.json", {key: value for key, value in report.items() if key not in {"rows", "scores"}})


if __name__ == "__main__":
    main()
