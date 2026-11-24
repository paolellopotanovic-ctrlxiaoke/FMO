"""Full-cohort fixed-heavy QM hydrogen normalization; diagnostic, not free energy."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import multiprocessing
from pathlib import Path
import re
import sys
import time

import numpy as np
from rdkit import Chem

from audit_cycle_ligand_chemistry import geometry_screen, stereochemistry
from fmotpiefinal.core.parser import parse_pie
from fmotpiefinal.core.runner import charge_aware_pair_outliers
from fmotpiefinal.core.solution_cycle import HARTREE_TO_KCAL_MOL
from run_cdk8_bound_qm import replace_coordinates
from run_cdk8_solution_cycle import digest, group, metrics, write_json
from run_cohort_solution_cycle import execute
from run_local_sidechain_transfer import directory_key


def xyz_rows(deck):
    return group(deck, "FMOXYZ")[1].strip().splitlines()


def hydrogen_only_deck(deck):
    if ("RUNTYP=ENERGY" not in deck or not re.search(r"NFRAG=1\s+NBODY=1", deck)
            or re.search(r"\b(?:MODGRD|NZVAR|IACTAT|IFREEZ)\s*=|\$STATPT", deck)):
        raise ValueError("unconstrained single-fragment energy deck required")
    active = [str(index+1) for index, row in enumerate(xyz_rows(deck)) if float(row.split()[1]) == 1.]
    if not active:
        raise ValueError("no hydrogen atoms to optimize")
    result = deck.replace("RUNTYP=ENERGY", "RUNTYP=OPTIMIZE", 1)
    control = group(result, "CONTRL")
    replacement = control[0].replace("$END", "\n NZVAR=0 $END")
    result = result[:control.start()]+replacement+result[control.end():]
    result = result.replace("NFRAG=1 NBODY=1 NACUT=0", "NFRAG=1 NBODY=1 NACUT=0 MODGRD=32", 1)
    active_lines, current = [], " IACTAT(1)="
    for index in active:
        if len(current)+len(index)+1 > 75:
            active_lines.append(current)
            current = " "
        current += index+","
    active_lines.append(current.rstrip(","))
    result += ("\n $STATPT METHOD=RFO NSTEP=500 OPTTOL=0.000100 TRMIN=0.001\n"
               +"\n".join(active_lines)+"\n $END\n")
    if "MODGRD=32" not in result or any(len(line) > 80 for line in result.splitlines()):
        raise ValueError("invalid constrained optimization deck")
    return result


def fixed_heavy_displacement(before, after):
    if len(before) != len(after):
        raise ValueError("atom inventory changed")
    distances = []
    for initial, final in zip(before, after, strict=True):
        left, right = initial.split(), final.split()
        if left[0].upper() != right[0].upper() or float(left[1]) != float(right[1]):
            raise ValueError("atom order or elements changed")
        if float(left[1]) != 1.:
            distances.append(float(np.linalg.norm(np.asarray(left[2:], float)-np.asarray(right[2:], float))))
    maximum = max(distances, default=0.)
    if maximum > 1.e-5:
        raise ValueError("fixed ligand heavy atoms moved")
    return maximum


def transplant_ligand(complex_deck, initial_ligand, final_rows):
    original = xyz_rows(complex_deck)
    ligand_rows = xyz_rows(initial_ligand)
    fixed_heavy_displacement(ligand_rows, final_rows)
    if original[-len(ligand_rows):] != ligand_rows:
        raise ValueError("ligand is not the exact complex suffix")
    return replace_coordinates(complex_deck, original[:-len(ligand_rows)]+final_rows)


def worker(payload):
    item, workspace = payload
    parent = Path(item["input"]).parent.parent
    try:
        key, optimized = execute((item, workspace))
        if optimized["status"] != "SUCCEEDED":
            raise ValueError(optimized.get("error", "H optimization failed"))
        rows = optimized["optimized_coordinates"]
        shift = fixed_heavy_displacement(xyz_rows(item["bound_ligand_deck"]), rows)
        molecule = Chem.MolFromMolBlock(item["molblock"], removeHs=False)
        initial = np.asarray([row.split()[2:] for row in xyz_rows(item["bound_ligand_deck"])], float)
        screen = geometry_screen(molecule, np.asarray([row.split()[2:] for row in rows], float), stereochemistry(molecule, initial))
        if not screen["passed_screen"]:
            raise ValueError("constrained QM geometry failed chemical screen")
        result = {"status": "HYDROGEN_VERIFIED", "hydrogen": optimized, "chemistry": screen,
                  "max_heavy_displacement_A": shift}
        write_json(parent/"hydrogen_summary.json", result)
        template_path = Path(item["complex_template"])
        if digest(template_path) != item["complex_template_sha256"]:
            raise ValueError("source complex changed")
        deck = transplant_ligand(template_path.read_text(), item["bound_ligand_deck"], rows)
        path = parent/"complex/input.inp"
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists() and path.read_text() != deck:
            raise ValueError("refusing to change the generated complex")
        path.write_text(deck)
        task = {"key": key+"/complex", "state": "complex", "input": str(path),
                "input_sha256": digest(path), "runtime": item["runtime"]}
        write_json(parent/"generated_complex_protocol.json", task)
        _, complex_result = execute((task, workspace))
        if complex_result["status"] != "SUCCEEDED":
            raise ValueError(complex_result.get("error", "complex failed"))
        fragment = int(re.search(r"NFRAG=(\d+)", deck)[1])
        pairs = [pair for pair in parse_pie(path.parent/"gamess.log") if fragment in (pair.fragment_i, pair.fragment_j)]
        outliers, ionic = charge_aware_pair_outliers(pairs, deck, 100.)
        if outliers:
            raise ValueError("changed complex failed original charge-aware PIEDA gate")
        result.update(status="SUCCEEDED", complex=complex_result, strong_ionic_pair_count=len(ionic))
    except Exception as error:
        result = {"status": "FAILED", "error": f"{type(error).__name__}: {error}"}
    write_json(parent/"summary.json", result)
    return item["name"], result


def evaluate(rows, frozen, prior, truth):
    if set(truth) != set(frozen["names"]) or set(rows) != set(frozen["names"]):
        raise ValueError("full matched cohort inventory required")
    complete = all(row["status"] == "SUCCEEDED" for row in rows.values())
    result = {"status": "COMPLETE_VERIFIED" if complete else "INCOMPLETE", "rows": rows,
              "n_expected": len(frozen["names"]), "n_succeeded": sum(row["status"] == "SUCCEEDED" for row in rows.values()),
              "production_default_changed": False, "scores": {}, "metrics": {}}
    if not complete:
        return result
    for name in frozen["names"]:
        old, row = prior["scores"][name], rows[name]
        bound = prior["rows"][frozen["source_case_tasks"][name]["ligand"]]["energy_hartree"]
        hydrogen_strain = (bound-row["hydrogen"]["energy_hartree"])*HARTREE_TO_KCAL_MOL
        residual = old["ligand_strain"]-hydrogen_strain
        if hydrogen_strain < -.01:
            raise ValueError("H optimization raised energy beyond printed precision")
        tpie = row["complex"]["tpie"]
        result["scores"][name] = {"matched_tpie": old["matched_tpie"],
            "original_total_strain_score": old["tpie_plus_ligand_strain"],
            "hydrogen_relaxation_energy": hydrogen_strain, "residual_local_strain": residual,
            "original_tpie_plus_hydrogen_component": old["matched_tpie"]+hydrogen_strain,
            "original_tpie_plus_residual_strain": old["matched_tpie"]+residual,
            "hydrogen_normalized_tpie": tpie,
            "hydrogen_normalized_tpie_plus_residual_strain": tpie+residual}
    names = sorted(rows)
    for arm in frozen["arms"]:
        result["metrics"][arm] = metrics([truth[name] for name in names], [-result["scores"][name][arm] for name in names])
    result["warning"] = "Isolated-ligand H optimization is not bound-state optimization. Residual strain can depend on local-minimum basins; no entropy or free-energy guarantee."
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--evaluation", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=18)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("positive workers required")
    start = time.monotonic()
    protocol = json.loads(args.protocol.read_text())
    prior = json.loads(args.evaluation.read_text())
    cohort = json.loads(args.registry.read_text())["cohorts"][protocol["target"]]
    names = protocol["names"]
    if prior["status"] != "COMPLETE_VERIFIED" or set(prior["scores"]) != set(names) or set(cohort["names"]) != set(names):
        raise ValueError("complete verified source required")
    for path, fingerprint in cohort["input_sha256"].items():
        if digest(path) != fingerprint:
            raise ValueError("registered evidence changed")
    molecules = {molecule.GetProp("_Name").strip(): molecule for molecule in Chem.SDMolSupplier(cohort["sdf"], removeHs=False) if molecule is not None}
    if set(molecules) != set(names):
        raise ValueError("complete chemical inventory required")
    output = args.output.resolve()
    tasks = {}
    for name in names:
        source_key = protocol["case_tasks"][name]["ligand"]
        source = protocol["tasks"][source_key]
        ligand_path = Path(source["input"])
        if digest(ligand_path) != prior["rows"][source_key]["input_sha256"]:
            raise ValueError("source ligand changed")
        complex_key = protocol["case_tasks"][name]["complex"]
        complex_path = Path(protocol["tasks"][complex_key]["input"])
        if digest(complex_path) != prior["rows"][complex_key]["input_sha256"]:
            raise ValueError("source complex changed")
        parent = output/protocol["target"]/directory_key(name)
        parent.mkdir(parents=True, exist_ok=True)
        link = parent/"params"
        if link.exists() and link.resolve() != Path(protocol["parameters"]).resolve():
            raise ValueError("parameter link mismatch")
        if not link.exists():
            link.symlink_to(Path(protocol["parameters"]).resolve(), target_is_directory=True)
        path = parent/"hydrogen/input.inp"
        path.parent.mkdir(parents=True, exist_ok=True)
        deck = hydrogen_only_deck(ligand_path.read_text())
        if path.exists() and path.read_text() != deck:
            raise ValueError("refusing changed optimization input")
        path.write_text(deck)
        tasks[name] = {"name": name, "key": name+"/hydrogen", "state": "free_ligand_minimum",
            "input": str(path), "input_sha256": digest(path), "bound_ligand_deck": ligand_path.read_text(),
            "complex_template": str(complex_path), "complex_template_sha256": digest(complex_path),
            "runtime": source["runtime"], "molblock": Chem.MolToMolBlock(molecules[name])}
    frozen = {"schema": "fixed-heavy-hydrogen-normalized-strain-v1", "target": protocol["target"], "names": names,
        "tasks": tasks, "source_case_tasks": protocol["case_tasks"], "driver_sha256": digest(__file__),
        "dependencies": {name: digest(Path(__file__).with_name(name)) for name in ("run_cohort_solution_cycle.py", "audit_cycle_ligand_chemistry.py")},
        "protocol_sha256": digest(args.protocol), "evaluation_sha256": digest(args.evaluation),
        "registry_sha256": digest(args.registry), "labels_used_for_generation": False,
        "arms": ["matched_tpie", "original_total_strain_score", "original_tpie_plus_hydrogen_component",
                 "original_tpie_plus_residual_strain", "hydrogen_normalized_tpie", "hydrogen_normalized_tpie_plus_residual_strain"],
        "controls": "Same source k10 heavy coordinates, receptor/waters, charge, DFTB3/D3 and PCM IDISP0. Isolated H-only RFO500, OPTTOL1e-4, strict original SCF. No weight scan or label-dependent recovery.",
        "limitations": "Diagnostic not default: isolated optimized H may be incompatible with binding contacts; all failures retained; no partial-cohort metrics"}
    path = output/"frozen_protocol.json"
    if path.exists() and json.loads(path.read_text()) != frozen:
        raise ValueError("frozen experiment changed")
    write_json(path, frozen)
    rows = {}
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context("spawn")) as pool:
        futures = [pool.submit(worker, (task, str(args.workspace.resolve()))) for task in tasks.values()]
        for future in as_completed(futures):
            name, row = future.result()
            rows[name] = row
            write_json(output/"progress.json", {"n_expected": len(names), "n_finished": len(rows),
                "n_succeeded": sum(value["status"] == "SUCCEEDED" for value in rows.values()), "rows": rows})
            print(name, row["status"], row.get("error", ""), flush=True)
    from report_local_sidechain_transfer import load_truth
    result = evaluate(rows, frozen, prior, load_truth(args.workspace.resolve(), protocol["target"]))
    result["wall_seconds"] = time.monotonic()-start
    write_json(output/"evaluation.json", result)
    print(result["status"], result["metrics"], flush=True)


if __name__ == "__main__":
    main()
