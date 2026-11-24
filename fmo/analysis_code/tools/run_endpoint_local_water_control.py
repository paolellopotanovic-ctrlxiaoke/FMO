"""Matched fixed-count sampled-water control; not a default affinity score."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import time

import numpy as np

from analyze_cdk8_geometry_evidence import fragments, xyz
from analyze_openfe_endpoints import minimum_image, read_topology
from report_dense_endpoint_complex_tpie import is_water_fragment
from run_cdk8_bound_qm import replace_coordinates
from run_cdk8_solution_cycle import array, group
from run_endpoint_complex_tpie import execute
from run_matched_openfe import digest, write_json


def water_groups(topology):
    groups = {}
    covered = set()
    expected = {index for index, atom in enumerate(topology)
                if atom["class"] == 0 and atom["resname"] in {"HOH", "WAT", "TIP"}}
    for index in sorted(expected):
        if topology[index]["element"].upper() != "O":
            continue
        indices = [index, index + 1, index + 2]
        if not set(indices) <= expected or any(topology[other]["residue"] != topology[index]["residue"]
                                               for other in indices) or [topology[other]["element"].upper()
                                               for other in indices] != ["O", "H", "H"]:
            raise ValueError("invalid contiguous source water topology")
        groups[f'{topology[index]["residue"]}@{index}'] = indices
        covered.update(indices)
    if covered != expected:
        raise ValueError("incomplete source water topology")
    return groups


def select_water_geometry(coordinates, groups, ligand_polar, box, count):
    if count < 1 or len(groups) < count or len(ligand_polar) < 1:
        raise ValueError("complete nonempty water and ligand inventory required")
    ordered = sorted(groups)
    oxygen = np.asarray([coordinates[groups[key][0]] for key in ordered])
    displacement = minimum_image(oxygen[:, None] - ligand_polar[None], box)
    squared = np.sum(displacement ** 2, axis=-1)
    distances = squared.min(axis=1)
    chosen = sorted(range(len(ordered)), key=lambda index: (distances[index], ordered[index]))[:count]
    result = []
    for index in chosen:
        key = ordered[index]
        indices = groups[key]
        nearest = int(squared[index].argmin())
        origin = ligand_polar[nearest] + displacement[index, nearest]
        positions = np.array([origin] + [origin + minimum_image(coordinates[atom] - coordinates[indices[0]], box)
                                       for atom in indices[1:]])
        lengths = np.linalg.norm(positions[1:] - positions[0], axis=1)
        if len(indices) != 3 or not np.all((lengths > .7) & (lengths < 1.3)):
            raise ValueError("invalid sampled OHH water geometry")
        result.append({"saved_indices_OHH": indices, "residue": key,
                       "minimum_polar_distance_angstrom": float(np.sqrt(distances[index])),
                       "xyz_angstrom": positions.tolist()})
    return result


def main():
    from netCDF4 import Dataset

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--timeout", type=int, default=1800)
    args = parser.parse_args()
    if min(args.workers, args.timeout) < 1:
        parser.error("positive budgets required")
    started = time.monotonic()
    original_plan_path = args.source.resolve() / "frozen_protocol.json"
    original_plan = json.loads(original_plan_path.read_text())
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    tasks = []
    topology_cache = {}
    for task in original_plan["tasks"]:
        if task["iteration"] not in (400, 1200, 2000):
            continue
        source_path = Path(task["directory"]) / "input.inp"
        if digest(source_path) != task["input_sha256"]:
            raise ValueError("same-frame source changed")
        text = source_path.read_text()
        names, members = fragments(text)
        symbols, positions = xyz(text)
        charges = array(group(text, "FMO")[1], "ICHARG", "FRGNAM")
        water_slots = []
        for fragment, name in enumerate(names):
            if not is_water_fragment(name):
                continue
            atoms = [index - 1 for index in members[fragment]]
            if sorted(symbols[index] for index in atoms) != ["H", "H", "O"] or int(charges[fragment]) != 0:
                raise ValueError("only validated neutral OHH water fragments supported")
            water_slots.append((fragment, atoms))
        ligand_polar = positions[[index - 1 for index in members[-1] if symbols[index - 1] in {"N", "O", "S"}]]
        trajectory = task["trajectory"]
        if trajectory not in topology_cache:
            setup = sorted(Path(trajectory).parents[1].glob("*/hybrid_system.pdb"))[0]
            topology = read_topology(setup)
            topology_cache[trajectory] = water_groups(topology)
        with Dataset(trajectory) as data:
            iteration, replica = task["iteration"], task["replica"]
            if int(data["states"][iteration, replica]) != task["endpoint"]:
                raise ValueError("physical endpoint changed")
            raw = data["positions"][iteration, replica]
            if np.ma.getmaskarray(raw).any():
                raise ValueError("missing frame")
            waters = select_water_geometry(np.asarray(raw, float) * 10, topology_cache[trajectory], ligand_polar,
                                           np.asarray(data["box_vectors"][iteration, replica], float) * 10, len(water_slots))
        changed = positions.copy()
        slot_mapping = []
        for (fragment, atoms), water in zip(water_slots, waters, strict=True):
            water_positions = iter(water["xyz_angstrom"][1:])
            for index in atoms:
                changed[index] = water["xyz_angstrom"][0] if symbols[index] == "O" else next(water_positions)
            slot_mapping.append({"fragment_id": fragment + 1, "slot_label_not_molecule_identity": names[fragment], **water})
        rows = group(text, "FMOXYZ")[1].strip().splitlines()
        deck = replace_coordinates(text, [f"{row.split()[0]} {row.split()[1]} " + " ".join(f"{value:.10f}" for value in point)
                                         for row, point in zip(rows, changed, strict=True)])
        ligand_atoms = [index - 1 for index in members[-1]]
        receptor_atoms = [index - 1 for fragment in members[:-1] for index in fragment]
        contacts = np.linalg.norm(changed[receptor_atoms, None] - changed[None, ligand_atoms], axis=-1)
        receptor_heavy = [index for index, atom in enumerate(receptor_atoms) if symbols[atom] != "H"]
        ligand_heavy = [index for index, atom in enumerate(ligand_atoms) if symbols[atom] != "H"]
        geometry = {"minimum_all_atom_contact_angstrom": float(contacts.min()),
                    "minimum_heavy_contact_angstrom": float(contacts[np.ix_(receptor_heavy, ligand_heavy)].min())}
        directory = output / task["target"] / task["name"] / Path(trajectory).parent.name / f'iteration_{iteration}'
        directory.mkdir(parents=True)
        (directory / "input.inp").write_text(deck)
        write_json(directory / "water_mapping.json", slot_mapping)
        tasks.append({**task, "kind": "same_frame_fixed_count_local_water_control", "directory": str(directory),
                      "source": str(source_path), "source_sha256": digest(source_path),
                      "input_sha256": digest(directory / "input.inp"), "geometry": geometry})
    if len(tasks) != 36:
        raise ValueError("fixed 36-frame matched schedule required")
    plan = {"tasks": tasks, "fixed_iterations": [400, 1200, 2000], "driver_sha256": digest(__file__),
        "source_protocol_sha256": digest(original_plan_path), "labels_used": False, "main_pipeline_changed": False,
        "scope": "Matched same frames; only water slots replaced by fixed-count closest actual MD waters to ligand N/O/S. Same protein, ligand, caps, PCM and fragment count; not a binding free energy or default water policy."}
    write_json(output / "frozen_protocol.json", plan)
    results = []
    write_json(output / "status.json", {"status": "RUNNING", "n_requested": len(tasks)})
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(execute, task, args.timeout): task for task in tasks}
        for future in as_completed(futures):
            task = futures[future]
            try:
                row = future.result()
            except Exception as error:
                row = {**task, "result": {"status": "FAILED", "error": f"{type(error).__name__}: {error}"}}
                write_json(Path(task["directory"]) / "summary.json", row)
            results.append(row)
            write_json(output / "progress.json", results)
            print(len(results), len(tasks), task["target"], task["name"], row["result"]["status"], flush=True)
    succeeded = sum(row["result"]["status"] == "SUCCEEDED" for row in results)
    write_json(output / "status.json", {"status": "COMPLETE" if succeeded == len(tasks) else "INCOMPLETE_WITH_FAILURES",
        "n_requested": len(tasks), "n_succeeded": succeeded, "wall_seconds_including_preflight": time.monotonic() - started,
        "cohort_ranking_available": False})


if __name__ == "__main__":
    main()
