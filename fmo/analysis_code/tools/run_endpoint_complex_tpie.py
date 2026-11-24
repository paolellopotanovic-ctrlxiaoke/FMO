"""Evaluate matched MD complexes with fixed FMO chemistry, not a free energy."""
from __future__ import annotations

import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import replace
import json
from pathlib import Path
import re
import time

import numpy as np

from analyze_endpoint_internal_energy import physical_atom_maps, verify_bonds
from analyze_openfe_endpoints import minimum_image, read_topology, unwrap_ligand
from run_cdk8_bound_qm import replace_coordinates
from run_cdk8_solution_cycle import group
from run_matched_openfe import digest, write_json


WATERS = {"HOH", "WAT", "TIP"}


def map_receptor(atoms, topology):
    lookup = defaultdict(list)
    for index, atom in enumerate(topology):
        if atom["class"] == 0:
            lookup[(tuple(atom["residue"].split()), atom["resname"], atom["name"])].append(index)
    assigned = []
    modes = []
    for atom in atoms:
        candidates = lookup.get(((atom.chain, atom.residue_number), atom.residue_name, atom.atom_name), [])
        if len(candidates) > 1:
            raise ValueError("ambiguous receptor identity")
        if not candidates or atom.residue_name in WATERS:
            candidates = [index for index, other in enumerate(topology)
                          if other["class"] == 0 and other["element"].upper() == atom.element.upper()
                          and (other["resname"] in WATERS if atom.residue_name in WATERS
                               else atom.record == "HETATM" or other["resname"] == atom.residue_name)
                          and np.linalg.norm(np.asarray(other["xyz"]) - atom.xyz) < .003]
            mode = "verified_initial_geometry"
        else:
            mode = "residue_atom_identity"
        if len(candidates) != 1:
            raise ValueError(f"unique receptor mapping required: {atom.chain}/{atom.residue_number}/{atom.atom_name}")
        index = candidates[0]
        if topology[index]["element"].upper() != atom.element.upper():
            raise ValueError("mapped receptor element mismatch")
        assigned.append(index)
        modes.append(mode)
    if len(set(assigned)) != len(assigned):
        raise ValueError("receptor mapping is not one-to-one")
    return np.asarray(assigned, dtype=int), modes


def whole_components(coordinates, bonds, box, ligand_center, anchor_indices):
    adjacency = defaultdict(list)
    for first, second in bonds:
        adjacency[first].append(second)
        adjacency[second].append(first)
    result = np.array(coordinates, dtype=float, copy=True)
    remaining = set(range(len(coordinates)))
    anchors = set(anchor_indices)
    while remaining:
        seed = min(remaining)
        remaining.remove(seed)
        component = [seed]
        pending = [seed]
        while pending:
            first = pending.pop()
            for second in adjacency[first]:
                if second in remaining:
                    remaining.remove(second)
                    result[second] = result[first] + minimum_image(coordinates[second] - coordinates[first], box)
                    component.append(second)
                    pending.append(second)
        selected = sorted(anchors.intersection(component)) or component
        center = result[selected].mean(axis=0)
        displacement = center - ligand_center
        result[component] += minimum_image(displacement, box) - displacement
    if not np.isfinite(result).all():
        raise ValueError("nonfinite unwrapped receptor")
    for first, second in bonds:
        reference = minimum_image(coordinates[second] - coordinates[first], box)
        if not np.allclose(result[second] - result[first], reference, atol=1e-4):
            raise ValueError("periodic bond cycle inconsistent")
    return result


def deck_geometry(source, receptor, ligand, keys):
    from fmotpiefinal.core.fmo_input import build_deck_atoms
    from fmotpiefinal.core.fragmentation import make_facio_like_fragments
    from fmotpiefinal.core.structure import Protein, residue_key

    pocket = Protein(tuple(atom for atom in receptor.atoms if residue_key(atom) in keys))
    prepared, fragments, metadata = make_facio_like_fragments(pocket, ligand,
        reference_protein=receptor, use_hop=True, merge_disulfides=True, repair_peptide_numbering=True)
    deck_atoms = build_deck_atoms(prepared, ligand, fragments)
    original = group(source, "FMOXYZ")[1].strip().splitlines()
    if len(original) != len(deck_atoms) or any(row.split()[0].upper() != item.atom.element.upper()
                                               for row, item in zip(original, deck_atoms, strict=True)):
        raise ValueError("regenerated FMO atom order/inventory changed")
    signature = [(item.atom.element, item.atom.chain, item.atom.residue_number, item.atom.residue_name,
                  item.atom.atom_name, item.fragment_id) for item in deck_atoms]
    fragment_signature = [(fragment.name, fragment.kind, fragment.atom_indices, fragment.charge,
                           fragment.multiplicity) for fragment in fragments]
    rows = [f"{row.split()[0]} {row.split()[1]} " + " ".join(f"{value:.10f}" for value in item.atom.xyz)
            for row, item in zip(original, deck_atoms, strict=True)]
    return replace_coordinates(source, rows), signature, fragment_signature, metadata, deck_atoms


def geometry_contacts(deck_atoms, ligand_count):
    receptor = deck_atoms[:-ligand_count]
    ligand = deck_atoms[-ligand_count:]
    distances = np.linalg.norm(np.asarray([item.atom.xyz for item in receptor])[:, None]
                               - np.asarray([item.atom.xyz for item in ligand])[None], axis=-1)
    receptor_heavy = [index for index, item in enumerate(receptor) if item.atom.element.upper() != "H"]
    ligand_heavy = [index for index, item in enumerate(ligand) if item.atom.element.upper() != "H"]
    first, second = np.unravel_index(distances.argmin(), distances.shape)
    return {"minimum_all_atom_contact_angstrom": float(distances.min()),
            "minimum_heavy_contact_angstrom": float(distances[np.ix_(receptor_heavy, ligand_heavy)].min()),
            "closest_receptor_atom": {"name": receptor[first].atom.atom_name,
                "residue": receptor[first].atom.residue_name, "number": receptor[first].atom.residue_number},
            "closest_ligand_atom_index": int(second),
            "closest_atom_in_cap": receptor[first].atom.residue_name in {"ACE", "NME"}}


def freeze(args, *, iterations=(400, 1200, 2000), driver_path=None):
    from gufe import LigandAtomMapping
    from netCDF4 import Dataset
    from openmm import app
    from fmotpiefinal.core.bound_geometry import load_bound_ligand
    from fmotpiefinal.core.minimization import _water_compatible_pdb
    from fmotpiefinal.core.structure import Protein, read_pdb, residue_key
    from fmotpiefinal.pipeline.target_stages import receptor_at_minimized_geometry

    root = args.workspace.resolve() / "FMOtpiefinal_HIV1lab/runs"
    registry_path = root / "graph_sidechain_transfer_20261021/frozen_protocol.json"
    cohorts = json.loads(registry_path.read_text())["cohorts"]
    plan_path = args.run.resolve() / "final_endpoint_audit_20261023/qm/frozen_protocol.json"
    previous = json.loads(plan_path.read_text())
    if previous["labels_used"]:
        raise ValueError("label-independent source required")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    tasks = []
    provenance = {str(registry_path): digest(registry_path), str(plan_path): digest(plan_path)}
    audits = {}
    for target in ("BACE", "CDK8"):
        cohort = cohorts[target]
        receptor = read_pdb(cohort["protein"])
        provenance[cohort["protein"]] = digest(cohort["protein"])
        topology_root = args.run.resolve() / f"production/{target}/edge_000/bound"
        topology_paths = sorted((topology_root / "shared").glob("*/hybrid_system.pdb"))
        topology = read_topology(topology_paths[0])
        indices, modes = map_receptor(receptor.atoms, topology)
        setups = []
        signature = lambda atoms: [(atom["name"], atom["element"], atom["class"]) for atom in atoms]
        for path in topology_paths:
            other = read_topology(path)
            if signature(other) != signature(topology) or any(other[index] != topology[index] for index in indices):
                raise ValueError("repeat receptor identity or atom-index signature differs")
            setup_path = path.with_name("hybrid_positions.npy")
            setups.append((path, np.load(setup_path) * 10, np.asarray([atom["xyz"] for atom in other])))
            provenance[str(path)] = digest(path)
            provenance[str(setup_path)] = digest(setup_path)
        mapping_path = topology_root / "mapping.json"
        mapping = LigandAtomMapping.from_json(mapping_path)
        provenance[str(mapping_path)] = digest(mapping_path)
        ligand_indices, molecules = physical_atom_maps(topology, mapping)
        ligand_bonds = {endpoint: verify_bonds(topology_paths[0], selected, molecules[endpoint])
                        for endpoint, selected in ligand_indices.items()}
        original_pdb = app.PDBFile(str(_water_compatible_pdb(cohort["protein"])))
        original_atoms = list(original_pdb.topology.atoms())
        if len(original_atoms) != len(receptor.atoms) or any(atom.element.symbol.upper() != expected.element.upper()
                or atom.name != expected.atom_name for atom, expected in zip(original_atoms, receptor.atoms, strict=True)):
            raise ValueError("bond topology atom order mismatch")
        receptor_bonds = [(first.index, second.index) for first, second in original_pdb.topology.bonds()]
        keys = set(cohort["pocket_keys"])
        anchors = [index for index, atom in enumerate(receptor.atoms)
                   if atom.atom_name == "CA" and residue_key(atom) in keys]
        if len(anchors) < 3:
            raise ValueError("three pocket receptor anchors required")
        audits[target] = {"n_receptor_atoms": len(indices), "mapping_modes": modes,
                         "source_to_saved_indices": indices.tolist(), "bonds": receptor_bonds,
                         "fixed_pocket_keys": sorted(keys), "n_retained_source_waters": len({residue_key(atom)
                            for atom in receptor.atoms if atom.residue_name in WATERS}),
                         "warning": "Tracks original water molecules, not exchanging occupancy; no bulk-water insertion."}
        verified_trajectories = set()
        source_cache = {}
        for original in previous["tasks"]:
            if original["target"] != target or original["leg"] != "bound" or original["iteration"] not in iterations:
                continue
            name, endpoint = original["name"], original["endpoint"]
            source = (root / f"main_response_regression_20261021/CDK8/sidechain_response/qm/{name}/input.inp"
                      if target == "CDK8" else Path(cohort["sources"][name]))
            pose = (root / f"main_response_regression_20261021/CDK8/sidechain_response/minimized/{name}/minimized_complex.pdb"
                    if target == "CDK8" else Path(cohort["complexes"][name]))
            if name not in source_cache:
                source_receptor = receptor_at_minimized_geometry(receptor, pose) if target == "CDK8" else receptor
                ligand = load_bound_ligand(cohort["sdf"], name, pose)
                text = source.read_text()
                regenerated, signature, fragments, metadata, atoms = deck_geometry(text, source_receptor, ligand, keys)
                source_positions = np.asarray([row.split()[2:] for row in group(text, "FMOXYZ")[1].strip().splitlines()], float)
                reproduced_positions = np.asarray([item.atom.xyz for item in atoms])
                if not np.allclose(source_positions, reproduced_positions, atol=1e-7, rtol=0):
                    raise ValueError("source cap/fragment coordinates cannot be reproduced")
                for path in (source, pose, Path(cohort["sdf"])):
                    provenance[str(path)] = digest(path)
                source_cache[name] = ligand, text, signature, fragments, metadata
            ligand, text, signature, fragments, metadata = source_cache[name]
            with Dataset(original["trajectory"]) as data:
                iteration, replica = original["iteration"], original["replica"]
                if int(data["states"][iteration, replica]) != endpoint:
                    raise ValueError("physical endpoint changed")
                raw = data["positions"][iteration, replica]
                if np.ma.getmaskarray(raw).any():
                    raise ValueError("missing same-frame positions")
                if original["trajectory"] not in verified_trajectories:
                    saved_indices = np.asarray(data["analysis_particle_indices"][:])
                    if not np.array_equal(np.sort(saved_indices), np.arange(len(topology))) or any(
                            not np.allclose(setup[saved_indices], positions, atol=.002)
                            for path, setup, positions in setups):
                        raise ValueError("saved topology permutation mismatch")
                    with Dataset(Path(original["trajectory"]).with_name("checkpoint.chk")) as checkpoint:
                        if not np.allclose(data["positions"][0], checkpoint["positions"][0][:, saved_indices], atol=1e-6):
                            raise ValueError("saved/checkpoint coordinate order mismatch")
                    verified_trajectories.add(original["trajectory"])
                raw = np.asarray(raw, dtype=float) * 10
                box = np.asarray(data["box_vectors"][iteration, replica], dtype=float) * 10
                ligand_xyz = unwrap_ligand(raw[ligand_indices[endpoint]], ligand_bonds[endpoint], box)
                if not np.allclose(ligand_xyz, original["xyz_angstrom"], atol=1e-4):
                    raise ValueError("frozen endpoint ligand coordinates mismatch")
                receptor_xyz = whole_components(raw[indices], receptor_bonds, box, ligand_xyz.mean(axis=0), anchors)
            current_receptor = Protein(tuple(replace(atom, xyz=tuple(point))
                                             for atom, point in zip(receptor.atoms, receptor_xyz, strict=True)))
            current_ligand = replace(ligand, atoms=tuple(replace(atom, xyz=tuple(point))
                                           for atom, point in zip(ligand.atoms, ligand_xyz, strict=True)))
            deck, current_signature, current_fragments, current_metadata, current_atoms = deck_geometry(
                text, current_receptor, current_ligand, keys)
            if signature != current_signature or fragments != current_fragments or metadata["hop_pairs_gamess_indices"] != current_metadata["hop_pairs_gamess_indices"]:
                raise ValueError("same-frame geometry changed frozen fragmentation/charge/HOP")
            path = output / target / name / Path(original["trajectory"]).parent.name / f"iteration_{iteration}"
            path.mkdir(parents=True)
            (path / "input.inp").write_text(deck)
            geometry = geometry_contacts(current_atoms, len(ligand.atoms))
            write_json(path / "geometry.json", geometry)
            tasks.append({"target": target, "name": name, "kind": "same_frame_complex", "iteration": iteration,
                "replica": replica, "endpoint": endpoint, "trajectory": original["trajectory"],
                "directory": str(path), "input_sha256": digest(path / "input.inp"), "source": str(source),
                "source_sha256": digest(source), "runtime": original["runtime"], "parameters": original["parameters"],
                "geometry": geometry, "ligand_fragment": len(current_fragments)})
    if len(tasks) != 12 * len(iterations):
        raise ValueError("complete predefined same-frame schedule required")
    plan = {"tasks": tasks, "source_sha256": provenance, "mapping_audits": audits,
        "fixed_iterations": list(iterations),
        "driver_sha256": digest(driver_path or __file__), "labels_used": False, "main_pipeline_changed": False,
        "scope": "Same-frame receptor/physical endpoint ligand/original retained waters; regenerated caps, frozen chemistry and PCM. Not full explicit-solvent QM or binding free energy."}
    write_json(output / "frozen_protocol.json", plan)
    return plan


def execute(task, timeout):
    from fmotpiefinal.core.runner import run_gamess

    path = Path(task["directory"]) / "input.inp"
    if digest(path) != task["input_sha256"] or digest(task["source"]) != task["source_sha256"]:
        raise ValueError("frozen QM input changed")
    geometry = task["geometry"]
    if geometry["minimum_all_atom_contact_angstrom"] < .5 or geometry["minimum_heavy_contact_angstrom"] < 1.2:
        result = {"status": "FAILED_GEOMETRY", "error": "severe same-frame or cap contact; no atom removed"}
    else:
        result = run_gamess(path, path.parent / "run", gamess_dir=task["runtime"]["gamess_dir"],
            library_dir=task["runtime"]["library_dir"], ligand_fragment=task["ligand_fragment"],
            parameter_dir=task["parameters"], ncores=1, timeout_seconds=timeout)
    row = {**task, "result": result}
    write_json(path.parent / "summary.json", row)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("workspace", "run", "output"):
        parser.add_argument(f"--{name}", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout", type=int, default=1800)
    parser.add_argument("--freeze-only", action="store_true")
    parser.add_argument("--execute-frozen", action="store_true")
    args = parser.parse_args()
    if min(args.workers, args.timeout) < 1 or (args.freeze_only and args.execute_frozen):
        parser.error("positive budgets and mutually exclusive modes required")
    output = args.output.resolve()
    if args.execute_frozen:
        plan = json.loads((output / "frozen_protocol.json").read_text())
        if plan["driver_sha256"] != digest(__file__) or (output / "status.json").exists():
            raise ValueError("changed driver or already executed; preserve original run")
        for path, expected in plan["source_sha256"].items():
            if digest(path) != expected:
                raise ValueError(f"frozen source changed: {path}")
    else:
        plan = freeze(args)
    if args.freeze_only:
        print("FROZEN", len(plan["tasks"]), flush=True)
        return
    started = time.monotonic()
    rows = []
    write_json(output / "status.json", {"status": "RUNNING", "n_requested": len(plan["tasks"])})
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(execute, task, args.timeout): task for task in plan["tasks"]}
        for future in as_completed(futures):
            task = futures[future]
            try:
                row = future.result()
            except Exception as error:
                row = {**task, "result": {"status": "FAILED", "error": f"{type(error).__name__}: {error}"}}
            rows.append(row)
            write_json(output / "progress.json", rows)
            print(len(rows), len(plan["tasks"]), task["target"], task["name"], row["result"]["status"], flush=True)
    succeeded = sum(row["result"]["status"] == "SUCCEEDED" for row in rows)
    write_json(output / "status.json", {"status": "COMPLETE" if succeeded == len(rows) else "INCOMPLETE_WITH_FAILURES",
        "n_requested": len(plan["tasks"]), "n_succeeded": succeeded, "wall_seconds": time.monotonic() - started,
        "cohort_ranking_available": False, "rows": rows})


if __name__ == "__main__":
    main()
