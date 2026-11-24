"""Measure sampled water contacts omitted by a fixed QM water inventory."""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np

from analyze_openfe_endpoints import distances, read_topology
from analyze_cdk8_geometry_evidence import fragments, xyz
from run_endpoint_complex_tpie import WATERS
from run_matched_openfe import digest, write_json


def coverage_counts(nearby, represented):
    nearby, represented = set(nearby), set(represented)
    return {"n_nearby_water_oxygens": len(nearby),
            "n_nearby_represented_water_oxygens": len(nearby & represented),
            "n_nearby_omitted_water_oxygens": len(nearby - represented),
            "nearby_omitted_saved_atom_indices": sorted(nearby - represented)}


def main():
    from netCDF4 import Dataset
    from fmotpiefinal.core.structure import read_pdb, residue_key

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    args = parser.parse_args()
    root = args.run.resolve()
    plan_path = root / "frozen_protocol.json"
    plan = json.loads(plan_path.read_text())
    registry_paths = [path for path in plan["source_sha256"] if path.endswith("graph_sidechain_transfer_20261021/frozen_protocol.json")]
    if len(registry_paths) != 1:
        raise ValueError("unique frozen cohort registry required")
    registry = Path(registry_paths[0])
    if digest(registry) != plan["source_sha256"][str(registry)]:
        raise ValueError("registry changed")
    cohorts = json.loads(registry.read_text())["cohorts"]
    report = {"status": "RUNNING", "labels_used": False, "main_pipeline_changed": False,
        "protocol_sha256": digest(plan_path), "driver_sha256": digest(__file__),
        "cutoff_angstrom": 3.5, "selection": "All frozen frames, water oxygen to ligand N/O/S distance; no frame filtering.",
        "limitations": ["Distance contacts are not hydrogen bonds or energetic importance.",
            "Omitted means absent as explicit water in this fixed QM pocket, not absent from PCM polarization.",
            "Counts do not establish that inserting or deleting waters improves affinity ranking."], "targets": {}}
    for target, audit in plan["mapping_audits"].items():
        receptor = read_pdb(cohorts[target]["protein"])
        if digest(cohorts[target]["protein"]) != plan["source_sha256"][cohorts[target]["protein"]]:
            raise ValueError("source receptor changed")
        keys = set(audit["fixed_pocket_keys"])
        represented = [saved for atom, saved in zip(receptor.atoms, audit["source_to_saved_indices"], strict=True)
                       if atom.residue_name in WATERS and atom.element.upper() == "O" and residue_key(atom) in keys]
        grouped = defaultdict(list)
        for task in plan["tasks"]:
            if task["target"] == target:
                grouped[task["trajectory"]].append(task)
        rows = []
        for trajectory, tasks in sorted(grouped.items()):
            setup_paths = sorted((Path(trajectory).parents[1]).glob("*/hybrid_system.pdb"))
            if not setup_paths:
                raise ValueError("hybrid topology missing")
            topology = read_topology(setup_paths[0])
            waters = np.array([index for index, atom in enumerate(topology)
                               if atom["class"] == 0 and atom["resname"] in WATERS and atom["element"].upper() == "O"], int)
            if not set(represented) <= set(waters.tolist()):
                raise ValueError("represented water identity differs from topology")
            with Dataset(trajectory) as data:
                for task in sorted(tasks, key=lambda item: (item["iteration"], item["endpoint"])):
                    path = Path(task["directory"]) / "input.inp"
                    if digest(path) != task["input_sha256"]:
                        raise ValueError("input changed")
                    text = path.read_text()
                    symbols, positions = xyz(text)
                    names, membership = fragments(text)
                    ligand_hetero = [index - 1 for index in membership[-1] if symbols[index - 1] in {"N", "O", "S"}]
                    if not ligand_hetero:
                        raise ValueError("water polar-contact audit requires ligand N/O/S atoms")
                    iteration, replica = task["iteration"], task["replica"]
                    if int(data["states"][iteration, replica]) != task["endpoint"]:
                        raise ValueError("physical endpoint changed")
                    raw = data["positions"][iteration, replica]
                    if np.ma.getmaskarray(raw).any():
                        raise ValueError("missing frame")
                    coordinates = np.asarray(raw, float) * 10
                    box = np.asarray(data["box_vectors"][iteration, replica], float) * 10
                    nearby = waters[distances(coordinates[waters], positions[ligand_hetero], box).min(axis=1) < 3.5]
                    rows.append({"name": task["name"], "trajectory": trajectory, "iteration": iteration,
                                 "endpoint": task["endpoint"], **coverage_counts(nearby.tolist(), represented)})
        by_ligand = {}
        for name in sorted({row["name"] for row in rows}):
            repeats = []
            for trajectory in grouped:
                frames = [row for row in rows if row["name"] == name and row["trajectory"] == trajectory]
                if sorted(row["iteration"] for row in frames) != plan["fixed_iterations"]:
                    raise ValueError("complete fixed water audit schedule required")
                repeats.append({"trajectory": trajectory, "n_frames": len(frames),
                    "fraction_frames_with_omitted_nearby_water": float(np.mean([row["n_nearby_omitted_water_oxygens"] > 0 for row in frames])),
                    "mean_nearby_count": float(np.mean([row["n_nearby_water_oxygens"] for row in frames])),
                    "mean_omitted_count": float(np.mean([row["n_nearby_omitted_water_oxygens"] for row in frames]))})
            by_ligand[name] = {"repeats": repeats}
        report["targets"][target] = {"n_explicit_qm_waters": len(represented), "ligands": by_ligand, "frames": rows}
        write_json(root / "water_coverage.json", report)
        print(target, json.dumps(by_ligand), flush=True)
    report["status"] = "COMPLETE"
    write_json(root / "water_coverage.json", report)


if __name__ == "__main__":
    main()
