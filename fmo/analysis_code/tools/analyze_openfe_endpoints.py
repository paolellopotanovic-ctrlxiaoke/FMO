"""Diagnose saved physical endpoints without selecting frames by affinity."""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np
import yaml

from run_matched_openfe import digest, write_json


def minimum_image(displacement, box):
    fractional = np.asarray(displacement) @ np.linalg.inv(box)
    return (fractional - np.rint(fractional)) @ box


def distances(first, second, box):
    return np.linalg.norm(minimum_image(first[:, None] - second[None, :], box), axis=-1)


def fit_coordinates(moving, reference, query):
    centered = moving - moving.mean(axis=0)
    target = reference - reference.mean(axis=0)
    left, _, right = np.linalg.svd(centered.T @ target)
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(left @ right)
    rotation = left @ correction @ right
    return (query - moving.mean(axis=0)) @ rotation + reference.mean(axis=0)


def fitted_rmsd(moving, reference):
    fitted = fit_coordinates(moving, reference, moving)
    return float(np.sqrt(np.mean(np.sum((fitted - reference) ** 2, axis=1))))


def torsion_degrees(coordinates, box):
    first, middle, last = minimum_image(np.diff(coordinates, axis=0), box)
    axis = middle / np.linalg.norm(middle)
    first = -first
    normal_first = first - np.dot(first, axis) * axis
    normal_last = last - np.dot(last, axis) * axis
    return float(np.degrees(np.arctan2(np.dot(np.cross(axis, normal_first), normal_last), np.dot(normal_first, normal_last))))


def circular_summary(angles):
    radians = np.radians(angles)
    resultant = np.mean(np.exp(1j * radians))
    return {"circular_mean_degrees": float(np.degrees(np.angle(resultant))),
            "resultant_length": float(abs(resultant)),
            "circular_std_degrees": float(np.degrees(np.sqrt(-2 * np.log(max(min(abs(resultant), 1), 1e-15)))))}


def summarize(values):
    values = np.asarray(values, dtype=float)
    return {"mean": float(values.mean()), "std": float(values.std()),
            "min": float(values.min()), "max": float(values.max())}


def read_topology(path):
    rows = [line for line in path.read_text().splitlines()
            if line.startswith(("ATOM  ", "HETATM"))]
    return [{"name": line[12:16].strip(), "resname": line[17:20].strip(),
             "residue": line[21:27].strip(), "element": line[76:78].strip(),
             "class": float(line[60:66]),
             "xyz": [float(line[30:38]), float(line[38:46]), float(line[46:54])]}
            for line in rows]


def endpoint_indices(atoms, endpoint):
    allowed = (0.25, 0.5) if endpoint == 0 else (0.5, 0.75)
    return np.array([index for index, atom in enumerate(atoms)
                     if atom["class"] in allowed and atom["element"] != "H"], dtype=int)


def ligand_bonds(path, ligand):
    lines = path.read_text().splitlines()
    serials = {int(line[6:11]): index for index, line in enumerate(
        line for line in lines if line.startswith(("ATOM  ", "HETATM")))}
    local = {int(index): offset for offset, index in enumerate(ligand)}
    bonds = set()
    for line in lines:
        if line.startswith("CONECT"):
            numbers = [int(line[start:start + 5]) for start in range(6, len(line), 5) if line[start:start + 5].strip()]
            for serial in numbers[1:]:
                first, second = serials[numbers[0]], serials[serial]
                if first in local and second in local:
                    bonds.add(tuple(sorted((local[first], local[second]))))
    return sorted(bonds)


def unwrap_ligand(coordinates, bonds, box):
    adjacency = defaultdict(list)
    for first, second in bonds:
        adjacency[first].append(second)
        adjacency[second].append(first)
    whole = np.array(coordinates, copy=True)
    visited = {0}
    pending = [0]
    while pending:
        first = pending.pop()
        for second in adjacency[first]:
            if second not in visited:
                whole[second] = whole[first] + minimum_image(coordinates[second] - coordinates[first], box)
                visited.add(second)
                pending.append(second)
    if len(visited) != len(coordinates):
        raise ValueError("endpoint heavy-atom bond graph is disconnected")
    return whole


def analyze_leg(root, first_iteration=400, final_iteration=2000, stride=40):
    from netCDF4 import Dataset, chartostring

    from gufe import LigandAtomMapping

    topologies = sorted((root / "shared").glob("*/hybrid_system.pdb"))
    atoms = read_topology(topologies[0])
    signature = [(atom["name"], atom["element"], atom["class"]) for atom in atoms]
    for path in topologies[1:]:
        other = read_topology(path)
        if [(atom["name"], atom["element"], atom["class"]) for atom in other] != signature:
            raise ValueError("repeat atom-index signatures differ; explicit setup mapping required")
    reference = np.array([atom["xyz"] for atom in atoms])
    waters = np.array([index for index, atom in enumerate(atoms)
                       if atom["resname"] in ("HOH", "WAT") and atom["element"] == "O"], dtype=int)
    protein = np.array([index for index, atom in enumerate(atoms)
                        if atom["class"] == 0 and atom["element"] not in ("H", "")
                        and atom["resname"] not in ("HOH", "WAT", "NA", "CL", "TLA", "UNK")], dtype=int)
    residue_groups = defaultdict(list)
    for index in protein:
        residue_groups[f'{atoms[index]["resname"]}:{atoms[index]["residue"]}'].append(int(index))
    physical = {endpoint: endpoint_indices(atoms, endpoint) for endpoint in (0, 10)}
    mapping = LigandAtomMapping.from_json(root / "mapping.json")
    for endpoint, component in ((0, mapping.componentA), (10, mapping.componentB)):
        if len(physical[endpoint]) != component.to_rdkit().GetNumHeavyAtoms():
            raise ValueError("endpoint ligand heavy-atom count differs from registered component")
    bonds = {endpoint: ligand_bonds(topologies[0], ligand) for endpoint, ligand in physical.items()}
    union_ligand = np.unique(np.concatenate(list(physical.values())))
    local_residues = {key: indices for key, indices in residue_groups.items()
                      if np.linalg.norm(reference[indices, None] - reference[union_ligand][None], axis=-1).min() < 6.0}
    chi1 = {}
    for key, indices in local_residues.items():
        named = {atoms[index]["name"]: index for index in indices}
        terminal = next((name for name in ("CG", "CG1", "OG", "OG1", "SG") if name in named), None)
        if terminal and all(name in named for name in ("N", "CA", "CB")):
            chi1[key] = [named[name] for name in ("N", "CA", "CB", terminal)]
    alpha_carbons = np.array([index for indices in local_residues.values() for index in indices if atoms[index]["name"] == "CA"], dtype=int)
    results = []
    skipped = []
    for history in sorted((root / "shared").glob("*/simulation_real_time_analysis.yaml")):
        timeline = yaml.safe_load(history.read_text())
        if not timeline or timeline[-1]["percent_complete"] != 100:
            skipped.append(str(history.parent))
            continue
        trajectory = history.parent / "simulation.nc"
        rows = {endpoint: [] for endpoint in physical}
        with Dataset(trajectory) as data:
            for endpoint, expected in ((0, "0.0"), (10, "1.0")):
                serialized = chartostring(data[f"thermodynamic_states/state{endpoint}"][:]).item()
                state = yaml.load(serialized, Loader=yaml.BaseLoader)
                if state["composable_states"][0]["parameters"]["lambda_torsions"] != expected:
                    raise ValueError("saved endpoint lambda does not match requested physical state")
            if data["positions"].shape[2] != len(atoms) or data["positions"].units != "nm":
                raise ValueError("full-system coordinates in nm required")
            if int(data.PositionInterval) != stride or final_iteration > int(data["last_iteration"][0]):
                raise ValueError("requested completed fixed sampling schedule unavailable")
            saved_indices = np.asarray(data["analysis_particle_indices"][:])
            if not np.array_equal(np.sort(saved_indices), np.arange(len(atoms))):
                raise ValueError("saved full-system atom indices must be a permutation")
            with Dataset(trajectory.with_name("checkpoint.chk")) as checkpoint:
                if not np.allclose(data["positions"][0], checkpoint["positions"][0][:, saved_indices], atol=1e-6):
                    raise ValueError("analysis/checkpoint atom-index mapping mismatch")
            setup_positions = np.load(topologies[0].with_name("hybrid_positions.npy")) * 10.0
            if not np.allclose(setup_positions[saved_indices], reference, atol=0.002):
                raise ValueError("PDB topology and saved selection order mismatch")
            for iteration in range(first_iteration, final_iteration + 1, stride):
                states = np.asarray(data["states"][iteration])
                for endpoint, ligand in physical.items():
                    replicas = np.flatnonzero(states == endpoint)
                    if len(replicas) != 1:
                        raise ValueError("endpoint must have exactly one physical replica")
                    replica = int(replicas[0])
                    raw = data["positions"][iteration, replica]
                    if np.ma.getmaskarray(raw).any():
                        raise ValueError("requested frame was not saved")
                    coordinates = np.asarray(raw) * 10.0
                    box = np.asarray(data["box_vectors"][iteration, replica]) * 10.0
                    ligand_xyz = coordinates[ligand]
                    ligand_xyz = unwrap_ligand(ligand_xyz, bonds[endpoint], box)
                    hetero = np.array([index for index in ligand if atoms[index]["element"] in ("N", "O", "S")], dtype=int)
                    water_distances = distances(coordinates[waters], coordinates[hetero], box).min(axis=1)
                    near_waters = waters[water_distances < 3.5]
                    radius = np.sqrt(np.mean(np.sum((ligand_xyz - ligand_xyz.mean(axis=0)) ** 2, axis=1)))
                    row = {"iteration": iteration, "production_ns": iteration * 0.0025,
                           "replica": replica, "ligand_self_fit_rmsd_angstrom": fitted_rmsd(ligand_xyz, reference[ligand]),
                           "heavy_atom_unweighted_radius_of_gyration_angstrom": float(radius),
                           "water_oxygen_within_3p5A_of_ligand_heteroatom": int(len(near_waters)),
                           "residue_contacts": {}, "sidechain_chi1_degrees": {}}
                    if len(protein):
                        if len(alpha_carbons) < 3:
                            raise ValueError("at least three local CA anchors required")
                        expected_offsets = reference[alpha_carbons] - reference[ligand[0]]
                        backbone_xyz = ligand_xyz[0] + expected_offsets + minimum_image(coordinates[alpha_carbons] - ligand_xyz[0] - expected_offsets, box)
                        fitted = fit_coordinates(backbone_xyz, reference[alpha_carbons], ligand_xyz)
                        row["local_receptor_fit_ligand_pose_rmsd_angstrom"] = float(np.sqrt(np.mean(np.sum((fitted - reference[ligand]) ** 2, axis=1))))
                        row["minimum_protein_ligand_heavy_distance_angstrom"] = float(distances(coordinates[protein], coordinates[ligand], box).min())
                        protein_hetero = np.array([index for index in protein if atoms[index]["element"] in ("N", "O", "S")], dtype=int)
                        row["distance_only_water_bridge_count"] = int(np.sum(distances(coordinates[near_waters], coordinates[protein_hetero], box).min(axis=1) < 3.5)) if len(near_waters) else 0
                        row["residue_contacts"] = {key: float(distances(coordinates[indices], coordinates[ligand], box).min()) for key, indices in local_residues.items()}
                        row["sidechain_chi1_degrees"] = {key: torsion_degrees(coordinates[indices], box) for key, indices in chi1.items()}
                    rows[endpoint].append(row)
        for endpoint, frames in rows.items():
            numeric_keys = [key for key, value in frames[0].items() if isinstance(value, (float, int)) and key not in ("iteration", "replica", "production_ns")]
            contacts = {key: {"minimum_distance_angstrom": summarize([frame["residue_contacts"][key] for frame in frames]),
                              "fraction_below_3p5A": float(np.mean([frame["residue_contacts"][key] < 3.5 for frame in frames]))}
                        for key in local_residues}
            results.append({"trajectory": str(trajectory), "history_sha256": digest(history),
                            "endpoint": endpoint, "frame_count": len(frames),
                            "heavy_atom_count": len(physical[endpoint]),
                            "summary": {key: summarize([frame[key] for frame in frames]) for key in numeric_keys},
                            "sidechain_chi1": {key: circular_summary([frame["sidechain_chi1_degrees"][key] for frame in frames]) for key in chi1},
                            "residue_contacts": contacts, "frames": frames})
    return {"topology_sources": {str(path): digest(path) for path in topologies},
            "repeat_atom_signatures_verified": True, "complete_repeats": results,
            "incomplete_repeats_skipped": skipped}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads((args.run / "frozen_protocol.json").read_text())
    report = {"schema": "openfe-endpoint-structure-diagnostics-v1", "targets": {},
              "analysis_driver_sha256": digest(__file__),
              "frozen_protocol_sha256": digest(args.run / "frozen_protocol.json"),
              "selection": "All completed repeats at snapshot time; production 1-5 ns inclusive every 0.1 ns; actual lambda endpoint, not fixed replica.",
              "affinity_used_for_frame_selection": False, "main_pipeline_changed": False,
              "limitations": ["Distance contacts are not hydrogen bonds or energetic contributions.",
                              "Frames are correlated; summaries are descriptive, not independent-sample uncertainty.",
                              "Hybrid dummy atoms are excluded from the selected endpoint ligand.",
                              "Self-fit RMSD measures internal shape; pose RMSD separately uses local receptor CA anchors.",
                              "No entropy, binding free energy or full-cohort ranking inferred from these summaries."]}
    for target in ("BACE", "CDK8"):
        root = args.run / "production" / target / "edge_000"
        first, second = plan["targets"][target]["edges"][0]
        report["targets"][target] = {"endpoint_names": {"0": first, "10": second},
            **{leg: analyze_leg(root / leg) for leg in ("solution", "bound")}}
    write_json(args.output, report)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
