"""Detect ligand motion associated with original sidechain overlap, without scores."""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import math

import numpy as np
from scipy.spatial import cKDTree

from .sidechain_policy import BACKBONE, RADII, canonical_name, key, residue_id
from .receptor_chemical_graph import CANONICAL


def add_displaced_clash_seeds(protein, molecules, bound_positions, seed_plan,
                              minimum_displacement_angstrom=0.3):
    if not math.isfinite(minimum_displacement_angstrom) or minimum_displacement_angstrom <= 0:
        raise ValueError("motion threshold must be finite and positive")
    if set(molecules) != set(bound_positions) or not molecules:
        raise ValueError("matching complete ligand identities required")
    result = deepcopy(seed_plan)
    nodes = {tuple(node["residue"]): node for node in result["nodes"]}
    atoms = [atom for atom in protein.atoms if atom.record == "ATOM" and atom.atom_name not in BACKBONE
        and atom.element.upper() != "H" and canonical_name(atom.residue_name) in CANONICAL - {"ALA", "GLY", "PRO"}]
    coordinates = np.asarray([atom.xyz for atom in atoms], dtype=float).reshape((-1, 3))
    tree = cKDTree(coordinates)
    selected = {tuple(residue) for residue in seed_plan["selected_residues"]}
    events = []
    for name in sorted(molecules):
        molecule = molecules[name]
        original = np.asarray(molecule.GetConformer().GetPositions(), dtype=float)
        bound = np.asarray(bound_positions[name], dtype=float)
        if original.shape != bound.shape or bound.shape != (molecule.GetNumAtoms(), 3):
            raise ValueError("ligand coordinates do not match molecular atom inventory")
        if not np.isfinite(bound).all() or not np.isfinite(original).all():
            raise ValueError("nonfinite ligand coordinates")
        for ligand_atom in molecule.GetAtoms():
            ligand_index = ligand_atom.GetIdx()
            element = ligand_atom.GetSymbol().upper()
            if ligand_atom.GetAtomicNum() == 1 or element not in RADII:
                continue
            displacement = float(np.linalg.norm(bound[ligand_index]-original[ligand_index]))
            if displacement < minimum_displacement_angstrom:
                continue
            for receptor_index in tree.query_ball_point(original[ligand_index], seed_plan["policy"]["contact_radius_angstrom"]):
                atom = atoms[receptor_index]
                if atom.element.upper() not in RADII:
                    continue
                distance = float(np.linalg.norm(original[ligand_index]-coordinates[receptor_index]))
                overlap = RADII[atom.element.upper()] + RADII[element] - distance
                polar = atom.element.upper() in {"N", "O", "S"} and element in {"N", "O", "S"}
                if overlap < seed_plan["policy"]["heavy_overlap_angstrom"] or (
                        polar and distance >= seed_plan["policy"]["polar_pair_min_distance_angstrom"]):
                    continue
                event = {"code": "ORIGINAL_HEAVY_OVERLAP_WITH_LIGAND_DISPLACEMENT", "ligand": name,
                    "receptor": residue_id(key(atom)), "receptor_atom": atom.atom_name,
                    "ligand_atom_zero_based": ligand_index, "original_distance_angstrom": distance,
                    "bound_distance_angstrom": float(np.linalg.norm(bound[ligand_index]-coordinates[receptor_index])),
                    "overlap_angstrom": overlap, "ligand_displacement_angstrom": displacement,
                    "interpretation": "co-occurring overlap and displacement; not causal proof or energy strain"}
                node = nodes[key(atom)]
                node["evidence"].append(event)
                events.append(event)
                if node["protection"]:
                    node["action"] = "REVIEW_PROTECTED_NETWORK_KEEP_FIXED"
                else:
                    selected.add(key(atom))
                    node["action"] = "RESTRAINED_SIDECHAIN_RELAX"
    result["selected_residues"] = sorted(map(list, selected))
    result["schema"] = "structure-driven-sidechain-motion-seeds-v1"
    result["motion_seed_rule"] = {"minimum_heavy_displacement_angstrom": minimum_displacement_angstrom,
        "requires_original_heavy_overlap": True, "scores_or_labels_used": False}
    result["motion_events"] = events
    result["parent_plan_sha256"] = seed_plan["plan_sha256"]
    result["action"] = "RESTRAINED_SIDECHAIN_RELAX" if selected else "KEEP_BASELINE_UNCHANGED"
    result.pop("plan_sha256")
    result["plan_sha256"] = hashlib.sha256(json.dumps(result, sort_keys=True, allow_nan=False).encode()).hexdigest()
    return result
