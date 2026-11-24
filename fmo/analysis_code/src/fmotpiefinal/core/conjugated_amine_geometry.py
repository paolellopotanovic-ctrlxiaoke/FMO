"""Label-free diagnostics for preparation-induced guanidino NH2 distortion."""
from __future__ import annotations

from itertools import permutations

import numpy as np


def conjugated_nh2_groups(molecule):
    groups = []
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.GetFormalCharge() != 0:
            continue
        hydrogens = sorted(neighbor.GetIdx() for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 1)
        heavy = [neighbor for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() > 1]
        if len(hydrogens) != 2 or len(heavy) != 1 or heavy[0].GetAtomicNum() != 6:
            continue
        center = heavy[0]
        nitrogens = sorted(neighbor.GetIdx() for neighbor in center.GetNeighbors() if neighbor.GetAtomicNum() == 7)
        if center.GetDegree() != 3 or len(nitrogens) != 3 or not any(bond.GetBondTypeAsDouble() == 2. for bond in center.GetBonds()):
            continue
        groups.append({"nitrogen": atom.GetIdx(), "center": center.GetIdx(), "hydrogens": hydrogens,
                       "other_nitrogens": [index for index in nitrogens if index != atom.GetIdx()]})
    return groups


def group_departure(positions, group):
    points = np.asarray(positions, dtype=float)
    center = points[group["center"]]
    neighbors = group["other_nitrogens"]
    normal = np.cross(points[neighbors[0]]-center, points[neighbors[1]]-center)
    if not np.isfinite(normal).all() or np.linalg.norm(normal) < 1.e-8:
        raise ValueError("undefined conjugated heavy-atom plane")
    normal /= np.linalg.norm(normal)
    departures = []
    for index in group["hydrogens"]:
        vector = points[index]-points[group["nitrogen"]]
        if not np.isfinite(vector).all() or np.linalg.norm(vector) < .6:
            raise ValueError("invalid conjugated N-H geometry")
        departures.append(float(np.degrees(np.arcsin(np.clip(abs(vector@normal)/np.linalg.norm(vector), 0., 1.)))))
    return max(departures)


def induced_conjugated_nh2(molecule, original, bound, *, departure_degrees=30., increase_degrees=20.):
    original, bound = np.asarray(original, float), np.asarray(bound, float)
    if (original.shape != (molecule.GetNumAtoms(), 3) or bound.shape != original.shape
            or not np.isfinite(original).all() or not np.isfinite(bound).all()
            or departure_degrees <= 0 or increase_degrees <= 0):
        raise ValueError("finite matched geometry and positive screening thresholds required")
    rows = []
    for group in conjugated_nh2_groups(molecule):
        before, after = group_departure(original, group), group_departure(bound, group)
        rows.append({**group, "original_departure_degrees": before, "bound_departure_degrees": after,
            "new_deformation_flag": before < departure_degrees and after > departure_degrees and after-before >= increase_degrees,
            "already_flagged_input": before >= departure_degrees})
    return rows


def reconstruct_flagged_nh2(molecule, positions, diagnostics):
    points = np.asarray(positions, float).copy()
    if points.shape != (molecule.GetNumAtoms(), 3) or not np.isfinite(points).all():
        raise ValueError("finite matched ligand geometry required")
    groups = {group["nitrogen"]: group for group in conjugated_nh2_groups(molecule)}
    for row in diagnostics:
        if not row["new_deformation_flag"]:
            continue
        group = groups.get(row["nitrogen"])
        if group is None or any(group[key] != row[key] for key in group):
            raise ValueError("diagnosed chemical group no longer matches molecule")
        nitrogen, center = group["nitrogen"], group["center"]
        axis = points[center]-points[nitrogen]
        axis /= np.linalg.norm(axis)
        other = points[group["other_nitrogens"][0]]-points[center]
        normal = np.cross(axis, other)
        if np.linalg.norm(normal) < 1.e-8:
            raise ValueError("undefined repair plane")
        normal /= np.linalg.norm(normal)
        perpendicular = np.cross(normal, axis)
        directions = [-.5*axis+np.sqrt(.75)*perpendicular, -.5*axis-np.sqrt(.75)*perpendicular]
        indices = group["hydrogens"]
        lengths = [float(np.linalg.norm(points[index]-points[nitrogen])) for index in indices]
        if min(lengths) < .6 or max(lengths) > 1.45:
            raise ValueError("invalid original N-H bond length")
        alternatives = [[points[nitrogen]+lengths[index]*directions[assignment[index]] for index in range(2)]
                        for assignment in permutations(range(2))]
        chosen = min(alternatives, key=lambda candidate: sum(float(np.sum((candidate[index]-points[atom])**2)) for index, atom in enumerate(indices)))
        for index, point in zip(indices, chosen, strict=True):
            points[index] = point
    return points
