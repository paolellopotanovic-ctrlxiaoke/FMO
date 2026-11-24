"""Label-independent acceptance evidence for a narrow preparation-artifact repair."""
from __future__ import annotations

from dataclasses import asdict

import numpy as np
from scipy.spatial import cKDTree

from .sidechain_policy import RADII, SidechainPolicy, amide_changes


def severe_overlaps(protein, molecule, positions, policy=None):
    policy = policy or SidechainPolicy()
    positions = np.asarray(positions, dtype=float)
    if positions.shape != (molecule.GetNumAtoms(), 3) or not np.isfinite(positions).all():
        raise ValueError("finite matching ligand atom coordinates required")
    atoms = [atom for atom in protein.atoms if atom.element.upper() in RADII]
    points = np.asarray([atom.xyz for atom in atoms], dtype=float).reshape((-1, 3))
    if not np.isfinite(points).all():
        raise ValueError("nonfinite receptor coordinates")
    tree = cKDTree(points)
    maximum_radius = max(RADII.values())
    overlaps = []
    for ligand_atom in molecule.GetAtoms():
        element = ligand_atom.GetSymbol().upper()
        if element not in RADII:
            continue
        index = ligand_atom.GetIdx()
        for receptor_index in tree.query_ball_point(positions[index], maximum_radius+RADII[element]):
            atom = atoms[receptor_index]
            distance = float(np.linalg.norm(points[receptor_index]-positions[index]))
            overlap = RADII[atom.element.upper()]+RADII[element]-distance
            polar = element in {"N", "O", "S"} and atom.element.upper() in {"N", "O", "S"}
            if overlap >= policy.heavy_overlap_angstrom and (not polar or distance < policy.polar_pair_min_distance_angstrom):
                overlaps.append({"receptor_atom_index": atom.index, "receptor_atom": atom.atom_name,
                    "residue": [atom.chain, atom.residue_number, atom.residue_name],
                    "ligand_atom_zero_based": index, "distance_A": distance, "overlap_A": overlap})
    return overlaps


def diagnose_case(molecule, original, baseline, candidate, baseline_receptor, candidate_receptor, policy=None):
    policy = policy or SidechainPolicy()
    baseline_amides = amide_changes(molecule, np.asarray(original), np.asarray(baseline), policy)
    candidate_amides = amide_changes(molecule, np.asarray(original), np.asarray(candidate), policy)
    if len(baseline_amides) != len(candidate_amides):
        raise ValueError("amide atom inventory changed")
    repairs = []
    new_flags = []
    for before, after in zip(baseline_amides, candidate_amides, strict=True):
        if before["atoms_zero_based"] != after["atoms_zero_based"]:
            raise ValueError("amide atom identities changed")
        if before["new_deformation_flag"]:
            improvement = before["bound_departure_deg"]-after["bound_departure_deg"]
            repairs.append({"atoms_zero_based": before["atoms_zero_based"],
                "baseline_departure_deg": before["bound_departure_deg"],
                "candidate_departure_deg": after["bound_departure_deg"],
                "improvement_deg": improvement,
                "repaired": improvement >= policy.amide_departure_increase_degrees and
                    after["bound_departure_deg"] < policy.amide_departure_degrees})
        elif after["new_deformation_flag"]:
            new_flags.append(after)
    baseline_overlaps = severe_overlaps(baseline_receptor, molecule, baseline, policy)
    candidate_overlaps = severe_overlaps(candidate_receptor, molecule, candidate, policy)
    return {"baseline_new_amide_flags": repairs, "new_candidate_amide_flags": new_flags,
        "baseline_severe_overlaps": baseline_overlaps, "candidate_severe_overlaps": candidate_overlaps,
        "severe_overlap_count_not_increased": len(candidate_overlaps) <= len(baseline_overlaps),
        "baseline_amides": baseline_amides, "candidate_amides": candidate_amides}


def decide_target(cases, expected_names, *, geometry_verified, qm_verified, policy=None):
    policy = policy or SidechainPolicy()
    if not expected_names or len(expected_names) != len(set(expected_names)) or set(cases) != set(expected_names):
        raise ValueError("complete unique case identities required")
    flags = [flag for case in cases.values() for flag in case["baseline_new_amide_flags"]]
    checks = {"baseline_preparation_artifact_detected": bool(flags),
        "all_triggered_amides_repaired": bool(flags) and all(flag["repaired"] for flag in flags),
        "no_new_amide_flags": not any(case["new_candidate_amide_flags"] for case in cases.values()),
        "no_case_increases_severe_overlaps": all(case["severe_overlap_count_not_increased"] for case in cases.values()),
        "complete_geometry_verified": geometry_verified is True,
        "complete_qm_verified": qm_verified is True}
    if not flags:
        action = "KEEP_BASELINE_MECHANISM_NOT_DETECTED"
    elif all(checks.values()):
        action = "ACCEPT_LOCAL_REPAIR_RESEARCH_CONTROL"
    else:
        action = "KEEP_BASELINE_REPAIR_REQUIRES_REVIEW"
    return {"action": action, "checks": checks, "policy": asdict(policy),
        "affinity_labels_used": False, "scores_used": False,
        "scope": "new k10-induced amide deformation only, not a general affinity-quality classifier",
        "production_default_changed": False}
