"""Topology-aware gross H-C-H angle screen, not a force-field strain model."""
from __future__ import annotations

from itertools import combinations

import numpy as np


def audit_hch_angles(molecule, positions, *, minimum_degrees=80.0, maximum_degrees=150.0):
    positions = np.asarray(positions, dtype=float)
    if positions.shape != (molecule.GetNumAtoms(), 3) or not np.isfinite(positions).all():
        raise ValueError("matching finite explicit-H ligand positions required")
    if not 0 < minimum_degrees < maximum_degrees < 180:
        raise ValueError("finite ordered interior angle thresholds required")
    angles = []
    for center in molecule.GetAtoms():
        if center.GetAtomicNum() != 6 or str(center.GetHybridization()) != "SP3":
            continue
        hydrogens = [atom.GetIdx() for atom in center.GetNeighbors() if atom.GetAtomicNum() == 1]
        for first, last in combinations(hydrogens, 2):
            first_vector = positions[first]-positions[center.GetIdx()]
            last_vector = positions[last]-positions[center.GetIdx()]
            denominator = np.linalg.norm(first_vector)*np.linalg.norm(last_vector)
            if denominator == 0:
                raise ValueError("coincident bonded hydrogen and carbon")
            angle = float(np.degrees(np.arccos(np.clip(np.dot(first_vector, last_vector)/denominator, -1., 1.))))
            angles.append({"atoms_zero_based": [first, center.GetIdx(), last], "angle_degrees": angle,
                "gross_distortion_flag": angle < minimum_degrees or angle > maximum_degrees})
    return {"schema": "gross-sp3-hch-geometry-screen-v1", "angles": angles,
        "n_gross_distortions": sum(row["gross_distortion_flag"] for row in angles),
        "thresholds_degrees": [minimum_degrees, maximum_degrees],
        "scope": "only gross SP3 carbon H-C-H distortions; no chemical-state, energy, chirality or all-angle validation"}


def rebuild_flagged_carbon_hydrogens(molecule, positions):
    from rdkit import Chem
    positions = np.asarray(positions, dtype=float)
    audit = audit_hch_angles(molecule, positions)
    centers = {row["atoms_zero_based"][1] for row in audit["angles"] if row["gross_distortion_flag"]}
    if not centers:
        return positions.copy(), []
    if any(atom.GetAtomicNum() == 1 and atom.GetIsotope() for atom in molecule.GetAtoms()):
        raise ValueError("isotopic H reconstruction requires an explicit validated mapping")
    working = Chem.Mol(molecule)
    for atom in working.GetAtoms():
        atom.SetIntProp("source_atom_index", atom.GetIdx())
    conformer = working.GetConformer()
    for index, point in enumerate(positions):
        conformer.SetAtomPosition(index, tuple(map(float, point)))
    rebuilt = Chem.AddHs(Chem.RemoveHs(working), addCoords=True)
    output = positions.copy()
    changed = []
    for center in centers:
        originals = sorted(atom.GetIdx() for atom in working.GetAtomWithIdx(center).GetNeighbors() if atom.GetAtomicNum() == 1)
        matches = [atom for atom in rebuilt.GetAtoms() if atom.GetAtomicNum() != 1 and atom.GetIntProp("source_atom_index") == center]
        if len(matches) != 1:
            raise ValueError("heavy atom identity changed during H reconstruction")
        generated = sorted(atom.GetIdx() for atom in matches[0].GetNeighbors() if atom.GetAtomicNum() == 1)
        if len(generated) != len(originals):
            raise ValueError("hydrogen count changed during reconstruction")
        for original_index, generated_index in zip(originals, generated, strict=True):
            output[original_index] = rebuilt.GetConformer().GetAtomPosition(generated_index)
            changed.append(original_index)
    if audit_hch_angles(molecule, output)["n_gross_distortions"]:
        raise ValueError("reconstruction did not resolve gross H-C-H geometry")
    return output, sorted(changed)
