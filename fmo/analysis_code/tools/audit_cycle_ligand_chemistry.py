"""Label-free bond, stereochemistry and H-parent checks on cycle geometries."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from rdkit import Chem

from fmotpiefinal.core.bound_geometry import load_bound_ligand
from run_cdk8_solution_cycle import digest, write_json


def stereochemistry(molecule, positions):
    molecule = Chem.Mol(molecule)
    positions = np.asarray(positions, dtype=float)
    if positions.shape != (molecule.GetNumAtoms(), 3) or not np.isfinite(positions).all():
        raise ValueError("invalid geometry or changed atom count")
    for index, point in enumerate(positions):
        molecule.GetConformer().SetAtomPosition(index, point)
    Chem.AssignStereochemistryFrom3D(molecule, replaceExistingTags=True)
    Chem.AssignStereochemistry(molecule, force=True, cleanIt=True)
    return {str(atom.GetIdx()): atom.GetProp("_CIPCode") for atom in molecule.GetAtoms() if atom.HasProp("_CIPCode")}


def geometry_screen(molecule, positions, reference_stereo):
    positions = np.asarray(positions, dtype=float)
    stereo = stereochemistry(molecule, positions)
    stretched = []
    for bond in molecule.GetBonds():
        first, second = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        distance = float(np.linalg.norm(positions[first]-positions[second]))
        hydrogen = any(molecule.GetAtomWithIdx(index).GetAtomicNum() == 1 for index in (first, second))
        if distance < .6 or distance > (1.45 if hydrogen else 2.3):
            stretched.append({"atoms": [first, second], "distance_A": distance})
    transfers = []
    heavy = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1]
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() != 1 or atom.GetDegree() != 1:
            continue
        index, parent = atom.GetIdx(), atom.GetNeighbors()[0].GetIdx()
        distances = np.linalg.norm(positions[heavy]-positions[index], axis=1)
        nearest = heavy[int(np.argmin(distances))]
        old_distance = float(np.linalg.norm(positions[parent]-positions[index]))
        if nearest != parent and old_distance > 1.25 and float(distances.min()) < 1.2:
            transfers.append({"hydrogen": index, "old_parent": parent, "new_nearest_heavy": nearest,
                              "old_bond_A": old_distance, "new_distance_A": float(distances.min())})
    return {"stereochemistry_3d": stereo, "stereo_changed": stereo != reference_stereo,
            "bond_distance_flags": stretched, "possible_hydrogen_transfer_flags": transfers,
            "passed_screen": stereo == reference_stereo and not stretched and not transfers,
            "scope": "Geometry screen only; no tautomer populations, bond-order reassignment, Hessian or pKa guarantee"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--evaluation", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--target", default="BACE")
    args = parser.parse_args()
    cohort = json.loads(args.registry.read_text())["cohorts"][args.target]
    frozen = json.loads(args.protocol.read_text())
    evaluation = json.loads(args.evaluation.read_text())
    molecules = {molecule.GetProp("_Name").strip(): molecule
                 for molecule in Chem.SDMolSupplier(cohort["sdf"], removeHs=False) if molecule is not None}
    if set(molecules) != set(cohort["names"]) or set(evaluation["scores"]) != set(cohort["names"]):
        raise ValueError("full chemical identity inventory required")
    rows = {}
    for name in cohort["names"]:
        molecule = molecules[name]
        reference = stereochemistry(molecule, molecule.GetConformer().GetPositions())
        bound = load_bound_ligand(cohort["sdf"], name, cohort["complexes"][name])
        free = evaluation["rows"][frozen["case_tasks"][name]["free_ligand_minimum"]]
        if free["status"] != "SUCCEEDED":
            raise ValueError("free minimum not verified")
        positions = [list(map(float, line.split()[2:])) for line in free["optimized_coordinates"]]
        rows[name] = {"input_stereochemistry_3d": reference,
            "bound_k10": geometry_screen(molecule, [atom.xyz for atom in bound.atoms], reference),
            "free_pcm_minimum": geometry_screen(molecule, positions, reference)}
    report = {"schema": "cycle-chemistry-geometry-screen-v1", "target": args.target, "n": len(rows),
        "n_flagged": sum(not all(row[state]["passed_screen"] for state in ("bound_k10", "free_pcm_minimum")) for row in rows.values()),
        "rows": rows, "driver_sha256": digest(__file__), "registry_sha256": digest(args.registry),
        "evaluation_sha256": digest(args.evaluation), "labels_used": False,
        "warning": "Passing these screens is not proof of the correct solution microstate or a true global minimum"}
    write_json(args.output, report)


if __name__ == "__main__":
    main()
