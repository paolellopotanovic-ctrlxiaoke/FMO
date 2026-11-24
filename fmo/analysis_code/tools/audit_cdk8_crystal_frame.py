"""Check a protein-aligned crystal pose with bijective element assignment."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from rdkit import Chem
from scipy.optimize import linear_sum_assignment

from run_cdk8_solution_cycle import digest, write_json


def align(reference, target):
    reference_center, target_center = reference.mean(axis=0), target.mean(axis=0)
    left, _, right = np.linalg.svd((reference-reference_center).T @ (target-target_center))
    parity = np.eye(3)
    parity[-1,-1] = np.linalg.det(left @ right)
    rotation = left @ parity @ right
    return rotation, reference_center, target_center


def main():
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--cif", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    base = args.workspace / "vendor/public_binding_free_energy_benchmark/fep_benchmark_inputs/structure_inputs/misc"
    protein, sdf = base / "cdk8_koehler_protein.pdb", base / "cdk8_koehler_ligands.sdf"
    structure = MMCIF2Dict(str(args.cif))
    reference_ca = {}
    for index, record in enumerate(structure["_atom_site.group_PDB"]):
        if record == "ATOM" and structure["_atom_site.auth_asym_id"][index] == "A" and structure["_atom_site.auth_atom_id"][index] == "CA":
            reference_ca[int(structure["_atom_site.auth_seq_id"][index])] = np.array([float(structure["_atom_site.Cartn_"+axis][index]) for axis in "xyz"])
    target_ca = {int(line[22:26]): np.array([float(line[start:stop]) for start,stop in ((30,38),(38,46),(46,54))])
        for line in protein.read_text().splitlines() if line.startswith("ATOM") and line[21] == "A" and line[12:16].strip() == "CA"}
    keys = sorted(reference_ca.keys() & target_ca.keys())
    if len(keys) < 100:
        raise ValueError("insufficient sequence-number matched CA atoms")
    reference, target = np.array([reference_ca[key] for key in keys]), np.array([target_ca[key] for key in keys])
    rotation, reference_center, target_center = align(reference, target)
    fitted = (reference-reference_center) @ rotation + target_center
    indices = [index for index, label in enumerate(structure["_atom_site.label_comp_id"]) if label == "50R"]
    crystal = np.array([[float(structure["_atom_site.Cartn_"+axis][index]) for axis in "xyz"] for index in indices])
    crystal = (crystal-reference_center) @ rotation + target_center
    symbols = [structure["_atom_site.type_symbol"][index] for index in indices]
    molecule = next(mol for mol in Chem.SDMolSupplier(str(sdf), removeHs=False) if mol.GetProp("_Name") == "22")
    positions = molecule.GetConformer().GetPositions()
    assignments = []
    for symbol in sorted(set(symbols)):
        first = [index for index,element in enumerate(symbols) if element == symbol]
        second = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetSymbol() == symbol]
        if len(first) != len(second):
            raise ValueError("ligand crystal element inventory differs")
        costs = np.linalg.norm(crystal[first,None,:] - positions[np.array(second)[None,:],:], axis=2)
        left, right = linear_sum_assignment(costs)
        assignments.extend({"crystal_atom": structure["_atom_site.auth_atom_id"][indices[first[first_index]]],
            "sdf_atom_zero_based": second[second_index], "element": symbol,
            "distance_angstrom": float(costs[first_index,second_index])} for first_index,second_index in zip(left,right,strict=True))
    result = {"reference_pdb": "5CEI", "reference_title": structure.get("_struct.title"),
        "n_matched_protein_ca": len(keys), "ca_alignment_rmsd_angstrom": float(np.sqrt(np.mean(np.sum((fitted-target)**2,axis=1)))),
        "ligand": "22", "crystal_component": "50R", "n_bijectively_matched_heavy_atoms": len(assignments),
        "ligand_heavy_rmsd_after_protein_alignment": float(np.sqrt(np.mean([record["distance_angstrom"]**2 for record in assignments]))),
        "assignments": assignments, "sources_sha256": {str(path): digest(path) for path in (protein,sdf,args.cif)},
        "warning": "Demonstrates input22 crystal coordinate consistency, not crystallographic validation of the other nine poses"}
    write_json(args.output,result)
    print(result["ca_alignment_rmsd_angstrom"],result["ligand_heavy_rmsd_after_protein_alignment"])


if __name__ == "__main__":
    main()
