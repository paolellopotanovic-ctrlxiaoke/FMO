#!/usr/bin/env python
"""Validate bound protein/protein, protein/RNA, or RNA/ligand complex geometry."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from rdkit import Chem
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

from .target_intake import digest


PROTEIN = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
           "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
RNA = {"A", "U", "G", "C", "RA", "RU", "RG", "RC"}
WATER = {"HOH", "WAT", "H2O"}
TYPES = {"protein_protein": ("protein", "protein"),
         "protein_rna": ("protein", "rna"), "rna_ligand": ("rna", "ligand")}


def parse_pdb(path: Path, chain: str) -> list[dict]:
    atoms = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line[:6].strip() not in ("ATOM", "HETATM") or line[21].strip() != chain:
            continue
        if line[16].strip():
            raise ValueError(f"alternate locations need explicit review in chain {chain}")
        element = line[76:78].strip().upper()
        if not element:
            raise ValueError(f"missing element field in chain {chain}")
        if element in ("H", "D") or line[17:20].strip() in WATER:
            continue
        coordinates = tuple(float(line[position:position + 8]) for position in (30, 38, 46))
        if not all(map(math.isfinite, coordinates)):
            raise ValueError("non-finite bound complex coordinates")
        atoms.append({"record": line[:6].strip(), "resname": line[17:20].strip().upper(),
                      "element": element, "coordinates": coordinates,
                      "residue": line[22:27].strip(), "atom": line[12:16].strip()})
    if not atoms:
        raise ValueError(f"no heavy atoms in chain {chain}")
    return atoms


def validate_partner(atoms: list[dict], partner_type: str, ligand_sdf: Path | None) -> dict:
    residues = {atom["resname"] for atom in atoms}
    if partner_type == "protein":
        if not residues <= PROTEIN or any(atom["record"] != "ATOM" for atom in atoms):
            raise ValueError(f"protein contains unsupported residues or HETATM: {sorted(residues - PROTEIN)}")
        if len(atoms) < 10 or len({atom["residue"] for atom in atoms}) < 3:
            raise ValueError("protein partner appears incomplete (<10 heavy atoms or <3 residues)")
    elif partner_type == "rna":
        if not residues <= RNA or any(atom["record"] != "ATOM" for atom in atoms):
            raise ValueError(f"RNA contains unsupported residues or HETATM: {sorted(residues - RNA)}")
        if len(atoms) < 5 or len({atom["residue"] for atom in atoms}) < 2:
            raise ValueError("RNA partner appears incomplete (<5 heavy atoms or <2 residues)")
    elif partner_type == "ligand":
        if ligand_sdf is None:
            raise ValueError("ligand requires a validated bound 3D SDF; PDB alone lacks chemistry")
        if (len(residues) != 1 or len({atom["residue"] for atom in atoms}) != 1 or
                any(atom["record"] != "HETATM" for atom in atoms)):
            raise ValueError("ligand must be one named HETATM residue in its own chain")
        molecules = list(Chem.SDMolSupplier(str(ligand_sdf), removeHs=False))
        if len(molecules) != 1 or molecules[0] is None or not molecules[0].GetNumConformers():
            raise ValueError("ligand SDF must contain exactly one bound 3D chemical structure")
        molecule = molecules[0]
        if not molecule.GetConformer().Is3D():
            raise ValueError("ligand SDF must be three-dimensional")
        heavy = [atom for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1]
        if len(heavy) != len(atoms):
            raise ValueError("ligand heavy atom count differs between PDB and SDF")
        for element in {atom.GetSymbol().upper() for atom in heavy} | {
                atom["element"] for atom in atoms}:
            sdf_indices = [atom.GetIdx() for atom in heavy if atom.GetSymbol().upper() == element]
            pdb_atoms = [atom for atom in atoms if atom["element"] == element]
            if len(sdf_indices) != len(pdb_atoms):
                raise ValueError("ligand heavy element composition differs between PDB and SDF")
            sdf_positions = molecule.GetConformer().GetPositions()[sdf_indices]
            pdb_positions = np.asarray([atom["coordinates"] for atom in pdb_atoms])
            distances = cdist(sdf_positions, pdb_positions)
            first, second = linear_sum_assignment(distances)
            if np.max(distances[first, second]) > .25:
                raise ValueError("ligand SDF coordinates not aligned with the bound complex PDB")
        return {"formal_charge": Chem.GetFormalCharge(molecule),
                "canonical_smiles": Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=True)}
    return {"n_residues": len({atom["residue"] for atom in atoms})}


def audit(complex_pdb: Path, kind: str, chain_a: str, chain_b: str,
          ligand_sdf: Path | None = None) -> dict:
    if kind not in TYPES or not chain_a or not chain_b or chain_a == chain_b or (
            kind != "rna_ligand" and ligand_sdf is not None):
        raise ValueError("distinct partner chains and a valid interface type are required")
    first, second = parse_pdb(complex_pdb, chain_a), parse_pdb(complex_pdb, chain_b)
    type_a, type_b = TYPES[kind]
    chemistry_a = validate_partner(first, type_a, None)
    chemistry_b = validate_partner(second, type_b, ligand_sdf)
    distances = cdist([item["coordinates"] for item in first],
                      [item["coordinates"] for item in second])
    nearest = float(np.min(distances))
    if nearest < 1.0 or nearest > 12.0:
        raise ValueError(f"unusable bound interface minimum heavy-atom distance {nearest:.3f} Å")
    return {"schema": "fmo-molecular-interface-intake-v1", "interface_type": kind,
            "complex_sha256": digest(complex_pdb),
            "ligand_sdf_sha256": digest(ligand_sdf) if ligand_sdf is not None else None,
            "partner_a": {"type": type_a, "chain": chain_a, "n_heavy_atoms": len(first), **chemistry_a},
            "partner_b": {"type": type_b, "chain": chain_b, "n_heavy_atoms": len(second), **chemistry_b},
            "minimum_cross_heavy_atom_distance_angstrom": nearest,
            "status": "BOUND_GEOMETRY_AUDITED_NEEDS_QM_FRAGMENT_MAPPING",
            "ranking_eligible": False,
            "next_steps": ["validate partner protonation, ions, waters, nucleic acid caps and QM method",
                           "construct a complete AB GAMESS FMO deck with explicit fragment-to-partner map",
                           "only then score all A–B PIEDA pairs with cli.py interface"]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--complex", type=Path, required=True)
    parser.add_argument("--interface-type", choices=tuple(TYPES), required=True)
    parser.add_argument("--chain-a", required=True)
    parser.add_argument("--chain-b", required=True)
    parser.add_argument("--ligand-sdf", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.complex, args.interface_type, args.chain_a, args.chain_b,
                   args.ligand_sdf)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": report["status"], "ranking_eligible": False}))


if __name__ == "__main__":
    main()
