"""Create a frozen, label-free sidechain decision graph for a prepared target."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from rdkit import Chem

from .core.bound_geometry import load_bound_ligand
from .core.sidechain_policy import explicit_nonpeptide_links, select_sidechains
from .core.structure import read_pdb


def audit_sidechains(protein_path, sdf_path, complexes):
    protein_path, sdf_path = Path(protein_path), Path(sdf_path)
    molecules = {}
    for molecule in Chem.SDMolSupplier(str(sdf_path), removeHs=False):
        if molecule is None:
            raise ValueError("invalid ligand record")
        name = molecule.GetProp("_Name").strip()
        if name in molecules:
            raise ValueError("duplicate ligand identity")
        molecules[name] = molecule
    if set(molecules) != set(complexes):
        raise ValueError("complete original/bound inventory mismatch")
    bound = {}
    for name in molecules:
        ligand = load_bound_ligand(sdf_path, name, Path(complexes[name]))
        if [atom.element.upper() for atom in ligand.atoms] != [atom.GetSymbol().upper() for atom in molecules[name].GetAtoms()]:
            raise ValueError("ligand atom identity/order mismatch")
        bound[name] = np.asarray([atom.xyz for atom in ligand.atoms])
    plan = select_sidechains(read_pdb(protein_path), molecules, bound,
        covalent_links=explicit_nonpeptide_links(protein_path))
    plan["source_sha256"] = {str(Path(path).resolve()): hashlib.sha256(Path(path).read_bytes()).hexdigest()
        for path in (protein_path, sdf_path, *complexes.values())}
    plan.pop("plan_sha256")
    plan["plan_sha256"] = hashlib.sha256(json.dumps(plan, sort_keys=True, allow_nan=False).encode()).hexdigest()
    return plan


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protein", required=True, type=Path)
    parser.add_argument("--poses", required=True, type=Path)
    parser.add_argument("--bound-complexes", required=True, type=Path,
                        help="JSON mapping ligand names to baseline k10 complex PDB paths")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    plan = audit_sidechains(args.protein, args.poses, json.loads(args.bound_complexes.read_text()))
    text = json.dumps(plan, indent=2, allow_nan=False) + "\n"
    if args.output.exists() and args.output.read_text() != text:
        raise ValueError("frozen target plan differs; use a versioned new output")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text)
    print(json.dumps({"action": plan["action"], "selected_residues": plan["selected_residues"]}))


if __name__ == "__main__":
    main()
