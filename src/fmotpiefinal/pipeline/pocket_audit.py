#!/usr/bin/env python
"""Audit whole-target all-atom versus heavy-atom shared-pocket membership."""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path


from fmotpiefinal.core.bound_geometry import load_bound_ligand, shared_pocket
from fmotpiefinal.core.structure import Ligand, Protein, read_pdb, residue_key, select_pocket

from .target_intake import digest


def audit_pockets(protein: Protein, ligands: list[Ligand], cutoff: float = 5.0) -> dict:
    if cutoff <= 0 or not ligands or len({ligand.name for ligand in ligands}) != len(ligands):
        raise ValueError("positive cutoff and a unique complete ligand cohort are required")
    all_atom = shared_pocket(protein, ligands, cutoff, include_waters=True)
    heavy = shared_pocket(protein, ligands, cutoff, include_waters=True, heavy_only=True)
    full_membership = {residue_key(atom) for atom in all_atom.atoms}
    heavy_membership = {residue_key(atom) for atom in heavy.atoms}
    direct_all = {residue_key(atom) for ligand in ligands
                  for atom in select_pocket(protein, ligand, cutoff).atoms}
    direct_heavy = {residue_key(atom) for ligand in ligands
                    for atom in select_pocket(protein, ligand, cutoff, heavy_only=True).atoms}
    differences = []
    for key in sorted(full_membership ^ heavy_membership):
        atoms = [atom for atom in protein.atoms if residue_key(atom) == key]
        distances = {}
        for name, receptor, ligand_filter in (
                ("all_atom", atoms, lambda atom: True),
                ("heavy_atom", [atom for atom in atoms if atom.element.upper() not in {"H", "D"}],
                 lambda atom: atom.element.upper() not in {"H", "D"})):
            candidates = [ligand_atom for ligand in ligands for ligand_atom in ligand.atoms
                          if ligand_filter(ligand_atom)]
            distances[name] = min((math.dist(atom.xyz, ligand_atom.xyz) for atom in receptor
                                    for ligand_atom in candidates), default=None)
        differences.append({"residue": key, "selected_all_atom": key in full_membership,
                            "selected_heavy_atom": key in heavy_membership,
                            "direct_contact_all_atom": key in direct_all,
                            "direct_contact_heavy_atom": key in direct_heavy,
                            "nearest_all_atom_angstrom": distances["all_atom"],
                            "nearest_heavy_atom_angstrom": distances["heavy_atom"]})
    return {"schema": "fmo-whole-target-pocket-audit-v1", "label_free": True,
            "cutoff_angstrom": cutoff, "n_ligands": len(ligands),
            "n_receptor_atoms_all_atom": len(all_atom.atoms),
            "n_receptor_atoms_heavy_atom": len(heavy.atoms),
            "n_residues_all_atom": len(full_membership),
            "n_residues_heavy_atom": len(heavy_membership),
            "membership_differences": differences}


def audit_run(run_dir: Path, output: Path, names: list[str], cutoff: float = 5.0) -> dict:
    target = run_dir / "frozen_protocol.json"
    protocol = json.loads(target.read_text())
    if protocol.get("ligand_names") != names:
        raise ValueError("frozen target cohort differs from requested pocket audit")
    protein_path = run_dir / f"{protocol['target_name'].lower()}_prepared.pdb"
    sdf_path = run_dir / "ligands_no_properties.sdf"
    protein = read_pdb(protein_path)
    ligands = [load_bound_ligand(sdf_path, name,
               run_dir / "minimized" / name / "minimized_complex.pdb") for name in names]
    report = audit_pockets(protein, ligands, cutoff)
    report["source_sha256"] = {"receptor": digest(protein_path), "ligands": digest(sdf_path),
                                "frozen_protocol": digest(target)}
    report["selection_production"] = "heavy_atom" if protocol.get("pocket_heavy_only") else "all_atom"
    report["ligand_names"] = names
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cutoff", type=float, default=5.)
    args = parser.parse_args()
    protocol = json.loads((args.run_dir / "frozen_protocol.json").read_text())
    if args.cutoff != protocol.get("cutoff_angstrom", 5.):
        raise ValueError("audit cutoff differs from the frozen scoring protocol")
    report = audit_run(args.run_dir, args.output, protocol["ligand_names"], args.cutoff)
    print(json.dumps({"n_ligands": report["n_ligands"],
                      "membership_differences": len(report["membership_differences"])}))


if __name__ == "__main__":
    main()
