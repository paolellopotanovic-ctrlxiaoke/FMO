"""Bound-pose loading and full-cohort shared-pocket selection."""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from .structure import (Ligand, Protein, read_pdb, read_sdf, residue_key,
                        select_pocket, bridge_one_residue_sequence_gaps)


def load_bound_ligand(sdf: str | Path, name: str, minimized_complex: str | Path) -> Ligand:
    ligand = next((item for item in read_sdf(sdf) if item.name == name), None)
    if ligand is None:
        raise ValueError(f"ligand {name!r} not found in {sdf}")
    complex_atoms = [atom for atom in read_pdb(minimized_complex).atoms
                     if atom.record == "HETATM" and atom.residue_name == "UNK"]
    if len(complex_atoms) != len(ligand.atoms) or any(
        left.element.upper() != right.element.upper()
        for left, right in zip(ligand.atoms, complex_atoms, strict=True)
    ):
        raise ValueError(f"complex pose does not match ligand SDF atom order: {name}")
    return replace(ligand, atoms=tuple(
        replace(atom, xyz=pose.xyz)
        for atom, pose in zip(ligand.atoms, complex_atoms, strict=True)
    ))


def shared_pocket(protein: Protein, ligands: list[Ligand], cutoff: float, *,
                  include_waters: bool = True, heavy_only: bool = False) -> Protein:
    """Union residues/waters over a target's entire predefined cohort."""
    if not ligands:
        raise ValueError("target must contain at least one ligand")
    selected = set()
    for ligand in ligands:
        selected.update(residue_key(atom) for atom in select_pocket(
            protein, ligand, cutoff, heavy_only=heavy_only).atoms)
    pocket = Protein(tuple(atom for atom in protein.atoms
                           if residue_key(atom) in selected and
                           (atom.record == "ATOM" or
                            (include_waters and atom.residue_name in {"HOH", "WAT", "TIP"}))))
    pocket, _ = bridge_one_residue_sequence_gaps(pocket, protein)
    return pocket
