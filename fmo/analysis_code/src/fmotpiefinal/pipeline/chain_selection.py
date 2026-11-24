#!/usr/bin/env python
"""Select one receptor chain from heavy-atom contacts, without affinity labels."""
from __future__ import annotations

from pathlib import Path
import math

import numpy as np
from rdkit import Chem
from scipy.spatial import cKDTree


def choose_receptor_chain(protein: Path, molecule: Chem.Mol) -> tuple[str, dict[str, float]]:
    if molecule.GetNumConformers() != 1:
        raise ValueError("crystal ligand requires exactly one bound conformer")
    conformer = molecule.GetConformer()
    ligand = np.array([list(conformer.GetAtomPosition(atom.GetIdx()))
                       for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1])
    if not len(ligand) or not np.isfinite(ligand).all():
        raise ValueError("crystal ligand has no valid heavy-atom coordinates")
    chains: dict[str, list[tuple[float, float, float]]] = {}
    residues: dict[str, dict[str, str]] = {}
    for line in protein.read_text().splitlines():
        if not line.startswith("ATOM  ") or line[76:78].strip().upper() in ("H", "D"):
            continue
        chains.setdefault(line[21], []).append(tuple(float(line[position:position + 8])
                                                  for position in (30, 38, 46)))
        residues.setdefault(line[21], {})[line[22:27]] = line[17:20]
    if not chains:
        raise ValueError("crystal receptor contains no protein heavy atoms")
    distances = {chain: round(float(cKDTree(coords).query(ligand)[0].min()), 4)
                 for chain, coords in sorted(chains.items())}
    touching = [chain for chain, distance in distances.items() if distance < 5.0]
    if (len(touching) == 2 and len(chains) == 2 and
            all(len(residues[chain]) >= 50 for chain in touching) and
            list(residues[touching[0]].values()) == list(residues[touching[1]].values()) and
            cKDTree(chains[touching[0]]).count_neighbors(cKDTree(chains[touching[1]]), 4.0) >= 80):
        return ",".join(sorted(touching)), distances
    if len(touching) != 1:
        raise ValueError(f"expected exactly one receptor chain within 5 A: {distances}")
    return touching[0], distances


def choose_cohort_chain(protein: Path, molecules: list[Chem.Mol]) -> str:
    if not molecules:
        raise ValueError("missing crystal poses for automatic receptor-chain selection")
    selected = {choose_receptor_chain(protein, molecule)[0] for molecule in molecules}
    if len(selected) != 1:
        raise ValueError(f"cohort crystal poses contact different receptor chains: {sorted(selected)}")
    return selected.pop()


def covalently_linked_chains(protein: Path, contact_chain: str) -> tuple[str, list[dict]]:
    """Keep chains bonded by an observed inter-chain cystine, not nearby crystal copies."""
    sulfurs: dict[str, list[tuple[str, tuple[float, float, float]]]] = {}
    for line in protein.read_text().splitlines():
        if line.startswith("ATOM  ") and line[12:16].strip() == "SG" and line[17:20] == "CYS":
            sulfurs.setdefault(line[21], []).append((line[22:27].strip(), tuple(
                float(line[position:position + 8]) for position in (30, 38, 46))))
    links = []
    for chain, atoms in sorted(sulfurs.items()):
        for partner, others in sorted(sulfurs.items()):
            if partner <= chain:
                continue
            for residue, xyz in atoms:
                for other_residue, other_xyz in others:
                    distance = math.dist(xyz, other_xyz)
                    if 1.8 <= distance <= 2.3:
                        links.append({"chain_a": chain, "residue_a": residue,
                                      "chain_b": partner, "residue_b": other_residue,
                                      "sg_distance_angstrom": round(distance, 4)})
    selected = set(contact_chain.split(","))
    while True:
        attached = selected | {chain for link in links for chain, partner in (
            (link["chain_a"], link["chain_b"]), (link["chain_b"], link["chain_a"]))
            if partner in selected}
        if selected == attached:
            break
        selected = attached
    return ",".join(sorted(selected)), [link for link in links if
           link["chain_a"] in selected and link["chain_b"] in selected]
