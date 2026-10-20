"""Minimal structure readers used before optional RDKit integration."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import re
import numpy as np


@dataclass(frozen=True)
class Atom:
    index: int
    element: str
    xyz: tuple[float, float, float]
    chain: str = ""
    residue_name: str = ""
    residue_number: str = ""
    record: str = "ATOM"
    atom_name: str = ""


@dataclass(frozen=True)
class Ligand:
    name: str
    atoms: tuple[Atom, ...]
    smiles: str | None = None
    formal_charge: int = 0
    properties: tuple[tuple[str, str], ...] = ()


@dataclass(frozen=True)
class Protein:
    atoms: tuple[Atom, ...]


def read_pdb(path: str | Path) -> Protein:
    atoms: list[Atom] = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            record = line[:6].strip()
            if record not in {"ATOM", "HETATM"}:
                continue
            try:
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            element = line[76:78].strip() or re.sub(r"[^A-Za-z]", "", line[12:16])[:1]
            atoms.append(Atom(len(atoms), element.upper(), xyz, line[21].strip(),
                              line[17:20].strip(), line[22:26].strip(), record,
                              line[12:16].strip()))
    if not atoms:
        raise ValueError(f"No ATOM/HETATM coordinates found in {path}")
    return Protein(tuple(atoms))


def _sdf_blocks(text: str) -> Iterable[str]:
    for block in text.split("$$$$"):
        if block.strip():
            yield block.strip("\n")


def read_sdf(path: str | Path) -> list[Ligand]:
    ligands: list[Ligand] = []
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    for block in _sdf_blocks(text):
        lines = block.splitlines()
        if len(lines) < 4:
            continue
        name = lines[0].strip() or f"ligand_{len(ligands)}"
        try:
            natoms = int(lines[3][0:3])
        except ValueError as exc:
            raise ValueError(f"Cannot parse SDF counts line in {path}: {lines[3]!r}") from exc
        atoms: list[Atom] = []
        atom_charges: dict[int, int] = {}
        for index, line in enumerate(lines[4:4 + natoms]):
            try:
                xyz = (float(line[0:10]), float(line[10:20]), float(line[20:30]))
                element = line[31:34].strip().upper()
            except (ValueError, IndexError) as exc:
                raise ValueError(f"Cannot parse SDF atom line in {path}: {line!r}") from exc
            atoms.append(Atom(index, element, xyz, record="HETATM"))
            charge_code = line[36:39].strip() if len(line) >= 39 else ""
            charge = {"1": 3, "2": 2, "3": 1, "5": -1, "6": -2, "7": -3}.get(charge_code, 0)
            if charge:
                atom_charges[index + 1] = charge
        properties: list[tuple[str, str]] = []
        for index, line in enumerate(lines):
            if line.startswith("M  CHG"):
                fields = line.split()
                count = int(fields[2])
                for pair_index in range(count):
                    atom_id = int(fields[3 + 2 * pair_index])
                    charge = int(fields[4 + 2 * pair_index])
                    atom_charges[atom_id] = charge
            elif line.startswith("> <") and ">" in line[3:]:
                prop_name = line[3:line.index(">", 3)].strip()
                values = []
                for value_line in lines[index + 1:]:
                    if not value_line.strip():
                        break
                    values.append(value_line.strip())
                properties.append((prop_name, "\n".join(values)))
        ligands.append(Ligand(name, tuple(atoms), formal_charge=sum(atom_charges.values()),
                              properties=tuple(properties)))
    if not ligands:
        raise ValueError(f"No SDF molecules found in {path}")
    return ligands


def residue_key(atom: Atom) -> str:
    return f"{atom.chain}:{atom.residue_name}:{atom.residue_number}"


def select_pocket(protein: Protein, ligand: Ligand, cutoff: float = 5.0) -> Protein:
    if cutoff <= 0:
        raise ValueError("cutoff must be positive")
    protein_xyz = np.asarray([atom.xyz for atom in protein.atoms], dtype=float)
    ligand_xyz = np.asarray([atom.xyz for atom in ligand.atoms], dtype=float)
    distances = np.linalg.norm(protein_xyz[:, None, :] - ligand_xyz[None, :, :], axis=-1)
    residue_minimum: dict[str, float] = {}
    for atom, row in zip(protein.atoms, distances, strict=True):
        if atom.record != "ATOM" and atom.residue_name not in {"HOH", "WAT", "TIP"}:
            continue
        key = residue_key(atom)
        residue_minimum[key] = min(residue_minimum.get(key, float("inf")), float(row.min()))
    selected_keys = {key for key, minimum in residue_minimum.items() if minimum <= cutoff}
    selected = [
        atom for atom in protein.atoms
        if residue_key(atom) in selected_keys
        and (atom.record == "ATOM" or atom.residue_name in {"HOH", "WAT", "TIP"})
    ]
    if not selected:
        raise ValueError(f"No protein atoms within {cutoff} Angstrom of {ligand.name}")
    return Protein(tuple(selected))


def bridge_one_residue_sequence_gaps(protein: Protein, reference: Protein) -> tuple[Protein, list[str]]:
    """Retain enclosed single-residue chain gaps before ACE/NME capping."""
    selected = {residue_key(atom) for atom in protein.atoms if atom.record == "ATOM"}
    reference_residues: dict[str, list[Atom]] = {}
    for atom in reference.atoms:
        if atom.record == "ATOM":
            reference_residues.setdefault(residue_key(atom), []).append(atom)

    def residue_number(key: str) -> int | None:
        try:
            return int(key.rsplit(":", 1)[-1])
        except ValueError:
            return None

    by_chain: dict[str, list[int]] = {}
    for key in selected:
        value = residue_number(key)
        if value is not None:
            by_chain.setdefault(key.split(":", 1)[0], []).append(value)
    bridges: list[str] = []
    for chain, values in by_chain.items():
        ordered = sorted(set(values))
        for left, right in zip(ordered, ordered[1:]):
            if right - left != 2:
                continue
            matches = [key for key in reference_residues
                       if key.startswith(f"{chain}:") and residue_number(key) == left + 1]
            if len(matches) != 1:
                continue
            selected.add(matches[0])
            bridges.append(matches[0])

    atoms = list(protein.atoms)
    known = {residue_key(atom) for atom in atoms}
    for key, residue_atoms in reference_residues.items():
        if key in selected and key not in known:
            atoms.extend(residue_atoms)
    return Protein(tuple(atoms)), sorted(bridges)

