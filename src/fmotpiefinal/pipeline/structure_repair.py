#!/usr/bin/env python
"""Run a label-free PDBFixer audit and optional hydrogen completion."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from pdbfixer import PDBFixer
from openmm.app import PDBFile


def _item_name(item) -> str:
    return str(getattr(item, "name", item))


def prepare(source: Path, output: Path, report_path: Path, ph: float) -> dict:
    fixer = PDBFixer(filename=str(source))
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    missing_residues = {
        str(key): [res.name for res in value]
        for key, value in fixer.missingResidues.items()
    }
    nonstandard = [
        {"residue": residue.name, "replacement": replacement}
        for residue, replacement in fixer.nonstandardResidues
    ]
    before_atoms = fixer.topology.getNumAtoms()
    if fixer.nonstandardResidues:
        fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    missing_atoms = {
        str(key): [_item_name(atom) for atom in value]
        for key, value in fixer.missingAtoms.items()
    }
    missing_terminals = {
        str(key): [_item_name(atom) for atom in value]
        for key, value in fixer.missingTerminals.items()
    }
    if fixer.missingResidues or fixer.missingAtoms or fixer.missingTerminals:
        # PDBFixer adds missing residues and terminal atoms through addMissingAtoms.
        fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=ph)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
    after_atoms = fixer.topology.getNumAtoms()
    report = {
        "source": str(source),
        "output": str(output),
        "pH": ph,
        "missing_residues": missing_residues,
        "missing_atoms": missing_atoms,
        "missing_terminals": missing_terminals,
        "nonstandard_residues": nonstandard,
        "input_atom_count": before_atoms,
        "pdbfixer_hydrogen_completed_atom_count": after_atoms,
        "added_atom_count": after_atoms - before_atoms,
        "policy": "audit then repair missing residues/atoms/nonstandard residues and add hydrogens at pH; no label-dependent edits",
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--ph", type=float, default=7.4)
    args = parser.parse_args()
    print(json.dumps(prepare(args.source, args.output, args.report, args.ph), indent=2))


if __name__ == "__main__":
    main()
