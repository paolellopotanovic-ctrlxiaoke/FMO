#!/usr/bin/env python
"""Freeze a label-independent, shared-pocket FMO affinity ranking protocol."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from rdkit import Chem
from scipy.spatial import cKDTree


REQUIRED = ("target", "ligand_name", "ligand_sdf", "protein_pdb", "protein_chain",
            "canonical_isomeric_smiles", "formal_charge")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def structural_sdf_digest(entries: dict[str, Chem.Mol]) -> str:
    structures = [(name, Chem.MolToMolBlock(molecule)) for name, molecule in sorted(entries.items())]
    return hashlib.sha256(json.dumps(structures, ensure_ascii=False).encode("utf-8")).hexdigest()


def protein_heavy_atoms(path: Path, chain: str) -> np.ndarray:
    chains = {value.strip() for value in chain.split(",") if value.strip()}
    if not chains:
        raise ValueError("protein_chain must name at least one chain")
    coordinates = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("ATOM  ") or line[21].strip() not in chains:
                continue
            if line[76:78].strip().upper() in ("H", "D") or line[12:16].strip().startswith("H"):
                continue
            coordinates.append([float(line[start:start + 8]) for start in (30, 38, 46)])
    if not coordinates or not np.isfinite(coordinates).all():
        raise ValueError(f"missing finite protein heavy atoms in {path} chain(s) {chain}")
    return np.asarray(coordinates)


def check_structure(row: dict, molecule: Chem.Mol) -> list[str]:
    failures = []
    if Chem.GetFormalCharge(molecule) != int(row["formal_charge"]):
        failures.append("cohort_sdf_charge_mismatch")
    if not any(atom.GetAtomicNum() == 1 for atom in molecule.GetAtoms()):
        failures.append("missing_explicit_hydrogens")
    if molecule.GetNumConformers() != 1 or not molecule.GetConformer().Is3D():
        failures.append("missing_3d_geometry")
    else:
        if not all(math.isfinite(value) for point in molecule.GetConformer().GetPositions()
                   for value in point):
            failures.append("nonfinite_coordinates")
    if Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=True) != row["canonical_isomeric_smiles"]:
        failures.append("cohort_sdf_identity_mismatch")
    return failures


def freeze(cohort: Path, output: Path) -> dict:
    with cohort.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not set(REQUIRED).issubset(reader.fieldnames or []):
            raise ValueError("missing required structural cohort columns")
        records = [{key: row[key] for key in REQUIRED} for row in reader]
    if not records or len({(row["target"], row["ligand_name"]) for row in records}) != len(records):
        raise ValueError("missing or duplicate case identities")
    sources = sorted({row["ligand_sdf"] for row in records})
    molecules = {}
    for source in sources:
        path = Path(source)
        if not path.is_file():
            raise FileNotFoundError(path)
        entries = {}
        for molecule in Chem.SDMolSupplier(str(path), removeHs=False):
            if molecule is None:
                raise ValueError(f"unreadable molecule in {source}")
            name = molecule.GetProp("_Name")
            if name in entries:
                raise ValueError(f"duplicate source ligand name: {source}/{name}")
            entries[name] = molecule
        molecules[source] = entries
    target_proteins = {}
    protein_coordinates = {}
    cases = []
    for record in records:
        target, name = record["target"], record["ligand_name"]
        protein = Path(record["protein_pdb"])
        if not protein.is_file():
            raise FileNotFoundError(protein)
        identity = (record["protein_pdb"], record["protein_chain"])
        if target in target_proteins and target_proteins[target] != identity:
            raise ValueError(f"multiple receptor inputs in one shared-pocket target: {target}")
        target_proteins[target] = identity
        if target not in protein_coordinates:
            protein_coordinates[target] = cKDTree(protein_heavy_atoms(protein, record["protein_chain"]))
        molecule = molecules[record["ligand_sdf"]].get(name)
        failures = (["missing_sdf_name"] if molecule is None else check_structure(record, molecule))
        distance = None
        if molecule is not None and molecule.GetNumConformers() == 1 and molecule.GetConformer().Is3D():
            indices = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1]
            if indices:
                minimum = protein_coordinates[target].query(molecule.GetConformer().GetPositions()[indices])[0]
                distance = float(min(minimum))
                if not math.isfinite(distance) or distance < 1.0:
                    failures.append("severe_protein_ligand_overlap")
                elif distance > 8.0:
                    failures.append("ligand_outside_protein_binding_site")
        cases.append({"target": target, "ligand": name, "ligand_sdf": record["ligand_sdf"],
                      "formal_charge": (None if molecule is None else int(Chem.GetFormalCharge(molecule))),
                      "minimum_source_heavy_contact_angstrom": distance,
                      "quality_gate": "PASS" if not failures else "REVIEW",
                      "reasons": failures})
    protein_hashes = {target: {"path": path, "chain": chain, "sha256": digest(Path(path))}
                      for target, (path, chain) in sorted(target_proteins.items())}
    report = {"n": len(cases), "n_ready": sum(row["quality_gate"] == "PASS" for row in cases),
              "n_review": sum(row["quality_gate"] != "PASS" for row in cases),
              "structural_cohort_sha256": hashlib.sha256(json.dumps(
                  records, sort_keys=True, ensure_ascii=False).encode("utf-8")).hexdigest(),
              "source_sdf_structures_sha256": {path: structural_sdf_digest(molecules[path])
                                               for path in sources},
              "protein_inputs": protein_hashes,
              "protocol": {
                  "case_selection": "all cases in each target, before shared-pocket selection",
                  "protein": "PDB2PQR pH 7.4 AMBER; fixed receptor coordinates and waters",
                  "pose": "OpenMM/OpenFF fixed-receptor restrained single-minimum, k=10 kJ/mol/A^2",
                  "pocket": "5.0 A target-wide union across all cohort ligand poses",
                  "qm": "GAMESS FMO2-DFTB3/PCM selected ligand-receptor PIEDA",
                  "score": "sum of ligand-receptor ES and SOLV from the same AB PIEDA log, coefficient 1",
                  "ligand_state": "preserve the explicit-H protonation and formal charge in each input SDF; no automatic state changes",
                  "mixed_charge_interpretation": "scores across charge states are diagnostic unless solution-state protonation populations are modeled at the assay pH",
                  "orientation": "smaller energy predicts stronger binding",
                  "water": "retain prepared crystallographic waters; never select water policy by label",
                  "fallback": "stop/review entire target on failure; never switch parameters for successful subset",
              },
              "cases": cases,
              "label_column_used": False,
              "ready_for_full_target_run": len(cases) == sum(row["quality_gate"] == "PASS" for row in cases)}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path(
        "protocol.json"))
    args = parser.parse_args()
    result = freeze(args.cohort, args.output)
    print(json.dumps({key: result[key] for key in ("n", "n_ready", "n_review", "ready_for_full_target_run")}))
    if not result["ready_for_full_target_run"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
