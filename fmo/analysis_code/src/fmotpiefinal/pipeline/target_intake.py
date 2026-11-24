#!/usr/bin/env python
"""Create a label-blind, fail-closed intake for new protein–ligand cohorts."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors

from .chain_selection import choose_cohort_chain, covalently_linked_chains
from .protocol import freeze
from fmotpiefinal.core.component_inventory import enforce_inventory


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_smiles(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not {"ligand_id", "smiles"}.issubset(reader.fieldnames or []):
            raise ValueError("CSV must have ligand_id and smiles columns")
        forbidden = {field for field in reader.fieldnames or [] if any(
            part in field.lower() for part in ("affinity", "experimental", "pki", "pkd", "delta_g"))}
        if forbidden:
            raise ValueError(f"label fields are not permitted in blind intake: {sorted(forbidden)}")
        rows = list(reader)
    if not rows or any(not row["ligand_id"].strip() or not row["smiles"].strip() for row in rows):
        raise ValueError("ligand ID and SMILES cannot be empty")
    if len({row["ligand_id"].strip() for row in rows}) != len(rows):
        raise ValueError("duplicate ligand IDs")
    ligands = []
    for row in rows:
        molecule = Chem.MolFromSmiles(row["smiles"].strip())
        if molecule is None or molecule.GetNumHeavyAtoms() < 1:
            raise ValueError(f"invalid SMILES for {row['ligand_id']}")
        ligands.append({"ligand_id": row["ligand_id"].strip(),
                        "canonical_smiles": Chem.MolToSmiles(molecule, isomericSmiles=True),
                        "formal_charge": Chem.GetFormalCharge(molecule),
                        "molecular_weight_g_mol": Descriptors.MolWt(molecule)})
    return ligands


def audit(target: str, ligands_csv: Path, output: Path, protein: Path | None = None,
          chain: str = "auto", poses: Path | None = None, component_decisions: Path | None = None) -> dict:
    output = output.resolve()
    if not target.strip() or not chain.strip() or poses is not None and protein is None:
        raise ValueError("target/chain must be nonempty; poses require a target protein structure")
    ligands = parse_smiles(ligands_csv)
    if protein is not None and not protein.is_file():
        raise FileNotFoundError(protein)
    if poses is not None and not poses.is_file():
        raise FileNotFoundError(poses)
    if protein is not None:
        enforce_inventory(protein, output / "raw_component_inventory.json", component_decisions)
    issues = []
    if protein is None:
        issues.append("TARGET_STRUCTURE_REQUIRED: provide experimental or validated predicted receptor PDB")
    if poses is None:
        issues.append("BOUND_POSES_REQUIRED: SMILES cannot define receptor-frame binding poses")
    if protein is None or poses is None:
        return {"schema": "fmo-new-target-intake-v1", "target": target, "n_ligands": len(ligands),
                "status": "NEEDS_STRUCTURES", "ranking_eligible": False,
                "ligand_csv_sha256": digest(ligands_csv),
                "protein_sha256": digest(protein) if protein is not None else None,
                "poses_sha256": None, "ligands": ligands, "blocking_issues": issues,
                "next_steps": ["Obtain an aligned receptor PDB and independently generated bound 3D pose per ligand",
                               "Verify atom identity, protonation, water policy and stereochemistry before k10/FMO"]}
    molecules = list(Chem.SDMolSupplier(str(poses), removeHs=False))
    if not molecules or any(molecule is None or not molecule.HasProp("_Name") for molecule in molecules):
        raise ValueError("SDF contains unreadable or unnamed 3D poses")
    mapping = {molecule.GetProp("_Name").strip(): molecule for molecule in molecules}
    if len(mapping) != len(molecules) or set(mapping) != {item["ligand_id"] for item in ligands}:
        raise ValueError("SDF identities must exactly match all ligand IDs with no duplicates")
    for item in ligands:
        name = item["ligand_id"]
        molecule = mapping[name]
        canonical = Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=True)
        if canonical != item["canonical_smiles"] or Chem.GetFormalCharge(molecule) != item["formal_charge"]:
            raise ValueError(f"3D ligand identity or formal charge differs from SMILES: {name}")
    detected_chain = choose_cohort_chain(protein, list(mapping.values()))
    if chain == "auto":
        chain = detected_chain
    elif chain != detected_chain:
        raise ValueError(f"selected receptor chain {chain} differs from bound-pose contact chain {detected_chain}")
    chain, covalent_links = covalently_linked_chains(protein, detected_chain)
    output.mkdir(parents=True, exist_ok=True)
    cohort = output / "structural_cohort.csv"
    columns = ("target", "ligand_name", "ligand_sdf", "protein_pdb", "protein_chain",
               "canonical_isomeric_smiles", "formal_charge")
    with cohort.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for item in ligands:
            writer.writerow({"target": target, "ligand_name": item["ligand_id"],
                             "ligand_sdf": str(poses.resolve()), "protein_pdb": str(protein.resolve()),
                             "protein_chain": chain, "canonical_isomeric_smiles": item["canonical_smiles"],
                             "formal_charge": item["formal_charge"]})
    preflight = freeze(cohort, output / "structure_preflight.json")
    return {"schema": "fmo-new-target-intake-v1", "target": target, "n_ligands": len(ligands),
            "status": "STRUCTURE_PREFLIGHT_PASS" if preflight["ready_for_full_target_run"]
            else "STRUCTURE_REVIEW_REQUIRED", "ranking_eligible": False,
            "ligand_csv_sha256": digest(ligands_csv), "protein_sha256": digest(protein),
            "poses_sha256": digest(poses), "cohort": str(cohort),
            "contact_chain": detected_chain, "protein_chains": chain,
            "interchain_disulfides": covalent_links,
            "n_structure_ready": preflight["n_ready"], "n_structure_review": preflight["n_review"],
            "ligands": ligands, "blocking_issues": [] if preflight["ready_for_full_target_run"] else [
                "Source ligand/receptor geometry or explicit H failed structural preflight"],
            "next_steps": ["Run validated target-level PDBFixer/PDB2PQR and water audit",
                           "Only after all poses pass, run common pocket, fixed k10 and audited full-cohort FMO"]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", required=True)
    parser.add_argument("--ligands", type=Path, required=True, help="CSV: ligand_id,smiles; no labels")
    parser.add_argument("--protein", type=Path)
    parser.add_argument("--component-decisions", type=Path)
    parser.add_argument("--chain", default="auto",
                        help="auto: unique protein chain within 5 A of every bound ligand; or explicit chain")
    parser.add_argument("--poses", type=Path, help="aligned 3D bound SDF with matching _Name IDs")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.target, args.ligands, args.output, args.protein, args.chain, args.poses,
                   component_decisions=args.component_decisions)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "intake.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "n_ligands": result["n_ligands"],
                      "ranking_eligible": result["ranking_eligible"]}))


if __name__ == "__main__":
    main()
