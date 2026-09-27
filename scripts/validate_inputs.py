#!/usr/bin/env python
"""Validate that all frozen cohort structures are locally runnable."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rdkit import Chem

from sophosqm_baseline.cohort import resolve_structure_path


EXPECTED_COUNTS = {"CDK2": 16, "JNK1": 21, "P38": 34, "TYK2": 16}
EXPECTED_COHORT_SHA256 = "19edfb8e2ccf3ab1029d848e4062d70cb35532d2b4fee63674f160b2a292997b"
EXPECTED_MSA_SHA256 = {
    "CDK2_1.csv": "980f0da21b711f9ddb2a5865bd428d25bdb7cafc6373a78ede893bc8086238f4",
    "JNK1_1.csv": "b0eecdfc9b96e9159e4278893cdcedf0eff7382154e8ecfae3410ef22c247898",
    "P38_1.csv": "b950710b09fc8fc03af213b9980201bf0a83fcd0bf21cfd31adca04757cb4f77",
    "TYK2_1.csv": "c7e2f9820cd37c2d90ebda67f02bccf2cf2ec87216ce51f98b19bc09aa318788",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", default="data/fmo87_cohort.csv")
    parser.add_argument("--structure-root", default="data/fmo87_structures")
    parser.add_argument("--msa-dir", default="data/boltz_msa")
    parser.add_argument("--output", default="runs/fmo87/input_validation.json")
    args = parser.parse_args()

    cohort_path = Path(args.cohort)
    structure_root = Path(args.structure_root)
    with cohort_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    counts = {target: sum(row["target"] == target for row in rows) for target in EXPECTED_COUNTS}
    if len(rows) != 87 or counts != EXPECTED_COUNTS:
        raise ValueError(f"invalid cohort counts: {counts}, n={len(rows)}")
    if _sha256(cohort_path) != EXPECTED_COHORT_SHA256:
        raise ValueError("frozen cohort SHA-256 mismatch")

    ligand_files: dict[str, tuple[Path, list[object]]] = {}
    protein_files: dict[str, Path] = {}
    for row in rows:
        ligand_value = row["ligand_sdf"]
        protein_value = row["protein_pdb"]
        if ligand_value not in ligand_files:
            path = resolve_structure_path(ligand_value, structure_root)
            molecules = [item for item in Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)]
            if any(item is None for item in molecules):
                raise ValueError(f"invalid SDF record in {path}")
            ligand_files[ligand_value] = (path, molecules)
        if protein_value not in protein_files:
            protein_files[protein_value] = resolve_structure_path(protein_value, structure_root)

        path, molecules = ligand_files[ligand_value]
        matched = [item for item in molecules if item.GetProp("_Name").strip() == row["ligand_name"]]
        if len(matched) != 1:
            raise ValueError(f"expected one SDF record named {row['ligand_name']!r} in {path}")
        molecule = matched[0]
        identity = Chem.MolToSmiles(
            Chem.RemoveHs(molecule), canonical=True, isomericSmiles=True
        )
        if identity != row["canonical_isomeric_smiles"]:
            raise ValueError(f"identity mismatch for {row['target']}_{row['ligand_name']}")
        if Chem.GetFormalCharge(molecule) != int(row["formal_charge"]):
            raise ValueError(f"formal charge mismatch for {row['target']}_{row['ligand_name']}")

    manifest = json.loads(cohort_path.with_suffix(".manifest.json").read_text(encoding="utf-8"))
    source_hashes = {}
    for recorded_name, expected_hash in manifest["input_sha256"].items():
        path = resolve_structure_path(recorded_name, structure_root)
        actual_hash = _sha256(path)
        if actual_hash != expected_hash:
            raise ValueError(f"source hash mismatch for {path}")
        source_hashes[path.name] = actual_hash

    output = {
        "status": "PASS",
        "cohort": str(cohort_path),
        "cohort_sha256": EXPECTED_COHORT_SHA256,
        "structure_root": str(structure_root),
        "n_compounds": len(rows),
        "counts": counts,
        "n_ligand_files": len(ligand_files),
        "n_protein_files": len(protein_files),
        "all_records_resolved": True,
        "all_identities_and_charges_match": True,
        "source_hashes_match_manifest": True,
        "source_sha256": source_hashes,
        "msa_sha256_match": True,
    }

    msa_hashes = {}
    for name, expected_hash in EXPECTED_MSA_SHA256.items():
        path = Path(args.msa_dir) / name
        actual_hash = _sha256(path)
        if actual_hash != expected_hash:
            raise ValueError(f"MSA hash mismatch for {path}")
        msa_hashes[name] = actual_hash
    output["msa_sha256"] = msa_hashes
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
