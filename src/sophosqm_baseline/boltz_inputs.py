"""Build the frozen 87-compound Boltz-2 input cohort."""
from __future__ import annotations

from pathlib import Path
import csv
import hashlib
import json

import yaml

from .cohort import resolve_structure_path


_AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "HID": "H", "HIE": "H", "HIP": "H", "CYX": "C",
}


def protein_sequence(pdb_path: str | Path, chain: str) -> str:
    residues: list[tuple[str, str]] = []
    with open(pdb_path, encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("ATOM  ") or line[21].strip() != chain:
                continue
            residue_name = line[17:20].upper()
            residue_number = line[22:27].strip()
            if line[12:16].strip().upper() == "CA" and residue_name in _AA3_TO_1:
                residues.append((residue_number, _AA3_TO_1[residue_name]))
    if not residues:
        raise ValueError(f"no standard protein residues found in chain {chain!r}")
    ordered: dict[str, str] = {}
    for number, code in residues:
        ordered.setdefault(number, code)
    return "".join(ordered.values())


def _sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_inputs(cohort_csv: str | Path, output_dir: str | Path,
                 msa_dir: str | Path, generate_msa: bool = False,
                 structure_root: str | Path | None = None) -> Path:
    output = Path(output_dir)
    msa_output = Path(msa_dir)
    output.mkdir(parents=True, exist_ok=True)
    msa_output.mkdir(parents=True, exist_ok=True)
    with open(cohort_csv, encoding="utf-8", newline="") as handle:
        cohort = list(csv.DictReader(handle))
    if len(cohort) != 87:
        raise ValueError(f"expected frozen 87-compound cohort, found {len(cohort)}")

    sequences: dict[str, str] = {}
    protein_paths: dict[str, str] = {}
    for row in cohort:
        target = row["target"]
        if target not in sequences:
            protein_path = row["protein_pdb"]
            if structure_root is not None:
                protein_path = str(resolve_structure_path(protein_path, structure_root))
            sequences[target] = protein_sequence(protein_path, row["protein_chain"])
            protein_paths[target] = protein_path

    if generate_msa:
        from boltz.main import compute_msa
        for target, sequence in sequences.items():
            expected = msa_output / f"{target}_1.csv"
            if expected.exists():
                continue
            compute_msa(
                data={f"{target}_1": sequence},
                target_id=target,
                msa_dir=msa_output,
                msa_server_url="https://api.colabfold.com",
                msa_pairing_strategy="greedy",
            )

    records = []
    for row in cohort:
        target = row["target"]
        ligand_name = row["ligand_name"].replace("/", "_")
        stem = f"{target}_{ligand_name}"
        schema = {
            "version": 1,
            "sequences": [
                {
                    "protein": {
                        "id": row["protein_chain"] or "A",
                        "sequence": sequences[target],
                        "msa": str((msa_output / f"{target}_1.csv").resolve()),
                    }
                },
                {"ligand": {"id": "B", "smiles": row["canonical_isomeric_smiles"]}},
            ],
            "properties": [{"affinity": {"binder": "B"}}],
        }
        input_path = output / f"{stem}.yaml"
        input_path.write_text(yaml.safe_dump(schema, sort_keys=False), encoding="utf-8")
        records.append({
            "target": target,
            "ligand_id": row["ligand_name"],
            "input": str(input_path),
            "input_sha256": _sha256(input_path),
            "protein_pdb": protein_paths[target],
            "protein_chain": row["protein_chain"],
            "protein_sequence_sha256": hashlib.sha256(sequences[target].encode()).hexdigest(),
            "msa": str((msa_output / f"{target}_1.csv").resolve()),
            "smiles": row["canonical_isomeric_smiles"],
            "experimental_dg_kcal_mol": float(row["experimental_dg_kcal_mol"]),
        })

    manifest_path = output / "manifest.json"
    manifest_path.write_text(json.dumps({
        "cohort": "Boltz2 Table 12 four-target FMO subset",
        "n_compounds": len(records),
        "counts": {target: sum(record["target"] == target for record in records)
                   for target in sequences},
        "msa_policy": "one shared mmseqs2/ColabFold MSA per target",
        "defaults": {
            "model": "boltz2",
            "seed": 0,
            "recycling_steps": 3,
            "sampling_steps": 200,
            "diffusion_samples": 1,
            "sampling_steps_affinity": 200,
            "diffusion_samples_affinity": 5,
        },
        "records": records,
    }, indent=2) + "\n", encoding="utf-8")
    return manifest_path
