"""Exact construction of the Boltz-2 four-target FMO benchmark cohort."""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import csv
import hashlib
import json


TARGETS = {
    "CDK2": ("cdk2", "cdk2_ligands.sdf", "cdk2_protein.pdb", "A"),
    "JNK1": ("jnk1", "jnk1_manual_flips_ligands.sdf", "jnk1_manual_flips_protein.pdb", "A"),
    "P38": ("p38", "p38_ligands.sdf", "p38_protein.pdb", "C"),
    "TYK2": ("tyk2", "tyk2_ligands.sdf", "tyk2_protein.pdb", "A"),
}


@dataclass(frozen=True)
class CohortRecord:
    target: str
    ligand_name: str
    ligand_sdf: str
    protein_pdb: str
    protein_chain: str
    canonical_isomeric_smiles: str
    formal_charge: int
    experimental_dg_kcal_mol: float
    source_record_count: int
    representative_rule: str


def file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_structure_path(path_value: str | Path, structure_root: str | Path) -> Path:
    """Resolve a frozen source path against a locally downloaded structure root.

    Older frozen cohorts record the path used when the cohort was first built
    (for example, under ``vendor/...``). Those paths are provenance, not a
    requirement that users recreate that checkout location. Input runners map
    the source basename to ``structure_root`` when the recorded path is absent.
    """
    recorded = Path(path_value)
    root = Path(structure_root)
    candidates = [recorded, root / recorded.name]
    if not recorded.is_absolute():
        candidates.append(root / recorded)
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        f"structure input {recorded.name!r} not found; checked "
        + ", ".join(str(item) for item in candidates)
    )


def build_cohort(root: str | Path) -> list[CohortRecord]:
    """Collapse repeated JNK1 ligand/tautomer poses to unique compounds.

    The public Ross four-target set has 104 SDF records. Exactly 17 JNK1
    records are duplicate canonical isomeric structures with alternate flip
    poses, leaving 16+21+34+16=87 unique neutral compounds. For a duplicate,
    the non-flip pose is retained; this is deterministic and preserves the
    unique compound rather than selecting using benchmark labels.
    """
    from rdkit import Chem

    root = Path(root)
    records: list[CohortRecord] = []
    for target, (slug, ligand_file, protein_file, chain) in TARGETS.items():
        ligand_path = root / ligand_file
        protein_path = root / protein_file
        identities: dict[str, list[object]] = defaultdict(list)
        for molecule in Chem.SDMolSupplier(str(ligand_path), removeHs=False, sanitize=True):
            if molecule is None:
                raise ValueError(f"invalid SDF record in {ligand_path}")
            identities[Chem.MolToSmiles(Chem.RemoveHs(molecule), canonical=True, isomericSmiles=True)].append(molecule)
        for identity, molecules in sorted(identities.items()):
            if len(molecules) > 1 and target != "JNK1":
                raise ValueError(f"unexpected duplicate compound in {target}: {identity}")
            non_flips = [molecule for molecule in molecules if "flip" not in molecule.GetProp("_Name").lower()]
            if len(molecules) > 1 and len(non_flips) != 1:
                raise ValueError(f"cannot deterministically select JNK1 representative for {identity}")
            molecule = non_flips[0] if len(molecules) > 1 else molecules[0]
            records.append(CohortRecord(
                target=target,
                ligand_name=molecule.GetProp("_Name").strip(),
                ligand_sdf=str(ligand_path),
                protein_pdb=str(protein_path),
                protein_chain=chain,
                canonical_isomeric_smiles=identity,
                formal_charge=Chem.GetFormalCharge(molecule),
                experimental_dg_kcal_mol=float(molecule.GetProp("r_exp_dg")),
                source_record_count=len(molecules),
                representative_rule="unique_canonical_smiles_prefer_non_flip" if len(molecules) > 1 else "unique_canonical_smiles",
            ))
    if [sum(r.target == t for r in records) for t in TARGETS] != [16, 21, 34, 16] or len(records) != 87:
        raise AssertionError("four-target cohort does not reconstruct 87 compounds")
    if any(record.formal_charge != 0 for record in records):
        raise AssertionError("Boltz-2 describes this cohort as neutral compounds")
    return records


def write_cohort(path: str | Path, records: list[CohortRecord], inputs: dict[str, str]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "target", "ligand_name", "ligand_sdf", "protein_pdb", "protein_chain",
        "canonical_isomeric_smiles", "formal_charge", "experimental_dg_kcal_mol",
        "source_record_count", "representative_rule",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for record in records:
            writer.writerow({field: getattr(record, field) for field in fields})
    counts = {target: sum(record.target == target for record in records) for target in TARGETS}
    manifest = {
        "cohort": "Boltz2 Table 12 four-target FMO subset reconstructed from Chen/Ross structures",
        "n_compounds": len(records),
        "counts": counts,
        "construction": "unique RDKit canonical isomeric SMILES; JNK1 alternate flip records collapsed",
        "input_sha256": inputs,
    }
    path.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
