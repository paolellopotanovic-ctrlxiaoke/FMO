"""Label-free warnings about transferability of a completed tPIE ranking."""
from __future__ import annotations

import json
import math
from pathlib import Path

from rdkit import Chem


def correlation(left: list[float], right: list[float]) -> float | None:
    if len(left) < 3:
        return None
    center_left = sum(left) / len(left)
    center_right = sum(right) / len(right)
    variance_left = sum((value - center_left) ** 2 for value in left)
    variance_right = sum((value - center_right) ** 2 for value in right)
    if variance_left == 0 or variance_right == 0:
        return None
    return sum((value - center_left) * (other - center_right)
               for value, other in zip(left, right, strict=True)) / math.sqrt(
                   variance_left * variance_right)


def audit(root: Path, ranking: dict, predictions: list[dict]) -> dict:
    """Report potential pitfalls without changing any ranking or selecting a protocol."""
    root = Path(root)
    report = {"schema": "fmotpiefinal-label-free-ranking-diagnostics-v1",
              "target": ranking["target"], "n": ranking["n"],
              "predictions_sha256": ranking["predictions_sha256"],
              "experimental_labels_used": False, "ranking_modified": False,
              "warnings": [], "per_ligand": []}
    intake_path = root / "intake/intake.json"
    if not intake_path.is_file():
        report["warnings"].append("INTAKE_NOT_AVAILABLE")
        return report
    intake = json.loads(intake_path.read_text(encoding="utf-8"))
    ligands = {entry["ligand_id"]: entry for entry in intake["ligands"]}
    ordered = sorted(ranking["ranked"], key=lambda entry: entry["ligand"])
    if len(ligands) != len(ordered) or set(ligands) != {entry["ligand"] for entry in ordered}:
        raise ValueError("intake ligand identities differ from QM ranking")
    from fmotpiefinal.scoring import digest
    report["intake_sha256"] = digest(intake_path)
    atom_counts = []
    charges = set()
    for entry in ordered:
        identity = entry["ligand"]
        molecule = Chem.MolFromSmiles(ligands[identity]["canonical_smiles"])
        if molecule is None:
            raise ValueError(f"invalid canonical SMILES: {identity}")
        atoms = molecule.GetNumHeavyAtoms()
        charge = Chem.GetFormalCharge(molecule)
        if charge != ligands[identity]["formal_charge"]:
            raise ValueError(f"intake charge inconsistent with SMILES: {identity}")
        atom_counts.append(atoms)
        charges.add(charge)
        report["per_ligand"].append({"ligand": identity, "heavy_atoms": atoms,
                                     "formal_charge": charge, "raw_score": entry["raw_score"]})
    report["heavy_atom_count_min_max"] = [min(atom_counts), max(atom_counts)]
    report["score_heavy_atom_pearson"] = correlation(
        [entry["raw_score"] for entry in ordered], atom_counts)
    if len(ordered) >= 8 and report["score_heavy_atom_pearson"] is not None and abs(
            report["score_heavy_atom_pearson"]) >= 0.7:
        report["warnings"].append("SCORE_TRACKS_HEAVY_ATOM_COUNT")
    if len(charges) > 1:
        report["warnings"].append("MIXED_LIGAND_FORMAL_CHARGES")
    fragment_counts = {entry.get("n_receptor_fragments") for entry in predictions}
    if len(fragment_counts) > 1:
        report["warnings"].append("NONUNIFORM_RECEPTOR_FRAGMENT_COUNTS")
    pose_audit_path = root / "qm/post_k10_geometry.json"
    if pose_audit_path.is_file():
        geometry = json.loads(pose_audit_path.read_text(encoding="utf-8"))
        report["post_k10_geometry_sha256"] = digest(pose_audit_path)
        cases = {entry["ligand"]: entry for entry in geometry["cases"]}
        if set(cases) != set(ligands):
            raise ValueError("post-k10 geometry cases differ from QM ranking")
        unresolved = []
        for entry in report["per_ligand"]:
            pose = cases[entry["ligand"]].get("minimized_pose")
            if not isinstance(pose, dict):
                unresolved.append(entry["ligand"])
            else:
                entry["minimized_pose_severe_contacts"] = pose.get("n_severe_contacts")
                if pose.get("n_severe_contacts") != 0:
                    unresolved.append(entry["ligand"])
        if unresolved:
            report["unresolved_minimized_poses"] = unresolved
            report["warnings"].append("MINIMIZED_POSE_NOT_CLEARED")
    else:
        report["warnings"].append("POST_K10_AUDIT_NOT_AVAILABLE")
    return report


def pocket_sensitivity(root: Path, ranking: dict) -> dict:
    """Keep previously frozen diagnostics untouched; audit pocket selection separately."""
    from fmotpiefinal.scoring import digest

    root = Path(root)
    pocket_path = root / "qm/pocket_selection_audit.json"
    result = {"schema": "fmotpiefinal-label-free-pocket-sensitivity-v1",
              "target": ranking["target"], "n": ranking["n"],
              "predictions_sha256": ranking["predictions_sha256"],
              "experimental_labels_used": False,
              "pocket_audit_available": pocket_path.is_file()}
    if pocket_path.is_file():
        audit = json.loads(pocket_path.read_text(encoding="utf-8"))
        names = audit.get("ligand_names")
        if (audit.get("n_ligands") != ranking["n"] or
                not isinstance(names, list) or len(names) != ranking["n"] or
                len(set(names)) != len(names) or set(names) != {
                    row["ligand"] for row in ranking["ranked"]}):
            raise ValueError("pocket membership audit does not match ranked cohort")
        result.update(pocket_audit_sha256=digest(pocket_path),
                      selected_residues=audit["n_residues_all_atom"],
                      heavy_selected_residues=audit["n_residues_heavy_atom"],
                      hydrogen_sensitive_residues=len(audit["membership_differences"]),
                      pocket_policy=audit["selection_production"])
    return result
