"""Label-free preparation diagnostics; never modify coordinates or ranking."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re

import numpy as np
from rdkit import Chem

from .core.bound_geometry import load_bound_ligand
from .core.structure import read_pdb, residue_key
from .pipeline.geometry_audit import audit_minimized_complex
from .scoring import digest, score_run, selected_predictions


def parse_propka_summary(text: str) -> dict[str, float]:
    if text.count("SUMMARY OF THIS PREDICTION") != 1:
        raise ValueError("expected one unambiguous PROPKA summary")
    table = text.split("SUMMARY OF THIS PREDICTION", 1)[1]
    table = table.split("Free energy of", 1)[0]
    records = {}
    for line in table.splitlines():
        match = re.fullmatch(r"\s*(ASP|GLU|HIS|LYS|ARG|CYS|TYR)\s+(-?\d+)\s+(\S)\s+"
                             r"([-+]?\d+(?:\.\d+)?)\s+[-+]?\d+(?:\.\d+)?\s*", line)
        if match is None:
            continue
        residue, number, chain, value = match.groups()
        key = f"{chain}:{residue}:{number}"
        if key in records:
            raise ValueError("duplicate PROPKA residue in summary")
        records[key] = float(value)
    if not records:
        raise ValueError("PROPKA summary contains no recognized residue pKa")
    return records


def audit_preparation(root: Path) -> dict:
    root = Path(root).resolve()
    protocol_path = root / "frozen_protocol.json"
    protocol = json.loads(protocol_path.read_text())
    names, target = protocol["ligand_names"], protocol["target_name"]
    if not names or len(names) != len(set(names)):
        raise ValueError("preparation audit requires unique complete ligand names")
    sdf = root / "ligands_no_properties.sdf"
    receptor_path = root / f"{target.lower()}_prepared.pdb"
    poses = Path(protocol.get("poses_root", str(root)))
    molecules = list(Chem.SDMolSupplier(str(sdf), removeHs=False))
    if not molecules or any(molecule is None for molecule in molecules):
        raise ValueError("unreadable prepared SDF")
    by_name = {molecule.GetProp("_Name").strip(): molecule for molecule in molecules}
    if len(molecules) != len(names) or set(by_name) != set(names):
        raise ValueError("SDF differs from complete preparation protocol")
    protein = read_pdb(receptor_path)
    report = {"schema": "fmotpiefinal-preparation-sensitivity-v1", "target": target, "n": len(names),
              "experimental_labels_used": False, "ranking_modified": False,
              "coordinates_modified": False, "warnings": [], "per_ligand": [],
              "protocol_sha256": digest(protocol_path), "sdf_sha256": digest(sdf),
              "prepared_receptor_sha256": digest(receptor_path)}
    all_coordinates = []
    for name in names:
        pose = poses / "minimized" / name / "minimized_complex.pdb"
        ligand = load_bound_ligand(sdf, name, pose)
        molecule = by_name[name]
        heavy = np.asarray([atom.GetAtomicNum() > 1 for atom in molecule.GetAtoms()])
        source = np.asarray(molecule.GetConformer().GetPositions())
        current = np.asarray([atom.xyz for atom in ligand.atoms])
        if not heavy.any() or not np.isfinite(current).all() or current.shape != source.shape:
            raise ValueError("invalid or incomplete prepared coordinates")
        displacement = np.linalg.norm(current[heavy] - source[heavy], axis=1)
        quality = audit_minimized_complex(pose)
        report["per_ligand"].append({"ligand": name,
            "formal_charge": Chem.GetFormalCharge(molecule), "n_heavy_atoms": int(heavy.sum()),
            "source_name_mentions_docked": "docked" in name.lower(),
            "source_to_minimized_heavy_rmsd_no_alignment_angstrom": float(np.sqrt(np.mean(displacement ** 2))),
            "source_to_minimized_heavy_max_displacement_angstrom": float(displacement.max()),
            "n_severe_contacts": quality["n_severe_contacts"],
            "n_warning_contacts": quality["n_warning_contacts"], "pose_sha256": digest(pose)})
        all_coordinates.extend(current[heavy])
    if any(case["source_name_mentions_docked"] for case in report["per_ligand"]):
        report["warnings"].append("SOURCE_NAMES_INCLUDE_DOCKED_POSES_VERIFY_POSE_PROVENANCE")
    if any(case["n_severe_contacts"] for case in report["per_ligand"]):
        report["warnings"].append("MINIMIZED_POSE_HAS_SEVERE_CONTACTS")
    metadata = root / f"{target.lower()}_prepared.preparation.json"
    if not metadata.exists():
        report["warnings"].append("PREPARATION_PH_NOT_AVAILABLE")
        return report
    preparation = json.loads(metadata.read_text())
    ph = float(preparation["ph"])
    if not math.isfinite(ph):
        raise ValueError("nonfinite preparation pH")
    report["preparation_ph"] = ph
    report["preparation_metadata_sha256"] = digest(metadata)
    log = root / f"{target.lower()}_prepared.log"
    if not log.exists():
        report["warnings"].append("PROPKA_SUMMARY_NOT_AVAILABLE")
        return report
    text = log.read_text(errors="replace")
    if "SUMMARY OF THIS PREDICTION" not in text:
        report["warnings"].append("PROPKA_SUMMARY_NOT_AVAILABLE")
        return report
    pka = parse_propka_summary(text)
    report["propka_log_sha256"] = digest(log)
    atoms_by_residue = {}
    for atom in protein.atoms:
        atoms_by_residue.setdefault(residue_key(atom), []).append(atom)
    ligand_coordinates = np.asarray(all_coordinates)
    nearby = []
    for key, value in pka.items():
        atoms = atoms_by_residue.get(key)
        if atoms is None:
            raise ValueError(f"PROPKA residue identity absent from prepared receptor: {key}")
        heavy = [atom for atom in atoms if atom.element not in {"H", "D"}]
        distance = min(float(np.min(np.linalg.norm(ligand_coordinates - atom.xyz, axis=1))) for atom in heavy)
        if distance > 5 or abs(value - ph) > 1:
            continue
        nearby.append({"residue": key, "predicted_pka": value, "distance_to_preparation_ph": abs(value - ph),
                       "minimum_ligand_heavy_distance_angstrom": distance,
                       "prepared_hydrogen_names": [atom.atom_name for atom in atoms if atom.element == "H"]})
    report["near_ligand_pka_within_one_ph_unit"] = nearby
    report["pka_review_rule"] = "residue within 5A of any ligand heavy atom AND abs(pKa-preparation_pH)<=1; warning only"
    if nearby:
        report["warnings"].append("NEAR_LIGAND_PROTONATION_SENSITIVE_AT_PREPARATION_PH")
    return report


def audit_selected_preparation(root: Path, ranking: dict) -> dict:
    root = Path(root).resolve()
    source = Path(ranking["selected_predictions"]).resolve()
    status, validated_source, predictions = selected_predictions(root)
    by_name = {row["ligand"]: row for row in predictions}
    ranked = ranking.get("ranked", [])
    if (not source.is_relative_to(root) or source != validated_source.resolve() or
            digest(source) != ranking["predictions_sha256"] or
            ranking["status_sha256"] != digest(root / "adaptive_status.json") or
            ranking["target"] != status["target"] or ranking["n"] != len(predictions) or
            len(ranked) != len(predictions) or {row["ligand"] for row in ranked} != set(by_name) or
            any(not math.isclose(row["raw_score"], -float(by_name[row["ligand"]]["tpie"]), abs_tol=1e-8)
                for row in ranked)):
        raise ValueError("selected preparation is outside validated run or has changed")
    protocol_path = source.parent / "frozen_protocol.json"
    if not protocol_path.exists():
        return {"schema": "fmotpiefinal-preparation-sensitivity-v1", "target": ranking["target"],
                "n": ranking["n"], "experimental_labels_used": False, "ranking_modified": False,
                "coordinates_modified": False, "warnings": ["SELECTED_PREPARATION_PROTOCOL_NOT_AVAILABLE"],
                "predictions_sha256": ranking["predictions_sha256"]}
    report = audit_preparation(source.parent)
    if report["target"] != ranking["target"] or report["n"] != ranking["n"]:
        raise ValueError("selected preparation differs from validated complete ranking")
    report["predictions_sha256"] = ranking["predictions_sha256"]
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--qm-dir", type=Path)
    source.add_argument("--run-dir", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.run_dir is not None:
        ranking_path = args.run_dir / "ranking/ranking.json"
        ranking = json.loads(ranking_path.read_text()) if ranking_path.exists() else score_run(args.run_dir)
        report = audit_selected_preparation(args.run_dir, ranking)
        destination = args.output or args.run_dir / "ranking/preparation_sensitivity.json"
    else:
        if args.output is None:
            parser.error("--qm-dir requires an explicit --output")
        report = audit_preparation(args.qm_dir)
        destination = args.output
    text = json.dumps(report, indent=2) + "\n"
    if destination.exists() and json.loads(destination.read_text()) != report:
        raise ValueError("existing preparation audit differs; use a fresh output")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(text)
    print(json.dumps({key: report[key] for key in ("target", "n", "warnings")}))


if __name__ == "__main__":
    main()
