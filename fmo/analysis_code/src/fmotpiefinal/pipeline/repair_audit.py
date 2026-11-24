"""Cohort-level audit and isolated candidate generation for geometry repairs."""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
from rdkit import Chem

from fmotpiefinal.core.bound_geometry import load_bound_ligand
from fmotpiefinal.core.repair_audit import (
    apply_experimental_nh2_repair,
    audit_local_strain_geometry,
)
from fmotpiefinal.pipeline.geometry_audit import element_from_pdb


def _digest(path: Path) -> str:
    import hashlib
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_pose_with_coordinates(source: Path, destination: Path, molecule,
                                 coordinates: np.ndarray) -> None:
    lines = source.read_text(encoding="utf-8").splitlines(keepends=True)
    ligand_line_indices = [index for index, line in enumerate(lines)
                           if line.startswith("HETATM") and line[17:20].strip() == "UNK"
                           and line[21].strip() == "X"]
    atoms = list(molecule.GetAtoms())
    if len(ligand_line_indices) != len(atoms) or coordinates.shape != (len(atoms), 3):
        raise ValueError("minimized PDB ligand atom inventory differs from SDF")
    for atom_index, line_index in enumerate(ligand_line_indices):
        line = lines[line_index]
        if element_from_pdb(line) != atoms[atom_index].GetSymbol().upper():
            raise ValueError("minimized PDB ligand atom order differs from SDF")
        x, y, z = coordinates[atom_index]
        lines[line_index] = line[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + line[54:]
    destination.parent.mkdir(parents=True, exist_ok=True)
    contents = "".join(lines)
    if destination.exists() and destination.read_text(encoding="utf-8") != contents:
        raise FileExistsError(f"refusing to replace a different repair candidate: {destination}")
    destination.write_text(contents, encoding="utf-8")


def audit_target(ligand_sdf: Path, minimized_root: Path, geometry_report: Path,
                 output: Path, *, local_strain_json: Path | None = None,
                 emit_experimental_candidates: bool = False,
                 target_name: str | None = None) -> dict:
    ligand_sdf = Path(ligand_sdf).resolve()
    minimized_root = Path(minimized_root).resolve()
    geometry_report = Path(geometry_report).resolve()
    output = Path(output).resolve()
    molecules = list(Chem.SDMolSupplier(str(ligand_sdf), removeHs=False))
    if not molecules or any(molecule is None for molecule in molecules):
        raise ValueError("repair audit requires a readable complete ligand SDF")
    names = [molecule.GetProp("_Name").strip() for molecule in molecules]
    if not all(names) or len(names) != len(set(names)):
        raise ValueError("repair audit requires unique nonempty ligand names")
    geometry = json.loads(geometry_report.read_text(encoding="utf-8"))
    geometry_cases = {row["ligand"]: row for row in geometry.get("cases", [])}
    if set(geometry_cases) != set(names):
        raise ValueError("geometry audit must cover exactly the full ligand cohort")
    strain_values = None
    if local_strain_json is not None:
        strain_values = json.loads(Path(local_strain_json).read_text(encoding="utf-8"))
        if not isinstance(strain_values, dict) or set(strain_values) != set(names):
            raise ValueError("local strain diagnostics must cover exactly the cohort")

    candidate_root = output.with_name(output.stem + "_experimental_nh2_repair_candidates")
    if emit_experimental_candidates:
        candidate_root.mkdir(parents=True, exist_ok=True)
    records = []
    for molecule, name in zip(molecules, names, strict=True):
        source_coordinates = np.asarray(molecule.GetConformer().GetPositions(), dtype=float)
        pose = minimized_root / name / "minimized_complex.pdb"
        if not pose.is_file():
            raise FileNotFoundError(f"missing minimized complex for {name}")
        bound_ligand = load_bound_ligand(ligand_sdf, name, pose)
        bound_coordinates = np.asarray([atom.xyz for atom in bound_ligand.atoms], dtype=float)
        if any(atom.GetSymbol().upper() != ligand_atom.element.upper()
               for atom, ligand_atom in zip(molecule.GetAtoms(), bound_ligand.atoms, strict=True)):
            raise ValueError(f"atom identity mismatch in minimized complex for {name}")
        pose_geometry = geometry_cases[name].get("minimized_pose")
        if pose_geometry is None:
            geometry_status = "INCOMPLETE"
        else:
            geometry_status = pose_geometry.get("quality_gate", "INCOMPLETE")
        strain = None if strain_values is None else float(strain_values[name])
        result = audit_local_strain_geometry(
            molecule, source_coordinates, bound_coordinates,
            geometry_status=geometry_status, local_strain_kcal_mol=strain)
        heavy_indices = [atom.GetIdx() for atom in molecule.GetAtoms()
                         if atom.GetAtomicNum() > 1]
        item = {"ligand": name, "atom_order_verified": True,
                "source_sdf_sha256": _digest(ligand_sdf),
                "minimized_pose_sha256": _digest(pose), "audit": result,
                "source_to_k10_heavy_rmsd_angstrom": float(np.sqrt(np.mean(np.sum(
                    (source_coordinates[heavy_indices] - bound_coordinates[heavy_indices]) ** 2,
                    axis=1))))}
        candidate = None
        if emit_experimental_candidates:
            destination = candidate_root / name / "minimized_complex.pdb"
            if result["repair_candidates"]:
                repaired, verification = apply_experimental_nh2_repair(
                    molecule, bound_coordinates, result)
                _write_pose_with_coordinates(pose, destination, molecule, repaired)
                item["experimental_repair"] = verification
                item["experimental_repair"]["candidate_pose_sha256"] = _digest(destination)
            else:
                destination.parent.mkdir(parents=True, exist_ok=True)
                if destination.exists() and destination.read_bytes() != pose.read_bytes():
                    raise FileExistsError(f"refusing to replace a different candidate pose: {destination}")
                shutil.copyfile(pose, destination)
            candidate = str(destination)
        if candidate:
            item["candidate_pose"] = candidate
        records.append(item)

    candidates = [row for row in records if row["audit"]["repair_candidates"]]
    report = {
        "schema": "fmotpiefinal-target-repair-audit-v1",
        "status": "COMPLETE" if len(records) == len(names) else "INCOMPLETE",
        "target": target_name or geometry.get("target"),
        "n_expected": len(names), "n_audited": len(records),
        "n_experimental_repair_candidates": len(candidates),
        "labels_used_for_trigger": False,
        "ranking_modified": False,
        "primary_score": "pure_tPIE",
        "local_strain_score_enabled": False,
        "candidate_generation_enabled": emit_experimental_candidates,
        "candidate_generation_status": (
            "CANDIDATES_ONLY_REQUIRES_SEPARATE_FULL_K10_AND_QM_RERUN"
            if emit_experimental_candidates and candidates else
            "NO_GEOMETRY_CANDIDATES" if emit_experimental_candidates else "NOT_REQUESTED"),
        "inputs": {"ligand_sdf": str(ligand_sdf), "ligand_sdf_sha256": _digest(ligand_sdf),
                   "minimized_root": str(minimized_root),
                   "geometry_report": str(geometry_report),
                   "geometry_report_sha256": _digest(geometry_report),
                   "local_strain_json": str(local_strain_json) if local_strain_json else None},
        "rules": {"conjugated_nh2_departure_degrees": 30.0,
                  "minimum_new_departure_degrees": 20.0,
                  "high_or_negative_strain_alone_repairs": False,
                  "target_or_affinity_labels_select_rules": False},
        "records": records,
    }
    if emit_experimental_candidates:
        report["candidate_root"] = str(candidate_root)
    output.parent.mkdir(parents=True, exist_ok=True)
    contents = json.dumps(report, indent=2) + "\n"
    if output.exists() and output.read_text(encoding="utf-8") != contents:
        raise FileExistsError(f"refusing to replace a different repair audit: {output}")
    output.write_text(contents, encoding="utf-8")
    return report
