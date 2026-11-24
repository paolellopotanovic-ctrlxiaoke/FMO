#!/usr/bin/env python
"""Audit external FMO poses, including receptor-water and hydrogen clashes.

The original external preflight only queried receptor heavy atoms.  This audit
also checks all prepared receptor atoms, because a ligand can overlap a
crystallographic water hydrogen while appearing acceptable against protein
heavy atoms alone.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from rdkit import Chem


VDW_RADII = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
    "P": 1.80, "S": 1.80, "CL": 1.75, "BR": 1.85, "I": 1.98,
}


def element_from_pdb(line: str) -> str:
    element = line[76:78].strip().upper()
    if element:
        return element
    name = "".join(char for char in line[12:16] if char.isalpha()).upper()
    return name[:2] if name[:2] in VDW_RADII else name[:1]


def receptor_atoms(path: Path) -> list[dict]:
    atoms = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            atoms.append({
                "line_number": line_number,
                "serial": line[6:11].strip(),
                "name": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "chain": line[21].strip(),
                "resid": line[22:26].strip(),
                "element": element_from_pdb(line),
                "water": line[17:20].strip() in {"HOH", "WAT", "H2O"},
                "coordinates": [float(line[start:start + 8]) for start in (30, 38, 46)],
                "line": line,
            })
    if not atoms:
        raise ValueError(f"no receptor atoms found in {path}")
    return atoms


def ligand_atoms(path: Path) -> list[tuple[str, Chem.Mol]]:
    molecules = list(Chem.SDMolSupplier(str(path), removeHs=False))
    if not molecules or any(molecule is None for molecule in molecules):
        raise ValueError(f"unreadable ligand SDF: {path}")
    return [(molecule.GetProp("_Name").strip(), molecule) for molecule in molecules]


def pair_record(ligand_name: str, atom_index: int, ligand_atom, receptor_index: int,
                receptor_atom: dict, distance: float) -> dict:
    ligand_element = ligand_atom.GetSymbol().upper()
    receptor_element = receptor_atom["element"]
    radius_sum = VDW_RADII.get(ligand_element, 1.70) + VDW_RADII.get(receptor_element, 1.70)
    return {
        "ligand": ligand_name,
        "ligand_atom_index": atom_index,
        "ligand_element": ligand_element,
        "receptor_atom_serial": receptor_atom["serial"],
        "receptor_atom_name": receptor_atom["name"],
        "receptor_residue": f"{receptor_atom['resname']}{receptor_atom['resid']}",
        "receptor_chain": receptor_atom["chain"],
        "receptor_element": receptor_element,
        "receptor_water": receptor_atom["water"],
        "distance_angstrom": float(distance),
        "vdw_ratio": float(distance / radius_sum),
    }


def pair_record_from_values(ligand_name: str, atom_index: int, ligand_element: str,
                            receptor_atom: dict, distance: float) -> dict:
    radius_sum = VDW_RADII.get(ligand_element, 1.70) + VDW_RADII.get(
        receptor_atom["element"], 1.70)
    return {
        "ligand": ligand_name,
        "ligand_atom_index": atom_index,
        "ligand_element": ligand_element,
        "receptor_atom_serial": receptor_atom["serial"],
        "receptor_atom_name": receptor_atom["name"],
        "receptor_residue": f"{receptor_atom['resname']}{receptor_atom['resid']}",
        "receptor_chain": receptor_atom["chain"],
        "receptor_element": receptor_atom["element"],
        "receptor_water": receptor_atom["water"],
        "distance_angstrom": float(distance),
        "vdw_ratio": float(distance / radius_sum),
    }


def contact_is_severe(ligand_element: str, receptor_element: str, distance: float,
                      ratio: float, severe_ratio: float) -> bool:
    # A generic vdW ratio incorrectly labels ordinary hydrogen bonds as clashes.
    # Use hard distance gates for hydrogen-involving pairs and the vdW ratio for
    # heavy-atom pairs.
    if ligand_element == "H" and receptor_element == "H":
        return distance < 0.80
    if "H" in {ligand_element, receptor_element}:
        return distance < 1.15
    return ratio <= severe_ratio


def minimized_atoms(path: Path, ligand_resname: str = "UNK",
                    ligand_chain: str = "X") -> tuple[list[dict], list[dict]]:
    """Split a minimized complex PDB into receptor and ligand atoms.

    The external minimizer writes the ligand as ``UNK`` on chain ``X``.  The
    parser intentionally keeps every other ATOM/HETATM record as receptor,
    including water and explicit receptor hydrogens.
    """
    receptor = []
    ligand = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            atom = {
                "line_number": line_number,
                "serial": line[6:11].strip(),
                "name": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "chain": line[21].strip(),
                "resid": line[22:26].strip(),
                "element": element_from_pdb(line),
                "water": line[17:20].strip() in {"HOH", "WAT", "H2O"},
                "coordinates": [float(line[start:start + 8]) for start in (30, 38, 46)],
            }
            if atom["resname"] == ligand_resname and atom["chain"] == ligand_chain:
                ligand.append(atom)
            else:
                receptor.append(atom)
    if not receptor or not ligand:
        raise ValueError(f"could not split minimized complex into receptor/ligand: {path}")
    return receptor, ligand


def audit_minimized_complex(path: Path, severe_ratio: float = 0.60,
                            warning_ratio: float = 0.75) -> dict:
    receptor, ligand = minimized_atoms(path)
    receptor_coordinates = np.asarray([atom["coordinates"] for atom in receptor], dtype=float)
    tree = cKDTree(receptor_coordinates)
    all_pairs = []
    for atom_index, atom in enumerate(ligand):
        distances, indices = tree.query(atom["coordinates"], k=min(8, len(receptor)))
        for distance, receptor_index in zip(np.atleast_1d(distances), np.atleast_1d(indices)):
            item = pair_record_from_values(path.parent.name, atom_index, atom["element"],
                                           receptor[int(receptor_index)], float(distance))
            warning = (item["vdw_ratio"] <= warning_ratio or
                       (item["ligand_element"] == "H" or item["receptor_element"] == "H") and
                       item["distance_angstrom"] < 1.60)
            if warning:
                all_pairs.append(item)
    severe = [item for item in all_pairs if contact_is_severe(
        item["ligand_element"], item["receptor_element"], item["distance_angstrom"],
        item["vdw_ratio"], severe_ratio)]
    return {
        "path": str(path),
        "n_receptor_atoms": len(receptor),
        "n_ligand_atoms": len(ligand),
        "n_warning_contacts": len(all_pairs),
        "n_severe_contacts": len(severe),
        "severe_contacts": sorted(severe, key=lambda item: item["vdw_ratio"]),
        "warning_contacts": sorted(all_pairs, key=lambda item: item["vdw_ratio"])[:20],
        "quality_gate": "REPAIR_REQUIRED" if severe else ("REVIEW" if all_pairs else "PASS"),
    }


def audit(protein: Path, ligands: Path, minimized_root: Path | None = None,
          severe_ratio: float = 0.60, warning_ratio: float = 0.75) -> dict:
    receptor = receptor_atoms(protein)
    receptor_coordinates = np.asarray([atom["coordinates"] for atom in receptor], dtype=float)
    tree = cKDTree(receptor_coordinates)
    heavy_indices = [index for index, atom in enumerate(receptor)
                     if atom["element"] not in {"H", "D"}]
    heavy_tree = cKDTree(receptor_coordinates[heavy_indices])
    records = []
    water_ids = defaultdict(set)
    for name, molecule in ligand_atoms(ligands):
        if molecule.GetNumConformers() != 1 or not molecule.GetConformer().Is3D():
            raise ValueError(f"ligand has no single 3D conformer: {name}")
        positions = np.asarray(molecule.GetConformer().GetPositions(), dtype=float)
        heavy_indices = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1]
        all_pairs = []
        for atom_index, atom in enumerate(molecule.GetAtoms()):
            distances, indices = tree.query(positions[atom_index], k=min(8, len(receptor)))
            distances = np.atleast_1d(distances)
            indices = np.atleast_1d(indices)
            for distance, receptor_index in zip(distances, indices):
                item = pair_record(name, atom_index, atom, int(receptor_index),
                                   receptor[int(receptor_index)], float(distance))
                warning = (item["vdw_ratio"] <= warning_ratio or
                           (item["ligand_element"] == "H" or item["receptor_element"] == "H") and
                           item["distance_angstrom"] < 1.60)
                if warning:
                    all_pairs.append(item)
                    if contact_is_severe(item["ligand_element"], item["receptor_element"],
                                         item["distance_angstrom"], item["vdw_ratio"],
                                         severe_ratio) and item["receptor_water"]:
                        water_ids[name].add((receptor[int(receptor_index)]["chain"],
                                             receptor[int(receptor_index)]["resid"],
                                             receptor[int(receptor_index)]["resname"]))
        severe = [item for item in all_pairs if contact_is_severe(
            item["ligand_element"], item["receptor_element"], item["distance_angstrom"],
            item["vdw_ratio"], severe_ratio)]
        warnings = all_pairs
        nearest_heavy = min(
            (heavy_tree.query(positions[index])[0] for index in
             [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1]),
            default=math.nan)
        minimization = None
        minimized_pose = None
        if minimized_root is not None:
            path = minimized_root / name / "minimization.json"
            if path.is_file():
                minimization = json.loads(path.read_text(encoding="utf-8"))
            pose_path = minimized_root / name / "minimized_complex.pdb"
            if pose_path.is_file():
                minimized_pose = audit_minimized_complex(
                    pose_path, severe_ratio=severe_ratio, warning_ratio=warning_ratio)
        records.append({
            "ligand": name,
            "formal_charge": int(Chem.GetFormalCharge(molecule)),
            "explicit_hydrogens": bool(any(atom.GetAtomicNum() == 1 for atom in molecule.GetAtoms())),
            "n_atoms": molecule.GetNumAtoms(),
            "n_heavy_atoms": molecule.GetNumHeavyAtoms(),
            "nearest_receptor_heavy_angstrom": float(nearest_heavy),
            "n_warning_contacts": len(warnings),
            "n_severe_contacts": len(severe),
            "severe_contacts": sorted(severe, key=lambda item: item["vdw_ratio"]),
            "warning_contacts": sorted(warnings, key=lambda item: item["vdw_ratio"])[:20],
            "minimization": minimization,
            "minimized_pose": minimized_pose,
            "quality_gate": "REPAIR_REQUIRED" if severe else ("REVIEW" if warnings else "PASS"),
        })
    water_residues = sorted({water for waters in water_ids.values() for water in waters})
    minimized_cases = [item["minimized_pose"] for item in records
                       if item["minimized_pose"] is not None]
    minimized_gate = None
    if minimized_cases:
        if len(minimized_cases) != len(records):
            minimized_gate = "INCOMPLETE"
        else:
            minimized_gate = ("REPAIR_REQUIRED" if any(item["quality_gate"] == "REPAIR_REQUIRED"
                                                       for item in minimized_cases)
                              else ("REVIEW" if any(item["quality_gate"] == "REVIEW"
                                                    for item in minimized_cases) else "PASS"))
    return {
        "protein": str(protein),
        "ligands": str(ligands),
        "n_ligands": len(records),
        "n_receptor_atoms": len(receptor),
        "n_receptor_waters": len({(a["chain"], a["resid"]) for a in receptor if a["water"]}),
        "thresholds": {"severe_vdw_ratio": severe_ratio, "warning_vdw_ratio": warning_ratio},
        "quality_gate": "REPAIR_REQUIRED" if any(r["quality_gate"] == "REPAIR_REQUIRED" for r in records)
        else ("REVIEW" if any(r["quality_gate"] == "REVIEW" for r in records) else "PASS"),
        "minimized_pose_gate": minimized_gate,
        "n_minimized_pose_audits": len(minimized_cases),
        "effective_quality_gate": (
            "INCOMPLETE" if len(minimized_cases) != len(records)
            else ("REPAIR_REQUIRED" if any(item["quality_gate"] == "REPAIR_REQUIRED"
                                           for item in minimized_cases)
                  else "PASS_AFTER_MINIMIZATION")
        ) if minimized_cases else "INCOMPLETE",
        "water_residues_with_severe_contacts": [
            {"chain": chain, "resid": resid, "resname": resname}
            for chain, resid, resname in water_residues
        ],
        "cases": records,
    }


def repair_waters(protein: Path, output: Path, water_residues: list[dict]) -> dict:
    blocked = {(item["chain"], item["resid"], item["resname"]) for item in water_residues}
    kept = removed = 0
    output.parent.mkdir(parents=True, exist_ok=True)
    with protein.open(encoding="utf-8") as source, output.open("w", encoding="utf-8") as destination:
        for line in source:
            if line.startswith(("ATOM  ", "HETATM")) and line[17:20].strip() in {"HOH", "WAT", "H2O"}:
                key = (line[21].strip(), line[22:26].strip(), line[17:20].strip())
                if key in blocked:
                    removed += 1
                    continue
            destination.write(line)
            if line.startswith(("ATOM  ", "HETATM")):
                kept += 1
    return {"input": str(protein), "output": str(output),
            "removed_water_atoms": removed, "retained_atom_records": kept,
            "removed_water_residues": water_residues}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protein", type=Path, required=True)
    parser.add_argument("--ligands", type=Path, required=True)
    parser.add_argument("--minimized-root", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repair-protein", type=Path)
    args = parser.parse_args()
    report = audit(args.protein, args.ligands, args.minimized_root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if args.repair_protein:
        repair = repair_waters(args.protein, args.repair_protein,
                               report["water_residues_with_severe_contacts"])
        report["repair"] = repair
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: report[key] for key in (
        "n_ligands", "n_receptor_atoms", "quality_gate",
        "minimized_pose_gate", "effective_quality_gate",
        "water_residues_with_severe_contacts")}, indent=2))
    if (report["quality_gate"] == "REPAIR_REQUIRED"
            and report.get("effective_quality_gate") != "PASS_AFTER_MINIMIZATION"
            and not args.repair_protein):
        raise SystemExit("structure quality gate requires repair")


if __name__ == "__main__":
    main()
