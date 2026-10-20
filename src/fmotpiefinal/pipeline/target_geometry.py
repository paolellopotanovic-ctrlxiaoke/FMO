#!/usr/bin/env python
"""Apply label-blind target-level water audit to an intake with all bound poses."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from .structure_treatment import apply_plan
from .geometry_audit import audit as audit_geometry
from .water_policy import decide
from .protocol import freeze
from .target_intake import digest


def audit(intake_dir: Path, output: Path) -> dict:
    output = output.resolve()
    intake_file = intake_dir / "intake.json"
    intake = json.loads(intake_file.read_text())
    if intake.get("schema") != "fmo-new-target-intake-v1" or intake.get("status") not in (
            "STRUCTURE_PREFLIGHT_PASS", "STRUCTURE_REVIEW_REQUIRED"):
        raise ValueError("complete new-target 3D intake is required")
    cohort = Path(intake["cohort"])
    with cohort.open(newline="", encoding="utf-8") as handle:
        records = list(csv.DictReader(handle))
    if len(records) != intake["n_ligands"] or not records:
        raise ValueError("intake cohort has changed")
    protein = Path(records[0]["protein_pdb"])
    poses = Path(records[0]["ligand_sdf"])
    if (len({row["protein_pdb"] for row in records}) != 1 or
            len({row["ligand_sdf"] for row in records}) != 1 or
            digest(protein) != intake["protein_sha256"] or digest(poses) != intake["poses_sha256"]):
        raise ValueError("intake structures changed or multiple receptors/pose sets were supplied")
    geometry = audit_geometry(protein, poses)
    if (geometry.get("n_ligands") != len(records) or
            {row["ligand"] for row in geometry["cases"]} !=
            {row["ligand_name"] for row in records}):
        raise ValueError("target geometry audit does not cover all intake ligands")
    treatment = decide(geometry)
    if len(treatment["records"]) != len(records):
        raise ValueError("treatment omitted a ligand")
    output.mkdir(parents=True, exist_ok=True)
    (output / "source_geometry.json").write_text(json.dumps(geometry, indent=2) + "\n")
    plan_file = output / "target_water_plan.json"
    plan_file.write_text(json.dumps(treatment, indent=2) + "\n")
    if any("REVIEW_POSTMINIMIZATION_GEOMETRY" in item["actions"]
           for item in treatment["records"]):
        raise ValueError("cannot automatically override a failed post-minimization geometry gate")
    if treatment["target_water_residues_to_remove"]:
        prepared_protein = output / "receptor_target_water_treated.pdb"
        apply_plan(protein, plan_file, prepared_protein)
    else:
        prepared_protein = protein
    repaired_geometry = audit_geometry(prepared_protein, poses)
    (output / "treated_geometry.json").write_text(json.dumps(repaired_geometry, indent=2) + "\n")
    if any(item["n_severe_contacts"] for item in repaired_geometry["cases"]):
        severe = [contact for item in repaired_geometry["cases"]
                  for contact in item["severe_contacts"]]
        if any(contact["receptor_water"] for contact in severe):
            status = "REVIEW_UNPLANNED_SEVERE_WATER_CONTACTS"
            next_step = "update the label-blind target water audit before minimization"
        else:
            status = "K10_REPAIR_REQUIRED"
            next_step = "run uniform fixed-receptor k10 on all ligands; forbid QM if any severe contact persists"
    else:
        status = "STRUCTURE_TREATMENT_COMPLETE"
        next_step = "audit PDBFixer/PDB2PQR protonation, minimize all poses with fixed k10, then run full-cohort FMO"
    treated_cohort = output / "treated_structural_cohort.csv"
    with treated_cohort.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows({**row, "protein_pdb": str(prepared_protein)} for row in records)
    preflight = freeze(treated_cohort, output / "treated_structure_preflight.json")
    if not preflight["ready_for_full_target_run"]:
        status = "REVIEW_PREFLIGHT_AFTER_WATER_TREATMENT"
        next_step = "inspect source coordinates and ligand chemistry before minimization"
    return {"schema": "fmo-new-target-geometry-gate-v1", "target": intake["target"],
            "status": status, "ranking_eligible": False,
            "intake_sha256": digest(intake_file), "source_protein_sha256": digest(protein),
            "ligand_sdf_sha256": digest(poses), "treated_protein_sha256": digest(prepared_protein),
            "target_water_plan_sha256": digest(plan_file), "n_ligands": len(records),
            "n_target_waters_to_remove": len(treatment["target_water_residues_to_remove"]),
            "n_residual_severe_contacts": sum(item["n_severe_contacts"]
                                               for item in repaired_geometry["cases"]),
            "treated_cohort": str(treated_cohort), "next_step": next_step}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--intake-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.intake_dir, args.output)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "geometry_gate.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": report["status"], "n_ligands": report["n_ligands"],
                      "ranking_eligible": report["ranking_eligible"]}))


if __name__ == "__main__":
    main()
