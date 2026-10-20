#!/usr/bin/env python
"""Create a label-free, fail-closed structure treatment plan from an audit."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def _water_key(item: dict) -> dict:
    residue = str(item.get("receptor_residue", ""))
    resname = str(item.get("receptor_resname", "HOH"))
    resid = str(item.get("receptor_resid", ""))
    if not resid and residue:
        prefix = residue.rstrip("0123456789")
        resid = residue[len(prefix):]
        if prefix:
            resname = prefix
    return {
        "chain": str(item.get("chain", item.get("receptor_chain", ""))),
        "resid": resid or str(item.get("resid", "")),
        "resname": resname,
    }


def decide(audit: dict) -> dict:
    records = audit.get("records", audit.get("cases"))
    if not isinstance(records, list) or not records:
        raise ValueError("structure audit must contain non-empty records")

    rows = []
    target_water_residues = {}
    for record in records:
        ligand = str(record.get("ligand", ""))
        if not ligand:
            raise ValueError("structure audit contains a record without ligand ID")
        severe = list(record.get("severe_contacts") or [])
        severe_water = [_water_key(contact) for contact in severe
                        if contact.get("receptor_water") is True]
        for water in severe_water:
            target_water_residues.setdefault(json.dumps(water, sort_keys=True), water)
        severe_nonwater = [contact for contact in severe
                           if contact.get("receptor_water") is not True]
        actions = []
        if not record.get("explicit_hydrogens", False):
            actions.append("REGENERATE_LIGAND_HYDROGENS")
        if severe_water:
            actions.append("REMOVE_ONLY_SEVERE_CONTACT_WATERS")
        if severe_nonwater:
            actions.append("REVIEW_NONWATER_SEVERE_CONTACT")
        minimized_pose = record.get("minimized_pose") or {}
        if minimized_pose.get("quality_gate") == "REPAIR_REQUIRED" or \
                minimized_pose.get("n_severe_contacts", 0) > 0:
            actions.append("REVIEW_POSTMINIMIZATION_GEOMETRY")
        if not actions and record.get("n_warning_contacts", 0):
            actions.append("RETAIN_AND_RECORD_WARNING_CONTACTS")
        if not actions:
            actions.append("RETAIN_STRUCTURE")
        if "REVIEW_NONWATER_SEVERE_CONTACT" in actions or \
                "REVIEW_POSTMINIMIZATION_GEOMETRY" in actions:
            status = "REVIEW"
        elif any(action.startswith("REGENERATE") or action.startswith("REMOVE")
                 for action in actions):
            status = "REPAIR_REQUIRED"
        else:
            status = "PASS"
        rows.append({
            "ligand": ligand,
            "status": status,
            "actions": actions,
            "remove_water_residues": sorted(severe_water,
                                             key=lambda item: (item["chain"], item["resid"])),
            "n_severe_contacts": len(severe),
            "n_severe_water_contacts": len(severe_water),
            "n_severe_nonwater_contacts": len(severe_nonwater),
        })

    statuses = {row["status"] for row in rows}
    overall = "REVIEW" if "REVIEW" in statuses else (
        "REPAIR_REQUIRED" if "REPAIR_REQUIRED" in statuses else "PASS")
    union_waters = sorted(target_water_residues.values(),
                          key=lambda item: (item["chain"], item["resid"], item["resname"]))
    if union_waters:
        target_action = "REMOVE_TARGET_UNION_OF_SEVERE_CONTACT_WATERS"
        target_status = "REPAIR_REQUIRED"
    else:
        target_action = "RETAIN_TARGET_WATERS"
        target_status = "PASS" if overall == "PASS" else overall
    return {
        "schema": "fmo-structure-treatment-plan-v1",
        "label_free": True,
        "overall_status": overall,
        "decision_rule": (
            "remove only receptor waters implicated in severe ligand contacts; "
            "regenerate ligand hydrogens without moving heavy atoms; "
            "never automatically move nonwater atoms or use affinity labels"
        ),
        "water_policy": "retain all waters unless a severe ligand-water contact is audited",
        "target_policy": "apply one shared receptor-water decision to every ligand in the target",
        "target_status": target_status,
        "target_actions": [target_action],
        "target_water_residues_to_remove": union_waters,
        "hydrogen_policy": "regenerate explicit ligand hydrogens while retaining source heavy-atom coordinates",
        "nonwater_severe_policy": "fail closed to REVIEW; no automatic heavy-atom repair",
        "records": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = decide(json.loads(args.audit.read_text(encoding="utf-8")))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"overall_status": report["overall_status"],
                      "n_records": len(report["records"]),
                      "output": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
