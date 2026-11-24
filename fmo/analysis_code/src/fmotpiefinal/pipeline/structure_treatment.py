#!/usr/bin/env python
"""Apply a frozen label-free structure treatment plan to a receptor PDB."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def apply_plan(protein: Path, plan: Path, output: Path) -> dict:
    report = json.loads(plan.read_text(encoding="utf-8"))
    if report.get("label_free") is not True:
        raise ValueError("treatment plan must be label-free")
    residues = {(str(item.get("chain", "")), str(item["resid"]), str(item["resname"]))
                for item in report.get("target_water_residues_to_remove", [])}
    output.parent.mkdir(parents=True, exist_ok=True)
    removed_lines = 0
    removed_residues = set()
    with protein.open(encoding="utf-8") as source, output.open("w", encoding="utf-8") as target:
        for line in source:
            if line.startswith(("ATOM  ", "HETATM")):
                resname = line[17:20].strip()
                chain = line[21].strip()
                resid = line[22:26].strip()
                key = (chain, resid, resname)
                matches = {candidate for candidate in residues
                           if candidate[1:] == (resid, resname) and
                           (not candidate[0] or candidate[0] == chain)}
                if resname in {"HOH", "WAT", "H2O"} and matches:
                    removed_lines += 1
                    removed_residues.update(matches)
                    continue
            target.write(line)
    missing = sorted({candidate for candidate in residues
                      if not any(found[1:] == candidate[1:] for found in removed_residues)})
    if missing:
        raise ValueError(f"planned water residues not found in receptor: {missing}")
    result = {
        "schema": "fmo-structure-treatment-application-v1",
        "label_free": True,
        "input": str(protein),
        "output": str(output),
        "input_sha256": digest(protein),
        "output_sha256": digest(output),
        "removed_water_residues": [
            {"chain": chain, "resid": resid, "resname": resname}
            for chain, resid, resname in sorted(removed_residues)
        ],
        "removed_atom_lines": removed_lines,
        "policy": report.get("decision_rule"),
    }
    (output.parent / "structure_treatment_application.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protein", type=Path, required=True)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(apply_plan(args.protein, args.plan, args.output), indent=2))


if __name__ == "__main__":
    main()
