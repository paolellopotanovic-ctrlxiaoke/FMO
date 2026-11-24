"""Audit QuanPol endpoint pilot logs and write a non-FEP report."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from run_matched_openfe import write_json


PATTERNS = {
    "bond_stretching_kcal_mol": r"BOND STRETCHING\s+ENERGY =\s+([-+0-9.Ee]+)",
    "bond_angle_bending_kcal_mol": r"BOND ANGLE BENDING\s+ENERGY =\s+([-+0-9.Ee]+)",
    "dihedral_rotation_kcal_mol": r"DIHEDRAL ROTATION\s+ENERGY =\s+([-+0-9.Ee]+)",
    "lj_repulsion_kcal_mol": r"LJ REPULSION\s+ENERGY =\s+([-+0-9.Ee]+)",
    "lj_dispersion_kcal_mol": r"LJ DISPERSION\s+ENERGY =\s+([-+0-9.Ee]+)",
    "charge_kcal_mol": r"CHARGE\s+ENERGY =\s+([-+0-9.Ee]+)",
    "qm_center_kcal_mol": r"QM CENTER\s+ENERGY =\s+([-+0-9.Ee]+)\s+KCAL/MOL",
    "potential_kcal_mol": r"POTENTIAL\s+ENERGY =\s+([-+0-9.Ee]+)",
    "total_kcal_mol": r"TOTAL\s+ENERGY =\s+([-+0-9.Ee]+)",
    "total_charge": r"TOTAL QM/MM CHARGE =\s+([-+0-9.Ee]+)",
}


def parse(path: Path):
    text = path.read_text(errors="replace")
    row = {key: float(match.group(1)) for key, pattern in PATTERNS.items()
           if (match := re.search(pattern, text))}
    for key in ("lj_repulsion_kcal_mol", "lj_dispersion_kcal_mol", "charge_kcal_mol"):
        pattern = PATTERNS[key] + r"\s+KCAL/MOL"
        matches = re.findall(pattern, text)
        if matches:
            row[key] = float(matches[0])
    row["status"] = "SUCCEEDED" if "EXECUTION OF GAMESS TERMINATED NORMALLY" in text else "FAILED"
    row["warnings_no_bond"] = len(re.findall(r"FORMS NO BOND", text))
    row["scf_unconverged"] = "SCF IS UNCONVERGED" in text
    if "potential_kcal_mol" in row and "total_kcal_mol" in row:
        row["total_potential_closure_kcal_mol"] = row["total_kcal_mol"] - row["potential_kcal_mol"]
    return row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = {}
    for leg, endpoint in (("bound", 0), ("bound", 10), ("solution", 0), ("solution", 10)):
        directory = args.root / leg / f"endpoint_{endpoint:02d}"
        path = directory / "matched2.stdout" if endpoint == 0 and leg == "bound" else directory / "matched.stdout"
        rows[f"{leg}_{endpoint}"] = {"log": str(path), **parse(path)}
    binding = {}
    for endpoint in (0, 10):
        binding[str(endpoint)] = rows[f"bound_{endpoint}"]["total_kcal_mol"] - rows[f"solution_{endpoint}"]["total_kcal_mol"]
    report = {
        "schema": "quanpol-endpoint-pilot-audit-v1",
        "status": "COMPLETE" if all(row["status"] == "SUCCEEDED" for row in rows.values()) else "INCOMPLETE",
        "rows": rows,
        "same_frame_endpoint_difference_kcal_mol": binding,
        "endpoint_double_difference_kcal_mol": binding["10"] - binding["0"],
        "interpretation": [
            "This is a single-frame finite-shell DFTB3/AMBER QuanPol potential-energy diagnostic, not FEP or a binding free energy.",
            "The bound and solution shells are different environments; their absolute energies cannot be compared as a calibrated affinity without Hamiltonian and boundary validation.",
            "All four runs have total-potential closure; SCF convergence and normal termination are required gates.",
            "FORMS NO BOND warnings include zero-force constraint bonds and unbonded ions; they are not by themselves missing-atom or failed-job evidence, but require topology/term interpretation.",
            "The result is useful to identify QM/MM interface and state-energy terms, but it must not replace tPIE or FEP in the ranking pipeline.",
        ],
        "next_gate": "validate FFDATA MM-only energy against the same-frame OpenMM endpoint Hamiltonian before any multi-frame QM/MM correction",
    }
    write_json(args.output, report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
