"""Report the QM/MM minus MM endpoint correction on the same pilot frame."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from run_matched_openfe import write_json
from report_quanpol_endpoint_pilot import parse


def first_potential(path: Path) -> float:
    row = parse(path)
    if row["status"] != "SUCCEEDED" or row["scf_unconverged"]:
        raise ValueError(f"failed or unconverged calculation: {path}")
    return row["potential_kcal_mol"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--qmmm", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.qmmm.read_text())
    qm = {key: row["potential_kcal_mol"] for key, row in source["rows"].items()
          if row["status"] == "SUCCEEDED" and not row["scf_unconverged"]}
    mm = {}
    for leg, endpoint in (("bound", 0), ("bound", 10), ("solution", 0), ("solution", 10)):
        mm[f"{leg}_{endpoint}"] = first_potential(args.root / leg / f"endpoint_{endpoint:02d}" / "mm_matched.stdout")
    qm_binding = {str(endpoint): qm[f"bound_{endpoint}"] - qm[f"solution_{endpoint}"] for endpoint in (0, 10)}
    mm_binding = {str(endpoint): mm[f"bound_{endpoint}"] - mm[f"solution_{endpoint}"] for endpoint in (0, 10)}
    qm_double = qm_binding["10"] - qm_binding["0"]
    mm_double = mm_binding["10"] - mm_binding["0"]
    result = {"schema": "same-frame-qmmm-mm-correction-v1", "qm_mm_total_potential_kcal_mol": qm,
              "mm_only_potential_kcal_mol": mm, "qmmm_binding_difference_kcal_mol": qm_binding,
              "mm_binding_difference_kcal_mol": mm_binding,
              "qmmm_double_difference_kcal_mol": qm_double,
              "mm_double_difference_kcal_mol": mm_double,
              "qm_minus_mm_double_difference_kcal_mol": qm_double - mm_double,
              "status": "DIAGNOSTIC_COMPLETE",
              "interpretation": [
                  "The MM-only values use the corrected LJSIGMA=1 QuanPol export.",
                  "This is a single-frame finite-shell diagnostic, not a free energy.",
                  "The previous -0.225 claim used mismatched 1-4 scaling and is superseded by this matched calculation.",
                  "This does not prove the QM correction is useless across frames; it establishes that multi-frame sampling is required before ranking claims.",
              ]}
    write_json(args.output, result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
