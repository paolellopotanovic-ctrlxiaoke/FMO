#!/usr/bin/env python
"""Merge and validate per-target outputs for the frozen FMO87 cohort."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


EXPECTED_COUNTS = {"CDK2": 16, "JNK1": 21, "P38": 34, "TYK2": 16}


def merge(output: Path, target_dirs: list[Path], kind: str) -> pd.DataFrame:
    frames = []
    for target in EXPECTED_COUNTS:
        path = next((item for item in target_dirs if item.name == target), None)
        if path is None:
            raise ValueError(f"missing {kind} target directory: {target}")
        filename = "fmo_features.csv" if kind == "FMO" else "predictions.csv"
        candidate = path / filename
        if not candidate.is_file():
            raise ValueError(f"missing {kind} output: {candidate}")
        frame = pd.read_csv(candidate, dtype={"ligand_name": str, "ligand_id": str})
        frames.append(frame)
    merged = pd.concat(frames, ignore_index=True)
    identity_columns = ["target", "ligand_name"] if kind == "FMO" else ["target", "ligand_id"]
    if merged.duplicated(identity_columns).any():
        duplicates = merged.loc[merged.duplicated(identity_columns, keep=False), identity_columns]
        raise ValueError(f"duplicate merged {kind} identities:\n{duplicates.to_string(index=False)}")
    counts = merged["target"].value_counts().to_dict()
    if counts != EXPECTED_COUNTS or len(merged) != 87:
        raise ValueError(f"invalid merged {kind} counts: {counts}")
    if kind == "FMO":
        statuses = merged["status"].value_counts(dropna=False).to_dict()
        if merged["status"].ne("SUCCEEDED").any():
            failures = merged.loc[merged["status"].ne("SUCCEEDED"), identity_columns + ["status"]]
            raise ValueError(f"failed FMO records:\n{failures.to_string(index=False)}")
        if not pd.to_numeric(merged["tpie_kcal_mol"], errors="coerce").notna().all():
            raise ValueError("nonfinite or missing FMO tPIE")
    else:
        if not merged["status"].fillna("").str.startswith("SUCCEEDED").all():
            failures = merged.loc[~merged["status"].fillna("").str.startswith("SUCCEEDED"), identity_columns + ["status"]]
            raise ValueError(f"failed Boltz2 records:\n{failures.to_string(index=False)}")
        numeric = pd.to_numeric(merged["affinity_pred_value"], errors="coerce")
        if not numeric.notna().all():
            raise ValueError("nonfinite or missing Boltz2 affinity")
    output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output, index=False)
    summary = {
        "kind": kind,
        "n_records": len(merged),
        "counts_by_target": counts,
        "status_counts": merged["status"].value_counts(dropna=False).astype(int).to_dict(),
        "output": str(output),
    }
    (output.parent / f"{kind.lower()}_merge_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return merged


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fmo-targets", nargs="*", type=Path)
    parser.add_argument("--boltz-targets", nargs="*", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if not args.fmo_targets and not args.boltz_targets:
        parser.error("at least one of --fmo-targets or --boltz-targets is required")
    if args.fmo_targets:
        merge(args.output_dir / "fmo_features.csv", args.fmo_targets, "FMO")
    if args.boltz_targets:
        merge(args.output_dir / "boltz_predictions.csv", args.boltz_targets, "Boltz2")


if __name__ == "__main__":
    main()
