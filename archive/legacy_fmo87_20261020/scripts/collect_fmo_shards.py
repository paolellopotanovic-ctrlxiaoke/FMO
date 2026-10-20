#!/usr/bin/env python
"""Collect parallel FMO shard outputs into auditable target directories."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))


EXPECTED_COUNTS = {"CDK2": 16, "JNK1": 21, "P38": 34, "TYK2": 16}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shard-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()

    shard_dirs = sorted(
        (path for path in args.shard_root.iterdir() if path.is_dir() and path.name.isdigit()),
        key=lambda path: int(path.name),
    )
    frames = []
    for shard in shard_dirs:
        features = shard / "fmo_features.csv"
        if not features.is_file():
            raise ValueError(f"missing shard features: {features}")
        frame = pd.read_csv(features, dtype={"ligand_name": str})
        frame["shard"] = int(shard.name)
        frames.append(frame)
    if not frames:
        raise ValueError(f"no numbered shard directories in {args.shard_root}")
    combined = pd.concat(frames, ignore_index=True)
    identity = ["target", "ligand_name"]
    if combined.duplicated(identity).any():
        duplicates = combined.loc[combined.duplicated(identity, keep=False), identity + ["shard"]]
        raise ValueError(f"duplicate shard identities:\n{duplicates.to_string(index=False)}")
    counts = combined["target"].value_counts().to_dict()
    complete = (
        len(combined) == 87
        and counts == EXPECTED_COUNTS
        and combined["status"].eq("SUCCEEDED").all()
        and combined["tpie_kcal_mol"].notna().all()
    )
    if args.require_complete and not complete:
        failures = combined.loc[combined["status"].ne("SUCCEEDED"), identity + ["status", "shard"]]
        raise ValueError(f"incomplete shard collection: {counts}; failures={failures.to_dict('records')}")

    args.output_root.mkdir(parents=True, exist_ok=True)
    for target in EXPECTED_COUNTS:
        rows = combined.loc[combined["target"] == target]
        target_dir = args.output_root / target
        if target_dir.exists():
            shutil.rmtree(target_dir)
        target_dir.mkdir(parents=True)
        rows.drop(columns="shard").to_csv(target_dir / "fmo_features.csv", index=False)
        for row in rows.itertuples(index=False):
            key = f"{row.target}_{row.ligand_name}"
            shutil.copytree(args.shard_root / str(row.shard) / "minimized" / key,
                            target_dir / "minimized" / key)
            shutil.copytree(args.shard_root / str(row.shard) / "fmo" / key,
                            target_dir / "fmo" / key)

    summary = {
        "n_compounds": len(combined),
        "counts": counts,
        "n_succeeded": int(combined["status"].eq("SUCCEEDED").sum()),
        "complete": complete,
        "n_shards": len(shard_dirs),
        "output_root": str(args.output_root),
    }
    (args.output_root / "collection_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
