#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from sophosqm_baseline.cohort import TARGETS, build_cohort, file_sha256, write_cohort


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structure-root", default="assets/structures/jacs_set")
    parser.add_argument("--output", default="data/fmo87_cohort.csv")
    args = parser.parse_args()
    records = build_cohort(args.structure_root)
    inputs = {filename: file_sha256(Path(args.structure_root) / filename)
              for _, ( _, filename, _, _) in TARGETS.items()}
    inputs.update({filename: file_sha256(Path(args.structure_root) / filename)
                   for _, (_, _, filename, _) in TARGETS.items()})
    write_cohort(args.output, records, inputs)
    counts = {target: sum(record.target == target for record in records) for target in TARGETS}
    print(counts, "total", len(records), "output", args.output)


if __name__ == "__main__":
    main()
