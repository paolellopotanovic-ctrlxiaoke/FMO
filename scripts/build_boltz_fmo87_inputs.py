#!/usr/bin/env python
from __future__ import annotations

import argparse

from sophosqm_baseline.boltz_inputs import build_inputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", default="data/fmo87_cohort.csv")
    parser.add_argument(
        "--structure-root",
        default="data/fmo87_structures",
        help="directory containing the frozen cohort protein PDB files",
    )
    parser.add_argument("--output", default="runs/fmo87/boltz_inputs")
    parser.add_argument("--msa-dir", default="data/boltz_msa")
    parser.add_argument("--generate-msa", action="store_true")
    args = parser.parse_args()
    print(build_inputs(
        args.cohort, args.output, args.msa_dir, args.generate_msa,
        structure_root=args.structure_root,
    ))


if __name__ == "__main__":
    main()
