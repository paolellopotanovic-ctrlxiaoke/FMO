#!/usr/bin/env python
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sophosqm_baseline.preparation import prepare_protein_pdb2pqr

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--drop-water", action="store_true",
                        help="remove crystallographic waters before PDB2PQR/PROPKA")
    args = parser.parse_args()
    print(json.dumps(prepare_protein_pdb2pqr(
        args.input, args.output, keep_water=not args.drop_water
    ), indent=2))


if __name__ == "__main__":
    main()
