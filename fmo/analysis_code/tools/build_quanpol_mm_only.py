"""Create MM-only QuanPol decks from the audited QM/MM FFDATA exports."""
from __future__ import annotations

import argparse
from pathlib import Path

from quanpol_constants import DEFAULT_SWRB_ANGSTROM


def make(source: Path, output: Path, atom_count: int):
    text = source.read_text()
    ffdata = text[text.index(" $FFDATA"):]
    capacity = atom_count + 200
    header = "\n".join([
        " $CONTRL RUNTYP=ENERGY COORD=FRAGONLY ICHARG=0 $END",
        " $SYSTEM TIMLIM=20 MWORDS=100 $END",
        " $QUANPO",
        f" NSTEP=0 NFFTYP=0 IDOPOL=0 MXLIST1={capacity} MXLIST2={capacity}",
        " LJSIGMA=1 WT14LJ=0.5 WT14CH=0.8333333333333333",
        " BUFWID1=1.0 BUFWID2=4.0 ISHIFT=0 ISWITCH=0",
        f" SWRA=1000 SWRB={DEFAULT_SWRB_ANGSTROM:g}",
        " $END",
        "",
    ])
    output.write_text(header + ffdata)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    for leg in ("bound", "solution"):
        for endpoint in (0, 10):
            directory = args.root / leg / f"endpoint_{endpoint:02d}"
            source = directory / f"matched_cdk8_{leg}_{endpoint:02d}.inp"
            if not source.exists():
                source = next(directory.glob("matched_*.inp"))
            ffdata = source.read_text()[source.read_text().index(" $FFDATA"):]
            atom_count = len(ffdata.split("COORDINATES", 1)[1].split("STOP", 1)[0].strip().splitlines())
            make(source, directory / f"mm_only_{leg}_{endpoint:02d}.inp", atom_count)
            print(leg, endpoint, atom_count, flush=True)


if __name__ == "__main__":
    main()
