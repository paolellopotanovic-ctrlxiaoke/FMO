#!/usr/bin/env python
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import urllib.request

from boltz.main import CCD_URL, download_boltz2

EXPECTED_SHA256 = {
    "boltz2_conf.ckpt": "090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1",
    "boltz2_aff.ckpt": "dcc5cd3722b1c9eaa34267e4ae32f55cbbf1963f4c19319381ccfa30fdd2ca9e",
    "ccd.pkl": "2d3b2f03a3c5665944adba51e33263511e51b21c9cd05d902f9c4b7c1e58d2f4",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", default="engines/boltz-cache")
    args = parser.parse_args()
    cache = Path(args.cache).resolve()
    cache.mkdir(parents=True, exist_ok=True)
    download_boltz2(cache)
    ccd = cache / "ccd.pkl"
    if not ccd.is_file():
        print(f"Downloading Boltz CCD dictionary to {ccd}")
        urllib.request.urlretrieve(CCD_URL, ccd)
    for name in ("mols.tar", "mols", *EXPECTED_SHA256):
        path = cache / name
        if not path.exists():
            raise FileNotFoundError(path)
    for name, expected_hash in EXPECTED_SHA256.items():
        actual_hash = sha256(cache / name)
        if actual_hash != expected_hash:
            raise ValueError(f"{name} SHA-256 mismatch: {actual_hash} != {expected_hash}")
    print(f"Boltz-2 engine cache is ready: {cache}")


if __name__ == "__main__":
    main()
