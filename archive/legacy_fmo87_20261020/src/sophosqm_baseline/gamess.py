"""GAMESS job preparation and execution without bundling GAMESS itself."""
from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
import hashlib
import json


@dataclass(frozen=True)
class GamessConfig:
    executable: str
    input_path: str
    method: str = "DFTB3"
    basis: str = "DFTB3/3OB-3-1"
    ncores: int = 1


def sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_job_manifest(output_dir: str | Path, config: GamessConfig, structure_manifest: str,
                       fragment_manifest: str, metadata: dict[str, object]) -> Path:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    payload = {"config": asdict(config), "structure_manifest": structure_manifest,
               "fragment_manifest": fragment_manifest, "metadata": metadata}
    path = output / "job_manifest.json"
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return path
