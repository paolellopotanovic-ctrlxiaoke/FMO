"""Create a reproducible GitHub source archive without local build/run artifacts."""
from __future__ import annotations

import argparse
import gzip
import hashlib
from pathlib import Path
import tarfile


EXCLUDED = {"__pycache__", ".pytest_cache", ".git", "build", "dist", "runs",
            "data", "datasets", "software", ".venv", "envs", "archive"}


def source_files(root: Path) -> list[Path]:
    selected = []
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root)
        if any(part in EXCLUDED or part.endswith(".egg-info") for part in relative.parts):
            continue
        if path.is_symlink():
            raise ValueError(f"release archive cannot contain symlinks: {relative}")
        if path.is_file() and path.suffix not in {".pyc", ".log"}:
            selected.append(path)
    return selected


def archive(root: Path, output: Path) -> str:
    if output.resolve().is_relative_to(root.resolve()):
        raise ValueError("write the archive outside the source directory")
    with output.open("wb") as destination:
        with gzip.GzipFile(filename="", mode="wb", fileobj=destination, mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode="w") as package:
                for path in source_files(root):
                    info = package.gettarinfo(str(path), "FMOtpiefinal/" + str(path.relative_to(root)))
                    info.uid = info.gid = info.mtime = 0
                    info.uname = info.gname = ""
                    info.mode = 0o755 if path.suffix == ".sh" else 0o644
                    with path.open("rb") as content:
                        package.addfile(info, content)
    return hashlib.sha256(output.read_bytes()).hexdigest()


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=root.parent / "FMOtpiefinal-0.1.0rc1-github.tar.gz")
    arguments = parser.parse_args()
    print(archive(root, arguments.output), arguments.output.name)
