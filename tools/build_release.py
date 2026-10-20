"""Build source and wheel distributions using the configured setuptools backend."""
from pathlib import Path
import os

from setuptools import build_meta


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1]
    os.chdir(root)
    destination = root / "dist"
    destination.mkdir(exist_ok=True)
    print(build_meta.build_sdist(str(destination)))
    print(build_meta.build_wheel(str(destination)))
