"""Check maintained-module closure and exported runtime fingerprints offline."""
from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path


def check(root: Path) -> dict:
    source = root / "src"
    manifest = json.loads((root / "SOURCE_PROVENANCE.json").read_text())
    for entry in manifest["files"]:
        path = root / entry["exported"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != entry["exported_sha256"]:
            raise ValueError(f"export fingerprint changed: {entry['exported']}")
    for name, expected in manifest.get("maintained_additions_sha256", {}).items():
        if hashlib.sha256((root / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"maintained addition fingerprint changed: {name}")
    modules = []
    for path in sorted((source / "fmotpiefinal").rglob("*.py")):
        relative = path.relative_to(source).with_suffix("")
        module = ".".join(relative.parts)
        package = module.rsplit(".", 1)[0]
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.ImportFrom):
                if node.level:
                    base = package.split(".")[:len(package.split(".")) - node.level + 1]
                    imported = ".".join(base + ([node.module] if node.module else []))
                else:
                    imported = node.module or ""
                names = [imported]
            elif isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            else:
                continue
            for imported in names:
                if imported.split(".")[0] in {"experiments", "validation", "sophosqm_baseline"}:
                    raise ValueError(f"research import in {path.name}: {imported}")
                if imported.startswith("fmotpiefinal."):
                    target = source.joinpath(*imported.split("."))
                    if not target.with_suffix(".py").is_file() and not (target / "__init__.py").is_file():
                        raise ValueError(f"missing runtime module in {path.name}: {imported}")
        if "/inspire/" in path.read_text() or "/global_user/" in path.read_text():
            raise ValueError(f"private absolute path in {path.name}")
        modules.append(module)
    return {"status": "PASS", "runtime_modules": len(modules),
            "verified_export_files": len(manifest["files"])}


if __name__ == "__main__":
    print(json.dumps(check(Path(__file__).resolve().parents[1]), indent=2))
