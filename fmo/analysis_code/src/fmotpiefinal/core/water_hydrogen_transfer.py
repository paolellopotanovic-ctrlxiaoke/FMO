"""Preserve validated source-water orientation after protein protonation."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def water_key(line: str) -> tuple[str, str, str]:
    return line[21], line[22:27], line[17:20]


def xyz(line: str) -> tuple[float, float, float]:
    return tuple(float(line[position:position + 8]) for position in (30, 38, 46))


def waters(path: Path) -> dict[tuple[str, str, str], dict[str, str]]:
    grouped: dict[tuple[str, str, str], dict[str, str]] = {}
    for line in path.read_text().splitlines(keepends=True):
        if not line.startswith("HETATM") or line[17:20] != "HOH":
            continue
        atom = line[12:16].strip()
        selected = grouped.setdefault(water_key(line), {})
        if atom in selected:
            raise ValueError(f"duplicate water atom {water_key(line)}:{atom}")
        selected[atom] = line
    return grouped


def valid_water(atoms: dict[str, str]) -> bool:
    if set(atoms) != {"O", "H1", "H2"} or any(
            atoms[name][76:78].strip().upper() != ("O" if name == "O" else "H")
            for name in ("O", "H1", "H2")):
        return False
    try:
        oxygen = xyz(atoms["O"])
        hydrogen_a, hydrogen_b = xyz(atoms["H1"]), xyz(atoms["H2"])
        left, right = math.dist(oxygen, hydrogen_a), math.dist(oxygen, hydrogen_b)
        cross = math.dist(hydrogen_a, hydrogen_b)
    except (ValueError, OverflowError):
        return False
    return all(math.isfinite(value) for value in (left, right, cross)) and (
        0.7 <= left <= 1.25 and 0.7 <= right <= 1.25 and 1.2 <= cross <= 1.95)


def preserve_source_waters(source: Path, prepared: Path) -> dict:
    original = waters(source)
    target = waters(prepared)
    if original and not target:
        raise ValueError("prepared receptor lost every source water")
    selected = {}
    skipped = {}
    for key, atoms in target.items():
        incoming = original.get(key)
        if incoming is None or not valid_water(incoming):
            skipped[":".join(part.strip() for part in key)] = "source_waters_missing_or_invalid"
            continue
        if not valid_water(atoms) or math.dist(xyz(incoming["O"]), xyz(atoms["O"])) > 0.01:
            raise ValueError(f"source/prepared water oxygen identity or geometry changed: {key}")
        selected[key] = incoming
    lines = []
    for line in prepared.read_text().splitlines(keepends=True):
        if line.startswith("HETATM") and line[17:20] == "HOH":
            name = line[12:16].strip()
            if name in ("H1", "H2") and water_key(line) in selected:
                line = line[:30] + selected[water_key(line)][name][30:54] + line[54:]
        lines.append(line)
    report = {"schema": "fmo-validated-source-water-hydrogens-v1",
              "source_sha256": digest(source), "prepared_pretransfer_sha256": digest(prepared),
              "n_source_waters": len(original), "n_prepared_waters": len(target),
              "n_preserved_waters": len(selected), "n_pqr_fallback_waters": len(skipped),
              "pqr_fallback_reasons": skipped}
    content = "".join(lines)
    if len(content) != len(prepared.read_text()):
        raise ValueError("source-water transfer changed PDB length")
    prepared.write_text(content)
    report["prepared_sha256"] = digest(prepared)
    report_file = prepared.with_suffix(".water_hydrogens.json")
    report_file.write_text(json.dumps(report, indent=2) + "\n")
    return report
