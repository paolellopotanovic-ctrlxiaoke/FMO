"""Parsers for GAMESS FMO PIEDA pair tables."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re


@dataclass(frozen=True)
class PairInteraction:
    fragment_i: int
    fragment_j: int
    value_kcal_mol: float
    electrostatic_kcal_mol: float | None = None
    exchange_repulsion_kcal_mol: float | None = None
    charge_transfer_kcal_mol: float | None = None
    dispersion_kcal_mol: float | None = None
    solvation_kcal_mol: float | None = None


_FRAGMENT_SECTION = re.compile(
    r"CONV\s*\n\s*={70,90}\s*\n(.*?)\n\s*Close fragment pairs",
    re.IGNORECASE | re.DOTALL,
)
_PIEDA_SECTION = re.compile(r"\s-{100,}\s*\n(.*?)\n\s*Total energy", re.I | re.S)
_FMO_PROPERTIES = re.compile(
    r"Two-body FMO properties\.\s*=+.*?\n(.*?)(?=\n\s*Total energy of the molecule:|\Z)",
    re.IGNORECASE | re.DOTALL,
)
_SELECTED_TIE = re.compile(
    r"^\s*(\d+)\([^)]*\):.*final TIE\s+([-+]?\d+(?:\.\d+)?(?:[EedD][-+]?\d+)?)",
    re.IGNORECASE | re.MULTILINE,
)


def _number(text: str) -> float:
    return float(text.strip().replace("D", "E").replace("d", "e"))


def _parse_pair_line(line: str) -> PairInteraction | None:
    if len(line) < 60:
        return None
    try:
        fragment_i = int(line[:5])
        fragment_j = int(line[5:10])
        values = {
            "total": _number(line[50:60]),
            "electrostatic": _number(line[60:70]),
            "exchange": _number(line[70:79]),
            "charge_transfer": _number(line[79:88]),
            "dispersion": _number(line[88:97]),
            "solvation": _number(line[97:106]) if len(line) >= 106 else 0.0,
        }
    except ValueError:
        return None
    return PairInteraction(
        fragment_i=fragment_i,
        fragment_j=fragment_j,
        value_kcal_mol=values["total"],
        electrostatic_kcal_mol=values["electrostatic"],
        exchange_repulsion_kcal_mol=values["exchange"],
        charge_transfer_kcal_mol=values["charge_transfer"],
        dispersion_kcal_mol=values["dispersion"],
        solvation_kcal_mol=values["solvation"],
    )


def parse_pie(path: str | Path) -> list[PairInteraction]:
    """Parse the fixed-width pair table in a GAMESS FMO PIEDA output.

    GAMESS PIEDA pair rows store the total interaction in columns 51-60,
    followed by ES, EX, CT+mix, DI, and solvent components. The displayed
    total is used once; component columns are not added a second time.
    """
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    sections = list(_PIEDA_SECTION.finditer(text))
    if sections:
        records = _parse_fixed_width_sections(sections)
    else:
        records = _parse_two_body_fmo_sections(text)
    if not records:
        raise ValueError(f"No GAMESS FMO PIEDA table found in {path}")

    seen: set[tuple[int, int]] = set()
    for record in records:
        pair = tuple(sorted((record.fragment_i, record.fragment_j)))
        if pair in seen:
            raise ValueError(f"Duplicate GAMESS PIEDA pair row for fragments {pair}")
        seen.add(pair)
    return records


def _parse_fixed_width_sections(sections: list[re.Match[str]]) -> list[PairInteraction]:
    records: list[PairInteraction] = []
    for section in sections:
        for line in section.group(1).splitlines():
            record = _parse_pair_line(line)
            if record is not None:
                records.append(record)
    return records


def _parse_two_body_fmo_sections(text: str) -> list[PairInteraction]:
    """Parse the standard GAMESS FMO table emitted by FMO2/PIEDA runs."""
    records: list[PairInteraction] = []
    for section in _FMO_PROPERTIES.finditer(text):
        for line in section.group(1).splitlines():
            fields = line.split()
            if len(fields) < 14 or not fields[0].isdigit() or not fields[1].isdigit():
                continue
            try:
                values = [float(value.replace("D", "E").replace("d", "e")) for value in fields[4:14]]
            except ValueError:
                continue
            records.append(PairInteraction(
                fragment_i=int(fields[0]), fragment_j=int(fields[1]),
                value_kcal_mol=values[4], electrostatic_kcal_mol=values[5],
                exchange_repulsion_kcal_mol=values[6], charge_transfer_kcal_mol=values[7],
                dispersion_kcal_mol=values[8], solvation_kcal_mol=values[9],
            ))
    return records


def parse_fragment_names(path: str | Path) -> dict[int, str]:
    """Read fragment IDs and names from GAMESS's converged FMO listing."""
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    sections = list(_FRAGMENT_SECTION.finditer(text))
    fragments: dict[int, str] = {}
    for section in sections:
        for line in section.group(1).splitlines():
            fields = line.split()
            if len(fields) >= 2 and fields[0].isdigit():
                fragments[int(fields[0])] = fields[1]
    if not fragments:
        for section in re.finditer(r"Fragment statistics\s*\n(.*?)(?=\n\s*Close fragment pairs)", text, re.I | re.S):
            for line in section.group(1).splitlines():
                fields = line.split()
                if len(fields) >= 2 and fields[0].isdigit() and fields[1].lower() != "name":
                    fragments[int(fields[0])] = fields[1]
    if not fragments:
        raise ValueError(f"GAMESS fragment listing in {path} contains no fragment IDs")
    return fragments


def aggregate_tpie(records: list[PairInteraction], ligand_fragment: int) -> float:
    """Sum ligand-protein pair totals, counting each unordered pair once."""
    selected = [
        record.value_kcal_mol
        for record in records
        if (record.fragment_i == ligand_fragment) != (record.fragment_j == ligand_fragment)
    ]
    if not selected:
        raise ValueError(f"No ligand-protein PIEDA records found for fragment {ligand_fragment}")
    return float(sum(selected))


def find_selected_tie(path: str | Path, ligand_fragment: int) -> float | None:
    """Return GAMESS's official selected-fragment TIE, if printed."""
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    matches = [
        _number(value)
        for fragment, value in _SELECTED_TIE.findall(text)
        if int(fragment) == ligand_fragment
    ]
    if len(matches) > 1:
        raise ValueError(f"multiple selected TIE values for fragment {ligand_fragment}")
    return matches[0] if matches else None
