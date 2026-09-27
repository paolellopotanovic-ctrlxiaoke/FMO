"""Render a conservative GAMESS FMO2/PIEDA input deck from prepared fragments."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import math

from .fragmentation import Fragment
from .structure import Atom, Ligand, Protein

_ELEMENT_Z = {
    "H": 1, "HE": 2, "LI": 3, "BE": 4, "B": 5, "C": 6, "N": 7,
    "O": 8, "F": 9, "NE": 10, "NA": 11, "MG": 12, "AL": 13,
    "SI": 14, "P": 15, "S": 16, "CL": 17, "AR": 18, "K": 19,
    "CA": 20, "BR": 35, "I": 53,
}
_ELEMENT_VALENCE = {
    "H": 1, "C": 4, "N": 5, "O": 6, "F": 7, "P": 5, "S": 6,
    "CL": 7, "BR": 7, "I": 7,
}

def _dftb_symbol(element: str) -> str:
    normalized = element.strip().upper()
    if normalized not in _ELEMENT_Z:
        raise ValueError(f"unsupported element {element!r}")
    return normalized[:1] + normalized[1:].lower()


def _format_namelist_array(name: str, values: list[str], *, width: int = 72,
                           comma: bool = True) -> list[str]:
    prefix = f"      {name}(1)="
    output: list[str] = []
    current = prefix
    for index, value in enumerate(values):
        token = value + ("," if comma and index < len(values) - 1 else "")
        separator = "" if current == prefix else " "
        if len(current) + len(separator) + len(token) > width and current != prefix:
            output.append(current)
            current = "                " + token
        else:
            current += separator + token
    output.append(current)
    return output


@dataclass(frozen=True)
class DeckAtom:
    atom: Atom
    fragment_id: int
    gamess_index: int


def _fragment_name(name: str, index: int) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9]", "", name)[:8]
    return cleaned or f"FRG{index:04d}"


def build_deck_atoms(protein: Protein, ligand: Ligand, fragments: list[Fragment]) -> list[DeckAtom]:
    if [fragment.fragment_id for fragment in fragments] != list(range(1, len(fragments) + 1)):
        raise ValueError("fragment IDs must be consecutive and one-based")
    if sum(fragment.kind == "ligand" for fragment in fragments) != 1:
        raise ValueError("exactly one ligand fragment is required")
    if any(fragment.multiplicity != 1 for fragment in fragments):
        raise ValueError("this RHF renderer supports only singlet fragments")
    if any(not fragment.atom_indices for fragment in fragments):
        raise ValueError("empty fragments are not allowed")
    protein_by_index = {atom.index: atom for atom in protein.atoms}
    ligand_by_index = {atom.index: atom for atom in ligand.atoms}
    if len(protein_by_index) != len(protein.atoms) or len(ligand_by_index) != len(ligand.atoms):
        raise ValueError("duplicate source atom indices")
    ligand_fragment = next(fragment for fragment in fragments if fragment.kind == "ligand")
    protein_fragment_ids = {
        atom_index: fragment.fragment_id
        for fragment in fragments if fragment.kind != "ligand"
        for atom_index in fragment.atom_indices
    }
    result: list[DeckAtom] = []
    for fragment in fragments:
        for atom_index in fragment.atom_indices:
            source = ligand_by_index if fragment.kind == "ligand" else protein_by_index
            if atom_index not in source:
                raise ValueError(f"fragment {fragment.name!r} references missing atom {atom_index}")
            atom = source[atom_index]
            result.append(DeckAtom(atom, fragment.fragment_id, len(result) + 1))
    if not any(item.fragment_id == ligand_fragment.fragment_id for item in result):
        raise ValueError("ligand fragment contains no atoms")
    if len(protein_fragment_ids) != sum(len(fragment.atom_indices) for fragment in fragments if fragment.kind != "ligand"):
        raise ValueError("duplicate atom index in protein fragments")
    if set(protein_fragment_ids) != set(protein_by_index):
        raise ValueError("all protein atoms must belong to a fragment")
    if set(ligand_fragment.atom_indices) != set(ligand_by_index) or len(ligand_fragment.atom_indices) != len(ligand.atoms):
        raise ValueError("all ligand atoms must belong to its fragment exactly once")
    return result


def validate_fragment_electron_counts(deck_atoms: list[DeckAtom], fragments: list[Fragment],
                                      fragment_bonds: list[tuple[int, int]] | None = None,
                                      charges: list[int] | None = None,
                                      multiplicities: list[int] | None = None,
                                      valence: bool = False) -> None:
    """Reject open-shell-invalid RHF fragments before GAMESS is launched."""
    atom_counts: dict[int, int] = {fragment.fragment_id: 0 for fragment in fragments}
    for item in deck_atoms:
        try:
            table = _ELEMENT_VALENCE if valence else _ELEMENT_Z
            atom_counts[item.fragment_id] += table[item.atom.element.upper()]
        except KeyError as error:
            raise ValueError(f"unsupported element {item.atom.element!r}") from error
    effective_charges = charges or [fragment.charge for fragment in fragments]
    effective_multiplicities = multiplicities or [fragment.multiplicity for fragment in fragments]
    atom_fragments = {item.gamess_index: item.fragment_id for item in deck_atoms}
    for fragment in fragments:
        electrons = atom_counts[fragment.fragment_id] - effective_charges[fragment.fragment_id - 1]
        for left, right in fragment_bonds or []:
            if atom_fragments.get(abs(left)) == fragment.fragment_id:
                electrons -= 1
            if atom_fragments.get(right) == fragment.fragment_id:
                electrons += 1
        if electrons < 1:
            raise ValueError(f"fragment {fragment.name!r} has non-positive electron count (electrons={electrons})")
        multiplicity = effective_multiplicities[fragment.fragment_id - 1]
        if multiplicity < 1 or (electrons + multiplicity) % 2 == 0:
            raise ValueError(
                f"fragment {fragment.name!r} has incompatible electrons/multiplicity "
                f"(electrons={electrons}, multiplicity={multiplicity})"
            )
        if multiplicity == 1 and electrons % 2:
            raise ValueError(
                f"fragment {fragment.name!r} has {electrons} electrons but singlet multiplicity; "
                "structure likely needs hydrogens/capping or a charge/multiplicity assignment"
            )


def render_fmo2_pieda(
    protein: Protein,
    ligand: Ligand,
    fragments: list[Fragment],
    *,
    title: str = "SophosQM-style FMO2 PIEDA",
    dftb_parameter_dir: str | Path,
    nacut: int = 0,
    fragment_bonds: list[tuple[int, int]] | None = None,
) -> str:
    """Create a syntax-level GAMESS FMO2 deck.

    Fragment coordinates and atom membership are explicit. Charges are taken from
    the fragment manifest; the caller is responsible for chemically validating
    protonation, caps, and residue boundary charges before production use.
    """
    deck_atoms = build_deck_atoms(protein, ligand, fragments)
    if not dftb_parameter_dir:
        raise ValueError("dftb_parameter_dir is required for reproducible GAMESS paths")
    max_fragment = max(fragment.fragment_id for fragment in fragments)
    names = [_fragment_name(fragment.name, fragment.fragment_id) for fragment in fragments]
    if len(set(names)) != len(names):
        raise ValueError("fragment names must be unique after GAMESS truncation")
    charges = [fragment.charge for fragment in fragments]
    indat_values = [0]
    atom_offset = 0
    for fragment in fragments:
        fragment_size = len(fragment.atom_indices)
        first_atom = atom_offset + 1
        last_atom = atom_offset + fragment_size
        indat_values.append(first_atom if first_atom == last_atom else first_atom)
        if last_atom > first_atom:
            indat_values.append(-last_atom)
        indat_values.append(0)
        atom_offset += fragment_size
    indat = [str(value) for value in indat_values]
    coordinates = []
    for item in deck_atoms:
        element = item.atom.element.upper()
        if element not in _ELEMENT_Z:
            raise ValueError(f"unsupported element {element!r} in atom {item.gamess_index}")
        x, y, z = item.atom.xyz
        if not all(math.isfinite(value) for value in item.atom.xyz):
            raise ValueError("coordinates must be finite")
        coordinates.append(f"{element:<2s} {_ELEMENT_Z[element]:3d} {x:14.8f} {y:14.8f} {z:14.8f}")
    species = sorted({item.atom.element.upper() for item in deck_atoms})
    fragment_bonds = fragment_bonds or []
    atom_fragments = {item.gamess_index: item.fragment_id for item in deck_atoms}
    for left, right in fragment_bonds:
        if abs(left) not in atom_fragments or right not in atom_fragments:
            raise ValueError("fragment bond references an atom outside the deck")
        if atom_fragments[abs(left)] == atom_fragments[right]:
            raise ValueError("fragment bond must connect two different fragments")
    atom_counts = {fragment.fragment_id: 0 for fragment in fragments}
    for item in deck_atoms:
        atom_counts[item.fragment_id] += _ELEMENT_VALENCE[item.atom.element.upper()]
    multiplicities = [fragment.multiplicity for fragment in fragments]
    validate_fragment_electron_counts(deck_atoms, fragments, fragment_bonds, charges, multiplicities,
                                      valence=True)
    if any(abs(left) < 1 or right < 1 or abs(left) == right for left, right in fragment_bonds):
        raise ValueError("fragment bond atom indices must be distinct positive integers")
    lines = [
        f"! {title}",
        (" $CONTRL RUNTYP=ENERGY SCFTYP=ROHF ISPHER=1 UNITS=ANGS MAXIT=200 $END"
         if any(value != 1 for value in multiplicities)
         else " $CONTRL RUNTYP=ENERGY SCFTYP=RHF ISPHER=1 UNITS=ANGS MAXIT=200 $END"),
        " $SYSTEM MWORDS=500 $END",
        " $SCF DIRSCF=.TRUE. NPUNCH=0 CONV=1.0E-8 DIIS=.FALSE. SOSCF=.TRUE. $END",
        " $BASIS GBASIS=DFTB $END",
        " $DFTB SCC=.TRUE. DFTB3=.TRUE. DAMPXH=.TRUE. DAMPEX=4.00 $END",
        *([" $DFTBSK"] + [
            f"   {first} {second} {dftb_parameter_dir}/{_dftb_symbol(first)}-{_dftb_symbol(second)}.skf"
            for first in species for second in species
        ] + [" $END"]),
        " $FMOPRP NPRINT=9 IPIEDA=1 MAXIT=200 CONV=1.0E-7 NGUESS=2 MODORB=3",
        " CNVDMP=50 MCONV(2)=530 MCONV(4)=530 $END",
        " $DFT DC=.TRUE. IDCVER=4 $END",
        " $PCM SOLVNT=WATER IEF=-10 ICOMP=0 ICAV=1 IDISP=1 IFMO=-1 $END",
        " $PCMCAV RADII=VANDW $END",
        " $TESCAV NTSALL=60 $END",
        f" $FMO NFRAG={max_fragment} NBODY=2 NACUT={nacut} MODMOL=1 MOLFRG(1)={max_fragment}",
    ]
    lines.extend(_format_namelist_array("ICHARG", [str(charge) for charge in charges]))
    if any(value != 1 for value in multiplicities):
        lines.extend(_format_namelist_array("MULT", [str(value) for value in multiplicities]))
        lines.extend(_format_namelist_array("SCFFRG", ["ROHF" if value != 1 else "RHF" for value in multiplicities]))
    lines.extend(_format_namelist_array("FRGNAM", names))
    lines.extend(_format_namelist_array("INDAT", indat, width=72, comma=False))
    homos = ([" $FMOHYB", " HOP_C 4 4",
           " 1 0  0.562060  0.000000  0.000000  0.827096",
           " 0 1  0.562060  0.779794  0.000000 -0.275699",
           " 0 1  0.562060 -0.389897  0.675322 -0.275698",
           " 0 1  0.562060 -0.389897 -0.675322 -0.275698",
           " $END"] if fragment_bonds else [])
    bond_lines = [f" {left} {right} HOP_C" for left, right in fragment_bonds]
    lines.extend([
        " $END",
        *homos,
        *( [" $FMOBND", *bond_lines, " $END"] if fragment_bonds else []),
        " $FMOXYZ",
        *coordinates,
        " $END",
        " $DATA",
        title,
        "C1",
        *(f"{element}-1 {atomic_number}" for element, atomic_number in sorted(
            {item.atom.element.upper(): _ELEMENT_Z[item.atom.element.upper()] for item in deck_atoms}.items()
        )),
        " $END",
        "",
    ])
    return "\n".join(lines)


def write_fmo2_pieda_input(path: str | Path, protein: Protein, ligand: Ligand,
                           fragments: list[Fragment], **kwargs: object) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(render_fmo2_pieda(protein, ligand, fragments, **kwargs), encoding="utf-8")
    return output
