"""Render a conservative GAMESS FMO2/PIEDA input deck from prepared fragments."""
from __future__ import annotations

from dataclasses import asdict, dataclass
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


@dataclass(frozen=True)
class FMODeckOptions:
    """Explicit GAMESS controls used by production and sensitivity protocols."""

    nacut: int = 0
    dampxh: bool = True
    dampex: float | None = 4.0
    dftb_dispersion: str | None = None
    dftb_modesd: int | None = None
    dftb_modgam: int | None = None
    modorb: int = 3
    use_pcm: bool = True
    pcm_solvent: str = "WATER"
    pcm_ief: int = -10
    pcm_icomp: int = 0
    pcm_icav: int = 1
    pcm_idisp: int = 1
    pcm_ifmo: int = -1
    pcm_cavity_radii: str = "VANDW"
    pcm_tessellation: int = 60

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


def _fragment_name(name: str, index: int) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9]", "", name)[:8]
    return cleaned or f"FRG{index:04d}"


def build_deck_atoms(protein: Protein, ligand: Ligand | None, fragments: list[Fragment]) -> list[DeckAtom]:
    if [fragment.fragment_id for fragment in fragments] != list(range(1, len(fragments) + 1)):
        raise ValueError("fragment IDs must be consecutive and one-based")
    if sum(fragment.kind == "ligand" for fragment in fragments) != (ligand is not None):
        raise ValueError("ligand fragments must match the supplied ligand")
    if any(fragment.multiplicity != 1 for fragment in fragments):
        raise ValueError("this RHF renderer supports only singlet fragments")
    if any(not fragment.atom_indices for fragment in fragments):
        raise ValueError("empty fragments are not allowed")
    protein_by_index = {atom.index: atom for atom in protein.atoms}
    ligand_by_index = {atom.index: atom for atom in ligand.atoms} if ligand is not None else {}
    if len(protein_by_index) != len(protein.atoms) or (ligand is not None and len(ligand_by_index) != len(ligand.atoms)):
        raise ValueError("duplicate source atom indices")
    ligand_fragment = next((fragment for fragment in fragments if fragment.kind == "ligand"), None)
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
    if ligand_fragment is not None and not any(item.fragment_id == ligand_fragment.fragment_id for item in result):
        raise ValueError("ligand fragment contains no atoms")
    if len(protein_fragment_ids) != sum(len(fragment.atom_indices) for fragment in fragments if fragment.kind != "ligand"):
        raise ValueError("duplicate atom index in protein fragments")
    if set(protein_fragment_ids) != set(protein_by_index):
        raise ValueError("all protein atoms must belong to a fragment")
    if ligand_fragment is not None and (set(ligand_fragment.atom_indices) != set(ligand_by_index) or len(ligand_fragment.atom_indices) != len(ligand_by_index)):
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


def infer_peptide_bonds(deck_atoms: list[DeckAtom], cutoff: float = 1.65) -> list[tuple[int, int]]:
    """Infer only protein inter-residue carbonyl-C to amide-N bonds.

    This is not a general bond perception algorithm. It is deliberately narrow:
    it only accepts C--N contacts across two protein residue fragments, which is
    the boundary needed by the FACIO-like peptide-chain preparation path.
    """
    bonds: list[tuple[int, int]] = []
    carbonyls = [item for item in deck_atoms
                 if item.atom.record == "ATOM" and item.atom.atom_name.upper() == "C"]
    nitrogens = [item for item in deck_atoms
                 if item.atom.record == "ATOM" and item.atom.atom_name.upper() == "N"]
    for carbon in carbonyls:
        for nitrogen in nitrogens:
            if carbon.fragment_id == nitrogen.fragment_id:
                continue
            if carbon.atom.chain != nitrogen.atom.chain:
                continue
            try:
                residue_gap = int(nitrogen.atom.residue_number) - int(carbon.atom.residue_number)
            except ValueError:
                continue
            if residue_gap != 1:
                continue
            distance = math.sqrt(sum((a - b) ** 2 for a, b in zip(carbon.atom.xyz, nitrogen.atom.xyz)))
            if distance <= cutoff:
                bonds.append((-carbon.gamess_index, nitrogen.gamess_index))
    return bonds


def infer_shifted_ca_c_bonds(deck_atoms: list[DeckAtom], cutoff: float = 1.8) -> list[tuple[int, int]]:
    """Find the C-alpha--carbonyl-C cuts used by shifted FACIO fragments."""
    alpha_atoms = [item for item in deck_atoms
                   if item.atom.record == "ATOM" and item.atom.atom_name.upper() == "CA"]
    carbonyl_atoms = [item for item in deck_atoms
                      if item.atom.record == "ATOM" and item.atom.atom_name.upper() == "C"]
    bonds: list[tuple[int, int]] = []
    for alpha in alpha_atoms:
        for carbonyl in carbonyl_atoms:
            if alpha.fragment_id == carbonyl.fragment_id:
                continue
            if alpha.atom.chain != carbonyl.atom.chain or alpha.atom.residue_number != carbonyl.atom.residue_number:
                continue
            distance = math.sqrt(sum((a - b) ** 2 for a, b in zip(alpha.atom.xyz, carbonyl.atom.xyz)))
            if distance <= cutoff:
                bonds.append((-alpha.gamess_index, carbonyl.gamess_index))
    return bonds


def render_fmo2_pieda(
    protein: Protein,
    ligand: Ligand | None,
    fragments: list[Fragment],
    *,
    title: str = "SophosQM-style FMO2 PIEDA",
    dftb_parameter_dir: str | Path,
    nacut: int = 0,
    fragment_bonds: list[tuple[int, int]] | None = None,
    deck_options: FMODeckOptions | None = None,
    full_energy: bool = False,
) -> str:
    """Create a syntax-level GAMESS FMO2 deck.

    Fragment coordinates and atom membership are explicit. Charges are taken from
    the fragment manifest; the caller is responsible for chemically validating
    protonation, caps, and residue boundary charges before production use.
    """
    deck_atoms = build_deck_atoms(protein, ligand, fragments)
    options = deck_options or FMODeckOptions(nacut=nacut)
    if nacut != 0 and deck_options is not None and nacut != options.nacut:
        raise ValueError("nacut cannot disagree with deck_options.nacut")
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
    dftb_tokens = ["SCC=.TRUE.", "DFTB3=.TRUE."]
    if options.dampxh:
        dftb_tokens.append("DAMPXH=.TRUE.")
    if options.dampex is not None:
        dftb_tokens.append(f"DAMPEX={options.dampex:.2f}")
    if options.dftb_dispersion:
        dftb_tokens.append(f"DISP={options.dftb_dispersion.upper()}")
    if options.dftb_modesd is not None:
        dftb_tokens.append(f"MODESD={options.dftb_modesd}")
    if options.dftb_modgam is not None:
        dftb_tokens.append(f"MODGAM={options.dftb_modgam}")
    dftb_line = " $DFTB " + " ".join(dftb_tokens) + " $END"
    if len(dftb_line) > 80:
        raise ValueError(
            "GAMESS $DFTB controls exceed the 80-column input limit; "
            "omit optional default controls such as DAMPEX=4.0"
        )
    fmoprp = (
        f" $FMOPRP NPRINT=9 IPIEDA={0 if full_energy else 1} MAXIT=200 CONV=1.0E-7 "
        f"NGUESS=2 MODORB={options.modorb}"
    )
    lines = [
        f"! {title}",
        (" $CONTRL RUNTYP=ENERGY SCFTYP=ROHF ISPHER=1 UNITS=ANGS MAXIT=200 $END"
         if any(value != 1 for value in multiplicities)
         else " $CONTRL RUNTYP=ENERGY SCFTYP=RHF ISPHER=1 UNITS=ANGS MAXIT=200 $END"),
        " $SYSTEM MWORDS=500 $END",
        " $SCF DIRSCF=.TRUE. NPUNCH=0 CONV=1.0E-8 DIIS=.FALSE. SOSCF=.TRUE. $END",
        " $BASIS GBASIS=DFTB $END",
        dftb_line,
        *([" $DFTBSK"] + [
            f"   {first} {second} {dftb_parameter_dir}/{_dftb_symbol(first)}-{_dftb_symbol(second)}.skf"
            for first in species for second in species
        ] + [" $END"]),
        fmoprp,
        " CNVDMP=50 MCONV(2)=530 MCONV(4)=530 $END",
        *([" $DFT DC=.TRUE. IDCVER=4 $END",
           (" $PCM "
            f"SOLVNT={options.pcm_solvent} IEF={options.pcm_ief} "
            f"ICOMP={options.pcm_icomp} ICAV={options.pcm_icav} "
            f"IDISP={options.pcm_idisp} IFMO={options.pcm_ifmo} $END"),
           f" $PCMCAV RADII={options.pcm_cavity_radii} $END",
           f" $TESCAV NTSALL={options.pcm_tessellation} $END"]
          if options.use_pcm else []),
        (f" $FMO NFRAG={max_fragment} NBODY={min(2, max_fragment)} NACUT={options.nacut}"
         if full_energy else
         f" $FMO NFRAG={max_fragment} NBODY=2 NACUT={options.nacut} MODMOL=1 MOLFRG(1)={max_fragment}"),
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
