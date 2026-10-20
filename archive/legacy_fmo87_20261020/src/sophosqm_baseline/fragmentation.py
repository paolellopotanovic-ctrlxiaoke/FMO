"""FACIO-like residue fragmentation contract."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from collections import defaultdict
import json

from .structure import Atom, Ligand, Protein, residue_key
import math

def make_facio_like_fragments(protein: Protein, ligand: Ligand,
                              include_waters: bool = True, *,
                              reference_protein: Protein | None = None) -> tuple[Protein, list[Fragment], dict]:
    """Build the published FACIO-like protein fragment pattern.

    A fragment contains one residue's side chain, C-alpha and backbone N/H,
    plus the carbonyl C/O of the preceding residue.  The resulting cut is the
    C-alpha--carbonyl-C bond.  Selected residues are split into contiguous
    chain segments; each segment is saturated with a full neutral ACE or NME
    residue at its N- or C-terminal boundary.  Internal cuts are represented
    by GAMESS HOP_C, not by extra boundary hydrogen atoms.
    """
    groups: dict[str, list[Atom]] = {}
    reference_groups: dict[str, list[Atom]] = {}
    for atom in (reference_protein or protein).atoms:
        reference_groups.setdefault(residue_key(atom), []).append(atom)
    for atom in protein.atoms:
        residue = atom.residue_name.upper()
        if atom.record == "HETATM" and residue in {"HOH", "WAT", "TIP"}:
            if include_waters:
                groups.setdefault(residue_key(atom), []).append(atom)
            continue
        if atom.record == "ATOM":
            groups.setdefault(residue_key(atom), []).append(atom)

    def sort_key(key: str) -> tuple[str, int]:
        chain, _, number = key.split(":", 2)
        try:
            value = int(number)
        except ValueError:
            value = 10**9
        return chain, value

    keys = sorted(groups, key=sort_key)
    protein_keys = [key for key in keys if groups[key][0].residue_name.upper() not in {"HOH", "WAT", "TIP"}]
    chain_keys: dict[str, list[str]] = {}
    for key in protein_keys:
        chain_keys.setdefault(key.split(":", 1)[0], []).append(key)

    prepared_atoms = list(protein.atoms)
    next_index = max((atom.index for atom in prepared_atoms), default=-1) + 1
    cap_records: list[dict] = []
    fragments: list[Fragment] = []
    moved_carbonyls: dict[str, list[int]] = {key: [] for key in keys}
    hop_pairs: list[tuple[int, int]] = []
    for chain in chain_keys.values():
        positions = {key: index for index, key in enumerate(chain)}
        for current, following in zip(chain, chain[1:]):
            try:
                consecutive = int(following.split(":", 2)[2]) == int(current.split(":", 2)[2]) + 1
            except ValueError:
                consecutive = False
            if not consecutive:
                continue
            current_atoms = groups[current]
            for atom in current_atoms:
                if atom.atom_name.upper() in {"C", "O", "OXT"}:
                    moved_carbonyls[following].append(atom.index)
            ca = next((atom for atom in current_atoms if atom.atom_name.upper() == "CA"), None)
            carbonyl = next((atom for atom in current_atoms if atom.atom_name.upper() == "C"), None)
            if ca is not None and carbonyl is not None:
                # Actual GAMESS indices are assigned after fragment ordering;
                # keep source indices here and translate below.
                hop_pairs.append((ca.index, carbonyl.index))

    # Build contiguous selected segments for terminal capping.
    segment_starts: set[str] = set()
    segment_ends: set[str] = set()
    for chain in chain_keys.values():
        for index, key in enumerate(chain):
            previous = chain[index - 1] if index else None
            following = chain[index + 1] if index + 1 < len(chain) else None
            def adjacent(left: str | None, right: str | None) -> bool:
                if left is None or right is None:
                    return False
                try:
                    return int(right.split(":", 2)[2]) == int(left.split(":", 2)[2]) + 1
                except ValueError:
                    return False
            if not adjacent(previous, key):
                segment_starts.add(key)
            if not adjacent(key, following):
                segment_ends.add(key)

    def normalized(vector: tuple[float, float, float]) -> tuple[float, float, float]:
        norm = math.sqrt(sum(value * value for value in vector)) or 1.0
        return tuple(value / norm for value in vector)

    def subtract(left: Atom | tuple[float, float, float], right: Atom | tuple[float, float, float]) -> tuple[float, float, float]:
        left_xyz = left.xyz if isinstance(left, Atom) else left
        right_xyz = right.xyz if isinstance(right, Atom) else right
        return tuple(left_xyz[i] - right_xyz[i] for i in range(3))

    def add(xyz: tuple[float, float, float], vector: tuple[float, float, float]) -> tuple[float, float, float]:
        return tuple(xyz[i] + vector[i] for i in range(3))

    def scale(vector: tuple[float, float, float], factor: float) -> tuple[float, float, float]:
        return tuple(factor * value for value in vector)

    def cross(left: tuple[float, float, float], right: tuple[float, float, float]) -> tuple[float, float, float]:
        return (left[1] * right[2] - left[2] * right[1],
                left[2] * right[0] - left[0] * right[2],
                left[0] * right[1] - left[1] * right[0])

    def adjacent_residue_atoms(anchor: Atom, offset: int) -> list[Atom]:
        try:
            number = int(anchor.residue_number) + offset
        except ValueError:
            return []
        matches = [atoms for key, atoms in reference_groups.items()
                   if key.split(":", 1)[0] == anchor.chain and key.rsplit(":", 1)[-1] == str(number)]
        return min(matches, key=lambda atoms: math.dist(anchor.xyz, next(
            (atom.xyz for atom in atoms if atom.element.upper() == "N" if offset > 0), atoms[0].xyz))
            if matches else float("inf")) if matches else []

    def new_atom(element: str, atom_name: str, xyz: tuple[float, float, float],
                 chain: str, residue_number: str, residue_name: str) -> Atom:
        nonlocal next_index
        atom = Atom(next_index, element.upper(), xyz, chain, residue_name,
                    residue_number, "HETATM", atom_name)
        next_index += 1
        return atom

    def methyl_hydrogens(carbon_xyz: tuple[float, float, float], parent_xyz: tuple[float, float, float],
                         residue_number: str, residue_name: str, chain: str) -> list[Atom]:
        parent_direction = normalized(subtract(parent_xyz, carbon_xyz))
        arbitrary = (1.0, 0.0, 0.0) if abs(parent_direction[0]) < 0.8 else (0.0, 1.0, 0.0)
        side = normalized(cross(parent_direction, arbitrary))
        up = normalized(cross(parent_direction, side))
        atoms = []
        for index, (side_factor, up_factor) in enumerate(((-0.866, 0.5), (0.866, 0.5), (0.0, -1.0)), 1):
            direction = normalized(tuple(
                -parent_direction[i] / 3.0
                + 0.9428090415820634 * (side_factor * side[i] + up_factor * up[i])
                for i in range(3)
            ))
            atoms.append(new_atom("H", f"HH3{index}", add(carbon_xyz, scale(direction, 1.09)), chain,
                                  residue_number, residue_name))
        return atoms

    def add_ace(n_atom: Atom, ca_atom: Atom) -> list[Atom]:
        nonlocal next_index
        residue_number = f"{10_000 + len(cap_records)}"
        predecessor = adjacent_residue_atoms(n_atom, -1)
        predecessor_ca = next((atom for atom in predecessor if atom.atom_name.upper() == "CA"), None)
        predecessor_carbonyl = next((atom for atom in predecessor if atom.atom_name.upper() == "C"), None)
        predecessor_oxygen = next((atom for atom in predecessor if atom.atom_name.upper() == "O"), None)
        if predecessor_ca is not None and predecessor_carbonyl is not None and predecessor_oxygen is not None:
            carbonyl = predecessor_carbonyl.xyz
            oxygen = predecessor_oxygen.xyz
            methyl = predecessor_ca.xyz
            reference_indices = [predecessor_ca.index, predecessor_carbonyl.index, predecessor_oxygen.index]
        else:
            carbonyl = add(n_atom.xyz, scale(normalized(subtract(n_atom.xyz, ca_atom.xyz)), 1.335))
            planar_side = normalized(cross(subtract(n_atom.xyz, ca_atom.xyz), subtract(carbonyl, ca_atom.xyz)))
            oxygen = add(carbonyl, scale(normalized(tuple(
                -0.5 * normalized(subtract(n_atom.xyz, carbonyl))[i] + 0.866 * planar_side[i] for i in range(3))), 1.229))
            methyl = add(carbonyl, scale(normalized(tuple(
                -0.5 * normalized(subtract(n_atom.xyz, carbonyl))[i] - 0.866 * planar_side[i] for i in range(3))), 1.520))
            reference_indices = []
        atoms = [new_atom("C", "C", carbonyl, n_atom.chain, residue_number, "ACE"),
                 new_atom("O", "O", oxygen, n_atom.chain, residue_number, "ACE"),
                 new_atom("C", "CH3", methyl, n_atom.chain, residue_number, "ACE")]
        atoms.extend(methyl_hydrogens(methyl, carbonyl, residue_number, "ACE", n_atom.chain))
        prepared_atoms.extend(atoms)
        cap_records.append({"type": "ACE", "anchor_atom": n_atom.index,
                            "residue_number": residue_number, "atom_indices": [atom.index for atom in atoms],
                            "retained_backbone_n_hydrogen": True,
                            "reference_atom_indices": reference_indices})
        return atoms

    def add_nme(c_atom: Atom, ca_atom: Atom) -> list[Atom]:
        nonlocal next_index
        residue_number = f"{20_000 + len(cap_records)}"
        to_ca = normalized(subtract(ca_atom.xyz, c_atom.xyz))
        oxygen = next((atom for atom in groups[residue_key(c_atom)] if atom.atom_name.upper() == "O"), None)
        if oxygen is None:
            raise ValueError(f"selected terminal residue lacks carbonyl oxygen: {residue_key(c_atom)}")
        to_oxygen = normalized(subtract(oxygen.xyz, c_atom.xyz))
        successor = adjacent_residue_atoms(c_atom, 1)
        successor_ca = next((atom for atom in successor if atom.atom_name.upper() == "CA"), None)
        successor_n = next((atom for atom in successor if atom.atom_name.upper() == "N"), None)
        successor_h = next((atom for atom in successor if atom.element.upper() == "H" and
                            math.dist(atom.xyz, successor_n.xyz) <= 1.25), None) if successor_n is not None else None
        if successor_ca is not None and successor_n is not None and successor_h is not None:
            nitrogen = successor_n.xyz
            hydrogen = successor_h.xyz
            methyl = successor_ca.xyz
            reference_indices = [successor_ca.index, successor_n.index, successor_h.index]
        else:
            nitrogen_direction = normalized(tuple(-(to_ca[i] + to_oxygen[i]) for i in range(3)))
            nitrogen = add(c_atom.xyz, scale(nitrogen_direction, 1.335))
            side = normalized(cross(nitrogen_direction, to_ca))
            methyl_direction = normalized(tuple(nitrogen_direction[i] + 0.8 * side[i] for i in range(3)))
            methyl = add(nitrogen, scale(methyl_direction, 1.450))
            hydrogen = add(nitrogen, scale(normalized(tuple(
                -nitrogen_direction[i] + 0.35 * side[i] for i in range(3))), 1.010))
            reference_indices = []
        atoms = [new_atom("N", "N", nitrogen, c_atom.chain, residue_number, "NME"),
                 new_atom("H", "H", hydrogen, c_atom.chain, residue_number, "NME"),
                 new_atom("C", "CH3", methyl, c_atom.chain, residue_number, "NME")]
        atoms.extend(methyl_hydrogens(methyl, nitrogen, residue_number, "NME", c_atom.chain))
        prepared_atoms.extend(atoms)
        cap_records.append({"type": "NME", "anchor_atom": c_atom.index,
                            "residue_number": residue_number, "atom_indices": [atom.index for atom in atoms],
                            "reference_atom_indices": reference_indices})
        return atoms

    cap_groups: dict[str, list[Atom]] = defaultdict(list)
    for atom in prepared_atoms:
        if atom.residue_name.upper() in {"ACE", "NME"}:
            cap_groups[f"{atom.chain}:{atom.residue_name}:{atom.residue_number}"].append(atom)
    cap_indices: dict[str, list[int]] = {}
    nme_replaced_oxt: dict[str, list[int]] = {}
    for key, atoms in cap_groups.items():
        cap_type = atoms[0].residue_name.upper()
        anchor_name = "N" if cap_type == "ACE" else "C"
        anchors = [atom for atom in prepared_atoms if atom.record == "ATOM" and atom.atom_name.upper() == anchor_name]
        anchor = min(anchors, key=lambda atom: math.dist(atom.xyz, atoms[0].xyz))
        cap_indices[residue_key(anchor)] = [atom.index for atom in atoms]
        cap_records.append({"type": cap_type, "anchor_atom": anchor.index,
                            "residue_number": atoms[0].residue_number,
                            "atom_indices": [atom.index for atom in atoms], "pre_existing": True})

    for fragment_id, key in enumerate(keys, start=1):
        atoms = groups[key]
        residue = atoms[0].residue_name.upper()
        water = residue in {"HOH", "WAT", "TIP"}
        moved = {index for values in moved_carbonyls.values() for index in values}
        oxt_replacements = nme_replaced_oxt.setdefault(key, [])
        if key in segment_ends:
            oxt_replacements.extend(
                atom.index for atom in atoms
                if atom.atom_name.upper() == "OXT" and atom.index not in oxt_replacements
            )
        replaced_oxt = set(oxt_replacements)
        atom_indices = [atom.index for atom in atoms
                        if atom.index not in moved and atom.index not in replaced_oxt]
        atom_indices.extend(moved_carbonyls[key])
        atom_indices.extend(cap_indices.get(key, []))
        if not water and key in segment_starts:
            n_atom = next((atom for atom in atoms if atom.atom_name.upper() == "N"), None)
            ca_atom = next((atom for atom in atoms if atom.atom_name.upper() == "CA"), None)
            if residue_key(n_atom) not in cap_indices:
                atom_indices.extend(atom.index for atom in add_ace(n_atom, ca_atom))
        if not water and key in segment_ends:
            c_atom = next((atom for atom in atoms if atom.atom_name.upper() == "C"), None)
            ca_atom = next((atom for atom in atoms if atom.atom_name.upper() == "CA"), None)
            if residue_key(c_atom) not in cap_indices:
                atom_indices.extend(atom.index for atom in add_nme(c_atom, ca_atom))
        charge = 0 if water else _prepared_residue_charge(residue, atoms)
        fragments.append(Fragment(fragment_id, key, "water" if water else "protein_residue",
                                  tuple(atom_indices), charge=charge))

    fragments.append(Fragment(len(fragments) + 1, ligand.name, "ligand",
                              tuple(atom.index for atom in ligand.atoms), charge=ligand.formal_charge))
    source_to_gamess: dict[int, int] = {}
    current_index = 1
    for fragment in fragments:
        for source_index in fragment.atom_indices:
            source_to_gamess[source_index] = current_index
            current_index += 1
    metadata = {
        "protocol": "SophosQM-like FACIO approximation",
        "fragment_pattern": "side chain + CA + backbone N/H + preceding carbonyl C/O",
        "cut_bond": "C-alpha to carbonyl C prime",
        "terminal_caps": "full neutral ACE/NME residues at selected segment ends",
        "hop_pairs_source_indices": [[-left, right] for left, right in hop_pairs],
        "cap_records": cap_records,
        "nme_replaced_oxt_atom_indices": nme_replaced_oxt,
        "protonation": "input PDB hydrogens retained; preparation backend recorded separately",
    }
    metadata["hop_pairs_gamess_indices"] = [[-source_to_gamess[left], source_to_gamess[right]]
                                              for left, right in hop_pairs
                                              if left in source_to_gamess and right in source_to_gamess]
    replaced_oxt_indices = {
        index for indices in nme_replaced_oxt.values() for index in indices
    }
    retained_indices = {
        atom_index
        for fragment in fragments
        for atom_index in fragment.atom_indices
    }
    capped_protein = Protein(tuple(
        atom for atom in prepared_atoms
        if atom.index in retained_indices and atom.index not in replaced_oxt_indices
    ))
    return capped_protein, fragments, metadata


def _prepared_residue_charge(residue: str, atoms: list[Atom]) -> int:
    if residue == "HIS":
        protonated_nitrogens = {
            atom.atom_name.upper() for atom in atoms if atom.atom_name.upper() in {"HD1", "HE2"}
        }
        if protonated_nitrogens == {"HD1", "HE2"}:
            return 1
        if protonated_nitrogens:
            return 0
        raise ValueError(f"ambiguous unprotonated histidine {residue_key(atoms[0])}")
    return {
        "ASP": -1, "GLU": -1, "ASH": 0, "GLH": 0,
        "LYS": 1, "ARG": 1, "HIP": 1, "HID": 0, "HIE": 0,
    }.get(residue, 0)


@dataclass(frozen=True)
class Fragment:
    fragment_id: int
    name: str
    kind: str
    atom_indices: tuple[int, ...]
    charge: int = 0
    multiplicity: int = 1


def write_fragment_manifest(path: str | Path, fragments: list[Fragment], metadata: dict) -> None:
    payload = {"metadata": metadata, "fragments": [asdict(fragment) for fragment in fragments]}
    Path(path).write_text(json.dumps(payload, indent=2), encoding="utf-8")
