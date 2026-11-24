"""Label-free, conservative selection of receptor sidechain degrees of freedom."""
from __future__ import annotations

from collections import defaultdict
from dataclasses import asdict, dataclass
import hashlib
import json
import math

import numpy as np
from scipy.spatial import cKDTree

from .receptor_chemical_graph import CANONICAL, VARIANTS, METALS, COORDINATION_LIMITS


BACKBONE = {"N", "CA", "C", "O", "OXT", "H", "H1", "H2", "H3", "HA", "HA2", "HA3"}
RADII = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80,
         "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98}


@dataclass(frozen=True)
class SidechainPolicy:
    contact_radius_angstrom: float = 4.0
    amide_contact_radius_angstrom: float = 3.5
    heavy_overlap_angstrom: float = 0.6
    polar_pair_min_distance_angstrom: float = 2.2
    amide_departure_degrees: float = 30.0
    amide_departure_increase_degrees: float = 20.0
    disulfide_distance_angstrom: float = 2.3

    def __post_init__(self):
        if any(not math.isfinite(value) or value <= 0 for value in asdict(self).values()):
            raise ValueError("sidechain policy thresholds must be finite and positive")


def key(atom):
    return atom.chain, atom.residue_number, atom.residue_name


def residue_id(residue):
    return f"{residue[0] or '_'}:{residue[2]}:{residue[1]}"


def canonical_name(residue):
    return VARIANTS.get(residue, residue)


def amide_changes(molecule, original, bound, policy):
    from rdkit import Chem
    from rdkit.Chem import rdMolTransforms
    conformers = []
    for positions in (original, bound):
        conformer = Chem.Conformer(molecule.GetNumAtoms())
        for index, position in enumerate(positions):
            conformer.SetAtomPosition(index, tuple(map(float, position)))
        conformers.append(conformer)
    records = []
    for carbon in molecule.GetAtoms():
        if carbon.GetAtomicNum() != 6:
            continue
        oxygens = [atom for atom in carbon.GetNeighbors() if atom.GetAtomicNum() == 8
                   and molecule.GetBondBetweenAtoms(carbon.GetIdx(), atom.GetIdx()).GetBondType() == Chem.BondType.DOUBLE]
        nitrogens = [atom for atom in carbon.GetNeighbors() if atom.GetAtomicNum() == 7]
        for oxygen in oxygens:
            for nitrogen in nitrogens:
                for substituent in nitrogen.GetNeighbors():
                    if substituent.GetIdx() == carbon.GetIdx():
                        continue
                    indices = [oxygen.GetIdx(), carbon.GetIdx(), nitrogen.GetIdx(), substituent.GetIdx()]
                    departures = []
                    for conformer in conformers:
                        angle = float(rdMolTransforms.GetDihedralDeg(conformer, *indices))
                        if not math.isfinite(angle):
                            raise ValueError("undefined amide dihedral; diagnose overlapping or collinear atoms")
                        departures.append(min(abs(angle), abs(180-abs(angle))))
                    triggered = departures[1] >= policy.amide_departure_degrees and departures[1]-departures[0] >= policy.amide_departure_increase_degrees
                    records.append({"atoms_zero_based": indices, "amide_atoms_zero_based": indices[:3],
                        "original_departure_deg": departures[0], "bound_departure_deg": departures[1],
                        "new_deformation_flag": triggered,
                        "interpretation": "geometry screen, not a rotamer probability or strain energy"})
    return records


def select_sidechains(protein, molecules, bound_positions, *, policy=None, covalent_links=()):
    policy = policy or SidechainPolicy()
    if not molecules or set(molecules) != set(bound_positions):
        raise ValueError("complete matching original/bound ligand identities required")
    groups = defaultdict(list)
    for atom in protein.atoms:
        if not np.isfinite(atom.xyz).all():
            raise ValueError("nonfinite receptor coordinates")
        groups[key(atom)].append(atom)
    protected = defaultdict(list)
    protection_edges = []
    sulfurs = [atom for atom in protein.atoms if atom.atom_name == "SG" and canonical_name(atom.residue_name) == "CYS"]
    for first_index, first in enumerate(sulfurs):
        for second in sulfurs[first_index+1:]:
            distance = math.dist(first.xyz, second.xyz)
            if distance <= policy.disulfide_distance_angstrom:
                for atom in (first, second):
                    protected[key(atom)].append("DISULFIDE_NETWORK")
                protection_edges.append({"kind": "disulfide_distance_candidate", "first": residue_id(key(first)),
                                         "second": residue_id(key(second)), "distance_angstrom": distance})
    metals = [atom for atom in protein.atoms if atom.element.upper() in METALS]
    for metal in metals:
        for atom in protein.atoms:
            if atom.element.upper() not in {"N", "O", "S"} or atom.record != "ATOM":
                continue
            distance = math.dist(atom.xyz, metal.xyz)
            if distance <= COORDINATION_LIMITS.get(metal.element.upper(), 3.0):
                protected[key(atom)].append("METAL_COORDINATION_CANDIDATE")
                protection_edges.append({"kind": "metal_distance_candidate", "first": residue_id(key(metal)),
                    "second": residue_id(key(atom)), "atom": atom.atom_name, "distance_angstrom": distance,
                    "confirmed_coordination_bond": False})
    for first, second in covalent_links:
        if first not in groups or second not in groups:
            raise ValueError("covalent link references absent residue")
        protected[first].append("EXPLICIT_NONPEPTIDE_COVALENT_LINK")
        protected[second].append("EXPLICIT_NONPEPTIDE_COVALENT_LINK")
        protection_edges.append({"kind": "protected_covalent_link", "first": residue_id(first), "second": residue_id(second)})
    receptor_sidechains = []
    for residue, atoms in groups.items():
        base = canonical_name(residue[2])
        if base not in CANONICAL or base in {"GLY", "ALA", "PRO"}:
            if base == "PRO":
                protected[residue].append("RING_COUPLED_TO_FIXED_BACKBONE")
            continue
        receptor_sidechains.extend(atom for atom in atoms if atom.record == "ATOM"
            and atom.atom_name not in BACKBONE and atom.element.upper() != "H")
    if not receptor_sidechains:
        sidechain_points = np.empty((0, 3))
    else:
        sidechain_points = np.asarray([atom.xyz for atom in receptor_sidechains])
    contact_edges, triggers, deformation = [], defaultdict(list), {}
    ligand_inventory = {}
    for name in sorted(molecules):
        molecule = molecules[name]
        original = np.asarray(molecule.GetConformer().GetPositions(), dtype=float)
        bound = np.asarray(bound_positions[name], dtype=float)
        if original.shape != bound.shape or bound.shape != (molecule.GetNumAtoms(), 3):
            raise ValueError(f"ligand atom count/order mismatch: {name}")
        if not np.isfinite(original).all() or not np.isfinite(bound).all():
            raise ValueError("nonfinite ligand coordinates")
        heavy_indices = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1]
        if not heavy_indices:
            raise ValueError("ligand contains no heavy atoms")
        ligand_inventory[name] = {"n_atoms": molecule.GetNumAtoms(), "formal_charge": sum(atom.GetFormalCharge() for atom in molecule.GetAtoms())}
        deformation[name] = amide_changes(molecule, original, bound, policy)
        changed = [record for record in deformation[name] if record["new_deformation_flag"]]
        for pose_type, positions in (("original", original), ("bound", bound)):
            tree = cKDTree(positions[heavy_indices])
            nearby = tree.query_ball_point(sidechain_points, policy.contact_radius_angstrom)
            minimum = {}
            for receptor_index, neighbors in enumerate(nearby):
                atom = receptor_sidechains[receptor_index]
                residue = key(atom)
                for neighbor in neighbors:
                    ligand_index = heavy_indices[neighbor]
                    element = molecule.GetAtomWithIdx(ligand_index).GetSymbol().upper()
                    distance = float(np.linalg.norm(np.asarray(atom.xyz)-positions[ligand_index]))
                    event = {"ligand": name, "pose": pose_type, "receptor": residue_id(residue),
                        "receptor_atom": atom.atom_name, "ligand_atom_zero_based": ligand_index,
                        "distance_angstrom": distance}
                    if residue not in minimum or distance < minimum[residue]["distance_angstrom"]:
                        minimum[residue] = event
                    if atom.element.upper() not in RADII or element not in RADII:
                        continue
                    overlap = RADII[atom.element.upper()] + RADII[element] - distance
                    polar_pair = atom.element.upper() in {"N", "O", "S"} and element in {"N", "O", "S"}
                    severe = overlap >= policy.heavy_overlap_angstrom and (not polar_pair or distance < policy.polar_pair_min_distance_angstrom)
                    if severe and pose_type == "bound":
                        triggers[residue].append({**event, "code": "PERSISTENT_BOUND_HEAVY_OVERLAP", "overlap_angstrom": overlap})
            contact_edges.extend({"kind": "sidechain_ligand_contact", **event} for event in minimum.values())
        for change in changed:
            for receptor_index, point in enumerate(sidechain_points):
                distances = np.linalg.norm(bound[change["amide_atoms_zero_based"]] - point, axis=1)
                distance = float(min(distances))
                if distance <= policy.amide_contact_radius_angstrom:
                    atom = receptor_sidechains[receptor_index]
                    triggers[key(atom)].append({"code": "NEW_LIGAND_AMIDE_DEFORMATION_NEAR_SIDECHAIN",
                        "ligand": name, "receptor_atom": atom.atom_name, "distance_angstrom": distance,
                        "amide": change})
    nodes = []
    selected = []
    for residue in sorted(groups):
        reasons = sorted(set(protected.get(residue, [])))
        evidence = triggers.get(residue, [])
        if evidence and not reasons:
            action = "RESTRAINED_SIDECHAIN_RELAX"
            selected.append(list(residue))
        elif evidence:
            action = "REVIEW_PROTECTED_NETWORK_KEEP_FIXED"
        else:
            action = "KEEP_FIXED_NO_TRIGGER"
        nodes.append({"id": residue_id(residue), "residue": list(residue),
            "n_atoms": len(groups[residue]), "canonical_name": canonical_name(residue[2]),
            "protection": reasons, "action": action, "evidence": evidence})
    payload = {"schema": "structure-driven-sidechain-policy-v1", "policy": asdict(policy),
        "selected_residues": selected, "nodes": nodes, "edges": protection_edges + contact_edges,
        "ligand_deformation": deformation, "ligand_inventory": ligand_inventory,
        "experimental_labels_used": False, "score_values_used_for_selection": False,
        "action": "RESTRAINED_SIDECHAIN_RELAX" if selected else "KEEP_BASELINE_UNCHANGED",
        "scope": "experimental local geometry rule, not a universal chemistry repair",
        "water_hydrogen_relaxation_enabled": False, "new_protonation_states_generated": False,
        "rotamer_library_used": False, "unknown_chemistry_parameterized": False,
        "production_default_changed": False,
        "future_ligands": "freeze and reuse this target plan; new out-of-region chemistry requires an explicit plan revision"}
    payload["plan_sha256"] = hashlib.sha256(json.dumps(payload, sort_keys=True, allow_nan=False).encode()).hexdigest()
    return payload


def explicit_nonpeptide_links(path):
    from openmm.app import PDBFile
    topology = PDBFile(str(path)).topology
    result = set()
    for first, second in topology.bonds():
        if first.residue == second.residue:
            continue
        peptide = {first.name, second.name} == {"C", "N"} and all(
            canonical_name(atom.residue.name) in CANONICAL | {"ACE", "NME"} for atom in (first, second))
        if peptide:
            continue
        residues = [(atom.residue.chain.id, atom.residue.id + atom.residue.insertionCode.strip(), atom.residue.name)
                    for atom in (first, second)]
        result.add(tuple(sorted(residues)))
    return sorted(result)
