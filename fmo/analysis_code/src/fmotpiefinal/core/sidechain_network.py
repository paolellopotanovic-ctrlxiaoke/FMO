"""One-hop geometric network control for a frozen local sidechain trigger plan."""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json

import numpy as np

from .receptor_chemical_graph import CANONICAL
from .sidechain_policy import BACKBONE, canonical_name, key, residue_id


def extend_one_contact_shell(protein, seed_plan, *, n_hops=1):
    if not isinstance(n_hops, int) or isinstance(n_hops, bool) or n_hops not in {1, 2}:
        raise ValueError("registered geometric controls support exactly one or two contact hops")
    result = deepcopy(seed_plan)
    seeds = {tuple(residue) for residue in seed_plan["selected_residues"]}
    node_by_key = {tuple(node["residue"]): node for node in result["nodes"]}
    binding_region = {edge["receptor"] for edge in seed_plan["edges"] if edge["kind"] == "sidechain_ligand_contact"}
    points = {}
    for residue, node in node_by_key.items():
        if node["id"] not in binding_region or node["protection"] or canonical_name(residue[2]) not in CANONICAL - {"ALA", "GLY", "PRO"}:
            continue
        atoms = [atom.xyz for atom in protein.atoms if key(atom) == residue
                 and atom.record == "ATOM" and atom.element.upper() != "H" and atom.atom_name not in BACKBONE]
        if atoms:
            points[residue] = np.asarray(atoms)
    selected = set(seeds)
    cutoff = seed_plan["policy"]["contact_radius_angstrom"]
    frontier = set(seeds)
    for hop in range(1, n_hops+1):
        next_frontier = set()
        for seed in sorted(frontier):
            if seed not in points:
                raise ValueError("selected seed missing from eligible binding-sidechain graph")
            for candidate, positions in points.items():
                if candidate == seed:
                    continue
                distance = float(np.min(np.linalg.norm(points[seed][:, None, :] - positions[None, :, :], axis=2)))
                if distance > cutoff:
                    continue
                result["edges"].append({"kind": "sidechain_geometric_contact_one_hop" if n_hops == 1 else "sidechain_geometric_contact_bounded_hops",
                    "first": residue_id(seed), "second": residue_id(candidate), "distance_angstrom": distance,
                    "hop": hop, "confirmed_hydrogen_bond": False, "confirmed_covalent_bond": False})
                if candidate not in selected:
                    next_frontier.add(candidate)
                selected.add(candidate)
                node = node_by_key[candidate]
                node["action"] = "RESTRAINED_SIDECHAIN_RELAX"
                node["evidence"].append({"code": "ONE_HOP_CONTACT_TO_GEOMETRY_TRIGGER" if n_hops == 1 else "BOUNDED_CONTACT_NETWORK_TO_GEOMETRY_TRIGGER",
                    "seed": residue_id(seed), "distance_angstrom": distance, "hop": hop})
        frontier = next_frontier
    result["seed_plan_sha256"] = seed_plan["plan_sha256"]
    result["seed_residues"] = sorted(map(list, seeds))
    result["selected_residues"] = sorted(map(list, selected))
    result["schema"] = f"structure-driven-sidechain-{n_hops}-hop-control-v2"
    result["network_extension"] = {"n_hops": n_hops, "cutoff_angstrom": cutoff,
        "binding_region_only": True, "protect_special_networks": True,
        "recursive_expansion": n_hops > 1, "unbounded_expansion": False,
        "evidence": "geometric adjacency, not a confirmed interaction or rotamer model"}
    result.pop("plan_sha256")
    result["plan_sha256"] = hashlib.sha256(json.dumps(result, sort_keys=True, allow_nan=False).encode()).hexdigest()
    return result
