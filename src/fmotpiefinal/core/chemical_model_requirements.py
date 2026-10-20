"""Generate per-component model requirements without guessing missing chemistry."""
from pathlib import Path

import numpy as np

from .component_inventory import METALS


def model_requirements(graph, parameter_dir=None):
    from .special_dispatch import choose_route
    components = [row for row in graph["residues"] if row["kind"] not in {"standard_residue", "water"}]
    all_atoms = [(row, atom) for row in graph["residues"] for atom in row["atoms"]]
    elements = sorted({atom["element"] for _, atom in all_atoms})
    missing_pairs = None
    if parameter_dir is not None:
        parameter_dir = Path(parameter_dir)
        if not parameter_dir.is_dir():
            raise FileNotFoundError(parameter_dir)
        missing_pairs = [f"{first.title()}-{second.title()}.skf" for first in elements for second in elements
                         if not (parameter_dir / f"{first.title()}-{second.title()}.skf").is_file()]
    rows = []
    for component in components:
        try:
            route = choose_route({**graph, "residues": [component]})
            reason = None
        except ValueError as error:
            route, reason = None, str(error)
        metal_atoms = [atom for atom in component["atoms"] if atom["element"] in METALS]
        contacts = []
        for metal in metal_atoms:
            for other_residue, other in all_atoms:
                if other_residue["model"] != component["model"] or other_residue["id"] == component["id"]:
                    continue
                if other["element"] not in {"N", "O", "S"}:
                    continue
                distance = float(np.linalg.norm(np.asarray(metal["xyz"]) - np.asarray(other["xyz"])))
                if distance <= 3.5:
                    contacts.append({"metal_atom": metal["name"], "candidate": other_residue["id"] + "/" + other["name"],
                                     "distance_angstrom": distance, "bond_assigned": False})
        heavy_names = [atom["name"] for atom in component["atoms"] if atom["element"] not in {"H", "D"}]
        requirements = ["source_atom_mapping_and_missing_atom_check", "explicit_formal_charge_and_protonation_contract",
                        "complete_MM_topology_and_charge_validation", "QM_region_and_fragment_boundary_validation",
                        "ordered_SKF_pair_and_final_convergence_validation", "matched_full_cohort_regression"]
        if metal_atoms:
            requirements += ["metal_oxidation_and_spin_if_QM", "coordination_model_with_donor_protonation",
                             "do_not_infer_covalent_metal_bonds_from_distance_only"]
        else:
            requirements += ["pinned_bond_orders_stereochemistry_and_covalent_attachment_check"]
        rows.append({"component": component["id"], "model": component["model"], "kind": component["kind"],
            "elements": component["elements"], "individual_mm_route": route, "route_failure": reason,
            "source_explicit_hydrogens": sum(atom["element"] in {"H", "D"} for atom in component["atoms"]),
            "duplicate_heavy_atom_names": sorted({name for name in heavy_names if heavy_names.count(name) > 1}),
            "candidate_metal_contacts": sorted(contacts, key=lambda row: row["distance_angstrom"]),
            "required_validation": requirements, "automatic_qm_approval": False})
    return {"schema": "chemical-model-requirements-v1", "source_sha256": graph["source_sha256"],
        "experimental_labels_used": False, "components": rows, "observed_elements": elements,
        "missing_ordered_skf_pairs_for_observed_source_elements": missing_pairs,
        "parameter_scope": "source inventory only; actual ligand-inclusive QM deck must be checked separately",
        "distance_contacts_are_not_bonds": True, "joint_component_route_not_implied_by_individual_routes": True,
        "unresolved_covalent_adduct_detection": "requires explicit LINK/CONECT/topology and ligand chemistry; proximity is not proof",
        "ranking_eligible": False}
