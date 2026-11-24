"""Label-free structure-repair decisions; local strain never selects a score."""
from __future__ import annotations

import numpy as np

from .conjugated_amine_geometry import induced_conjugated_nh2, reconstruct_flagged_nh2


SCHEMA = "fmotpiefinal-local-strain-repair-audit-v1"
HARD_GEOMETRY_STATES = {"REPAIR_REQUIRED", "INCOMPLETE", "REVIEW_REQUIRED"}
GEOMETRY_STATES = HARD_GEOMETRY_STATES | {"PASS", "REVIEW"}


def audit_local_strain_geometry(molecule, reference_coordinates, bound_coordinates, *,
                                geometry_status: str = "PASS",
                                local_strain_kcal_mol: float | None = None) -> dict:
    reference = np.asarray(reference_coordinates, dtype=float)
    bound = np.asarray(bound_coordinates, dtype=float)
    expected_shape = (molecule.GetNumAtoms(), 3)
    if (reference.shape != expected_shape or bound.shape != expected_shape or
            not np.isfinite(reference).all() or not np.isfinite(bound).all()):
        raise ValueError("finite, atom-matched reference and bound geometries required")
    if geometry_status not in GEOMETRY_STATES:
        raise ValueError("unknown geometry audit status")

    nh2 = induced_conjugated_nh2(molecule, reference, bound)
    new_flags = [row for row in nh2 if row["new_deformation_flag"]]
    input_flags = [row for row in nh2 if row["already_flagged_input"]]
    findings = []
    if geometry_status in HARD_GEOMETRY_STATES:
        findings.append({"rule": "hard_geometry_gate", "evidence": geometry_status,
                         "action": "STOP_AND_REPAIR_OR_REVIEW"})
    elif geometry_status == "REVIEW":
        findings.append({"rule": "geometry_warning_contacts", "evidence": geometry_status,
                         "action": "RECORD_WARNING_NO_AUTOMATIC_REPAIR"})
    if new_flags:
        findings.append({"rule": "k10_induced_conjugated_nh2_deformation",
                         "evidence": new_flags,
                         "action": "REGISTERED_EXPERIMENTAL_REPAIR_CANDIDATE"})
    if input_flags:
        findings.append({"rule": "input_conjugated_nh2_already_distorted",
                         "evidence": input_flags, "action": "REVIEW_INPUT_GEOMETRY"})

    strain = None
    if local_strain_kcal_mol is not None:
        value = float(local_strain_kcal_mol)
        if not np.isfinite(value):
            raise ValueError("local strain must be finite")
        strain = {"value_kcal_mol": value, "role": "DIAGNOSTIC_ONLY",
                  "negative_value_review": value < 0.0,
                  "automatic_geometry_repair_trigger": False,
                  "automatic_score_term": False}
        if value < 0.0:
            findings.append({"rule": "negative_local_strain", "evidence": value,
                             "action": "REVIEW_ENERGY_REFERENCE_AND_QM_MODEL_ONLY"})

    if geometry_status in HARD_GEOMETRY_STATES:
        decision = "REPAIR_REQUIRED_OR_REVIEW"
    elif (new_flags or input_flags or (strain and strain["negative_value_review"])):
        decision = "REVIEW_REQUIRED_NO_AUTOMATIC_SCORE_CHANGE"
    else:
        decision = "NO_REPAIR_TRIGGER_TPIE_REMAINS_PRIMARY"
    return {
        "schema": SCHEMA,
        "decision": decision,
        "geometry_status": geometry_status,
        "findings": findings,
        "conjugated_nh2": nh2,
        "local_strain": strain,
        "labels_used_for_trigger": False,
        "coordinates_modified": False,
        "ranking_modified": False,
        "primary_score": "pure_tPIE",
        "local_strain_score_enabled": False,
        "repair_candidates": [
            {"rule": "conjugated_nh2_planar_hydrogen_rebuild_v1",
             "groups": new_flags,
             "status": "EXPERIMENTAL_ONLY_REQUIRES_FULL_K10_AND_QM_RERUN"}
        ] if new_flags else [],
    }


def apply_experimental_nh2_repair(molecule, bound_coordinates, audit: dict) -> tuple[np.ndarray, dict]:
    if audit.get("schema") != SCHEMA or audit.get("labels_used_for_trigger") is not False:
        raise ValueError("verified label-free repair audit required")
    if audit.get("decision") != "REVIEW_REQUIRED_NO_AUTOMATIC_SCORE_CHANGE":
        raise ValueError("audit does not authorize an experimental NH2 candidate")
    if not audit.get("repair_candidates"):
        raise ValueError("no registered conjugated-NH2 repair candidate")
    before = np.asarray(bound_coordinates, dtype=float)
    if before.shape != (molecule.GetNumAtoms(), 3) or not np.isfinite(before).all():
        raise ValueError("finite atom-matched bound coordinates required")
    after = reconstruct_flagged_nh2(molecule, before, audit["conjugated_nh2"])
    hydrogen_indices = {
        index for row in audit["conjugated_nh2"] if row["new_deformation_flag"]
        for index in row["hydrogens"]
    }
    heavy_indices = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1]
    changed = {index for index in range(molecule.GetNumAtoms())
               if not np.array_equal(before[index], after[index])}
    if changed - hydrogen_indices or not np.array_equal(before[heavy_indices], after[heavy_indices]):
        raise ValueError("candidate repair changed atoms outside flagged hydrogens")
    if not np.isfinite(after).all():
        raise ValueError("candidate repair produced nonfinite coordinates")
    report = {"rule": "conjugated_nh2_planar_hydrogen_rebuild_v1",
              "changed_atom_indices_zero_based": sorted(changed),
              "changed_atom_elements": [molecule.GetAtomWithIdx(index).GetSymbol()
                                        for index in sorted(changed)],
              "heavy_atom_coordinates_unchanged": True,
              "bond_graph_or_formal_charge_modified": False,
              "requires_full_k10_and_qm_rerun": True,
              "production_default_changed": False}
    return after, report
