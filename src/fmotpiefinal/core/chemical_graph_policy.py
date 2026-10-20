"""Packaged label-free chemistry diagnostics, not affinity-based route selection."""
import hashlib
import json


def plan_actions(graph, retention=None):
    tasks = []

    def add(code, severity, evidence, action):
        tasks.append({"code": code, "severity": severity, "evidence": evidence, "action": action,
            "automatic_coordinate_edit": False, "ranking_benefit_demonstrated": False})

    for residue in graph["residues"]:
        kind = residue["kind"]
        template = residue.get("forcefield_template", {})
        if kind in {"canonical_protein", "water"} and template.get("status") != "MATCHED":
            add("INCOMPLETE_OR_UNSUPPORTED_STANDARD_TOPOLOGY", "BLOCK", residue,
                "Diagnose missing heavy atoms/H, terminal state and external bonds; run registered structure/protonation repair then re-audit. Do not mutate or delete to obtain a match.")
        if template.get("matching_ignores_external_bonds") is True:
            add("EXTERNAL_BOND_TEMPLATE_NOT_VALIDATED", "BLOCK", residue,
                "Require complete-topology template validation; isolated residue matching does not validate peptide/disulfide connections.")
        if kind == "nonstandard_component":
            add("NONSTANDARD_COMPONENT_REQUIRES_MODEL", "BLOCK", residue,
                "Resolve pinned CCD identity, missing atoms, protonation and covalent links; retain component and provide validated MD/fragment/QM models. CCD annotation is not parameterization.")
        if kind == "other_ion":
            add("ION_RETENTION_AND_CHARGE_CONTRACT", "BLOCK", residue,
                "Register whether this ion is retained in the target receptor; validate charge and fragment support before scoring.")
    for edge in graph["cross_residue_bonds"]:
        if edge["geometry_supported"] is False:
            add("INVALID_COVALENT_LINK_GEOMETRY", "BLOCK", edge,
                "Check chain breaks, missing sequence and CONECT provenance before loop/terminal rebuilding. Never bridge a long gap by renumbering alone.")
        elif edge["kind"] == "explicit_covalent_link":
            add("COVALENT_ADDUCT_MODEL_REQUIRED", "BLOCK", edge,
                "Require a registered covalent-adduct topology and fragment boundary model.")
    for atom in graph["unresolved_hydrogen_parent_atoms"]:
        add("UNRESOLVED_HYDROGEN_PARENT", "BLOCK", atom,
            "Resolve parent identity and protonation before constrained H relaxation; proximity alone may be ambiguous.")
    for edge in graph.get("protein_bond_geometry_screens", []):
        add("PROTEIN_BOND_LENGTH_SCREEN", "BLOCK" if edge["severe"] else "REVIEW", edge,
            "Check raw coordinates, atom identity and bond order, and local ligand/pocket proximity before proposing constrained geometry repair. The screening bound is not an optimization restraint.")
    for metal in graph["metal_coordination"]:
        if metal.get("registered_nonbonded_mm_model") is True:
            continue
        add("METAL_MODEL_CONTRACT_REQUIRED", "BLOCK", metal,
            "Preserve metal and candidate ligands; resolve oxidation/protonation, MD coordination model, all ordered SKF pairs and FMO fragment representation. Distance alone does not authorize a coordination bond.")
    if graph["elements_without_all_ordered_SKF_pairs"]:
        add("MISSING_QM_ELEMENT_PAIR_PARAMETERS", "BLOCK", graph["elements_without_all_ordered_SKF_pairs"],
            "Provide and validate every ordered element-pair parameter, or explicitly mark the chemistry unsupported. Never drop atoms to pass QM.")
    for contact in graph.get("water_hydrogen_clashes", []):
        add("WATER_HYDROGEN_STERIC_CLASH", "REVIEW", contact,
            "Repair only flagged water hydrogen orientations against validated source geometry or an explicit fixed-heavy MM model; retain all heavy atoms and audit the full target.")
    if retention is not None:
        if retention.get("status") == "UNRESOLVED_RAW_INTAKE":
            add("RAW_COMPONENT_INVENTORY_UNRESOLVED", "BLOCK", retention,
                "Resolve alternate models/locations and register source-to-prepared mapping before accepting component retention.")
        for component in retention.get("components", []):
            if component["n_matched"] < component["n_heavy_atoms"]:
                add("RAW_SPECIAL_COMPONENT_NOT_FULLY_ACCOUNTED", "BLOCK" if component["contains_metal"] else "REVIEW", component,
                    "Map source ligand/cofactor/ion to the declared prepared model or record a justified exclusion/mutation. Coordinate mismatch is not by itself proof of deletion.")
    summary = {severity: sum(task["severity"] == severity for task in tasks) for severity in ("BLOCK", "REVIEW")}
    payload = {"schema": "label-free-chemical-quality-policy-v1", "experimental_labels_used": False,
        "input_graph_sha256": hashlib.sha256(json.dumps(graph, sort_keys=True).encode()).hexdigest(),
        "scope": "toolkit_diagnostic_not_automatic_coordinate_repair", "tasks": tasks, "summary": summary,
        "quality_status": "BLOCKED" if summary["BLOCK"] else "REVIEW_REQUIRED" if summary["REVIEW"] else "PASS_CURRENT_CHEMISTRY_CHECKS",
        "passing_does_not_establish_affinity_accuracy": True,
        "target_receptor_policy": "One registered receptor inventory per target; union ligand-dependent structural flags before scoring the whole cohort.",
        "validated_new_score_routes": [], "production_score_changed": False,
        "research_routes": [{"code": "COUPLED_PROTON_WATER_STATES", "enabled_in_production": False,
            "reason": "Completed BACE proton-location and four-target water experiments do not establish a transferable ranking improvement."}]}
    payload["plan_sha256"] = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
    return payload
