from fmotpiefinal.core.chemical_model_requirements import model_requirements


def component(identifier, kind, atoms):
    return {"id": identifier, "model": "1", "kind": kind, "elements": {atom["element"]: 1 for atom in atoms}, "atoms": atoms}


def test_unknown_zinc_gets_specific_model_requirements_not_a_guessed_bond(tmp_path):
    graph = {"n_models": 1, "source_sha256": "a" * 64, "residues": [
        component("A:ZN:1", "metal", [{"name": "ZN", "element": "ZN", "xyz": [0, 0, 0]}]),
        component("A:HIS:2", "standard_residue", [{"name": "NE2", "element": "N", "xyz": [2, 0, 0]}])]}
    report = model_requirements(graph, tmp_path)
    zinc = report["components"][0]
    assert zinc["individual_mm_route"] is None
    assert zinc["candidate_metal_contacts"][0]["distance_angstrom"] == 2
    assert not zinc["candidate_metal_contacts"][0]["bond_assigned"]
    assert "metal_oxidation_and_spin_if_QM" in zinc["required_validation"]
    assert "Zn-N.skf" in report["missing_ordered_skf_pairs_for_observed_source_elements"]
    assert "N-Zn.skf" in report["missing_ordered_skf_pairs_for_observed_source_elements"]


def test_unknown_cofactor_is_parsed_without_inventing_charge_or_template():
    graph = {"n_models": 1, "source_sha256": "b" * 64, "residues": [
        component("A:XYZ:1", "nonstandard_component", [{"name": "C1", "element": "C", "xyz": [0, 0, 0]}])]}
    report = model_requirements(graph)
    row = report["components"][0]
    assert row["source_explicit_hydrogens"] == 0
    assert row["individual_mm_route"] is None
    assert "pinned_bond_orders_stereochemistry_and_covalent_attachment_check" in row["required_validation"]
    assert not row["automatic_qm_approval"]
