import pytest

from fmotpiefinal.core.counterion_preparation import assemble_counterions


def pdb_atom(serial, name, residue, number, element, coordinates):
    return (f"HETATM{serial:5d} {name:>4} {residue:>3} A{number:4d}    "
            f"{coordinates[0]:8.3f}{coordinates[1]:8.3f}{coordinates[2]:8.3f}  1.00  0.00          {element:>2}\n")


def test_chloride_restored_with_negative_charge_without_changing_water(tmp_path):
    raw, prepared = tmp_path / "raw.pdb", tmp_path / "prepared.pdb"
    raw.write_text(pdb_atom(1, "CL", "CL", 301, "Cl", (9., 0., 0.)) + "END\n")
    water = (pdb_atom(1, "O", "HOH", 1, "O", (0., 0., 0.)) +
             pdb_atom(2, "H1", "HOH", 1, "H", (.957, 0., 0.)) +
             pdb_atom(3, "H2", "HOH", 1, "H", (-.24, .927, 0.)))
    prepared.write_text(water + "END\n")
    output = tmp_path / "retained.pdb"
    result = assemble_counterions(raw, prepared, output)
    assert result["matched_ion_templates"]["A:CL:301"]["charge"] == pytest.approx(-1.)
    assert result["n_particles"] == 4
    assert output.read_text().startswith(water)
    with pytest.raises(FileExistsError):
        assemble_counterions(raw, prepared, output)


def test_blank_chain_ion_registration_matches_graph_identity(tmp_path):
    from fmotpiefinal.core.receptor_chemical_graph import audit
    from fmotpiefinal.core.chemical_graph_policy import plan_actions
    raw, prepared = tmp_path / "raw.pdb", tmp_path / "prepared.pdb"
    raw.write_text(pdb_atom(1, "CL", "CL", 2, "Cl", (9., 0., 0.)).replace(" A   2", "     2") + "END\n")
    prepared.write_text(pdb_atom(1, "O", "HOH", 1, "O", (0., 0., 0.)) +
                        pdb_atom(2, "H1", "HOH", 1, "H", (.957, 0., 0.)) +
                        pdb_atom(3, "H2", "HOH", 1, "H", (-.24, .927, 0.)) + "END\n")
    output = tmp_path / "retained.pdb"
    result = assemble_counterions(raw, prepared, output)
    assert result["matched_ion_templates"]["_:CL:2"]["charge"] == -1.
    assert plan_actions(audit(output))["quality_status"] == "PASS_CURRENT_CHEMISTRY_CHECKS"


def test_counterion_route_rejects_a_metal(tmp_path):
    raw, prepared = tmp_path / "raw.pdb", tmp_path / "prepared.pdb"
    raw.write_text(pdb_atom(1, "ZN", "ZN", 1, "Zn", (0., 0., 0.)) + "END\n")
    prepared.write_text("END\n")
    with pytest.raises(ValueError, match="other chemistry must not be dropped"):
        assemble_counterions(raw, prepared, tmp_path / "retained.pdb")


def test_mixed_na_ca_cl_pass_full_mm_and_graph_without_claiming_qm(tmp_path):
    from fmotpiefinal.core.special_dispatch import prepare_components
    from fmotpiefinal.core.receptor_chemical_graph import audit
    from fmotpiefinal.core.chemical_graph_policy import plan_actions
    raw, prepared = tmp_path / "raw.pdb", tmp_path / "prepared.pdb"
    raw.write_text(pdb_atom(1, "NA1", "CA", 301, "Na", (9., 0., 0.)) +
                   pdb_atom(2, "CA", "CA", 302, "Ca", (13., 0., 0.)) +
                   pdb_atom(3, "CL", "CL", 303, "Cl", (17., 0., 0.)) + "END\n")
    prepared.write_text(pdb_atom(1, "O", "HOH", 1, "O", (0., 0., 0.)) +
                        pdb_atom(2, "H1", "HOH", 1, "H", (.957, 0., 0.)) +
                        pdb_atom(3, "H2", "HOH", 1, "H", (-.24, .927, 0.)) + "END\n")
    output = tmp_path / "retained.pdb"
    result = prepare_components(raw, prepared, output)
    assert result["dispatch"]["route"] == "MONATOMIC_ION_MM"
    assert result["schema"] == "retained-monatomic-ion-model-v1"
    assert result["matched_ion_templates"]["A:NA:301"]["charge"] == pytest.approx(1.)
    assert result["matched_ion_templates"]["A:CA:302"]["charge"] == pytest.approx(2.)
    assert result["matched_ion_templates"]["A:CL:303"]["charge"] == pytest.approx(-1.)
    graph = audit(output)
    assert len(graph["metal_coordination"]) == 2
    assert all(row["registered_nonbonded_mm_model"] for row in graph["metal_coordination"])
    assert plan_actions(graph)["quality_status"] == "PASS_CURRENT_CHEMISTRY_CHECKS"
    assert result["dispatch"]["ranking_eligible"] is False


def test_explicit_metal_bond_cannot_use_free_ion_route(tmp_path):
    from fmotpiefinal.core.special_dispatch import prepare_components
    raw, prepared = tmp_path / "raw.pdb", tmp_path / "prepared.pdb"
    raw.write_text(pdb_atom(1, "CA", "CA", 301, "Ca", (9., 0., 0.)) +
                   pdb_atom(2, "C", "ALA", 1, "C", (7., 0., 0.)) + "CONECT    1    2\nEND\n")
    prepared.write_text("END\n")
    output = tmp_path / "retained.pdb"
    with pytest.raises(ValueError, match="explicit ion connectivity"):
        prepare_components(raw, prepared, output)
    assert not output.exists()


def test_na_oxygen_source_connectivity_is_recorded_as_nonbonded_not_deleted(tmp_path):
    raw, prepared = tmp_path / "raw.pdb", tmp_path / "prepared.pdb"
    water = (pdb_atom(2, "O", "HOH", 1, "O", (0., 0., 0.)) +
             pdb_atom(3, "H1", "HOH", 1, "H", (.957, 0., 0.)) +
             pdb_atom(4, "H2", "HOH", 1, "H", (-.24, .927, 0.)))
    raw.write_text(pdb_atom(1, "NA", "NA", 301, "Na", (-2.4, 0., 0.)) + water + "CONECT    1    2\nEND\n")
    prepared.write_text(water + "END\n")
    result = assemble_counterions(raw, prepared, tmp_path / "retained.pdb")
    assert len(result["source_coordination_contacts"]) == 1
    assert result["source_coordination_contacts"][0]["distance_angstrom"] == pytest.approx(2.4)
    assert result["coordination_bonds_added"] is False
