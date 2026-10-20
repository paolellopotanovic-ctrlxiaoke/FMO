import json
import xml.etree.ElementTree as ET

import pytest

from fmotpiefinal.core.special_chemistry import assemble_retained, parameter_paths, write_phospho_bonds


def atom_line(serial, name, residue, number, element, position):
    return (f"HETATM{serial:5d} {name:>4} {residue:>3} A{number:4d}    "
            f"{position[0]:8.3f}{position[1]:8.3f}{position[2]:8.3f}  1.00  0.00          {element:>2}\n")


def test_supplemental_polymer_bonds_and_charge(tmp_path):
    destination = write_phospho_bonds(tmp_path / "bonds.xml")
    residues = {row.attrib["name"]: row for row in ET.parse(destination).getroot()}
    assert set(residues) == {"TPO", "SEP", "PTR"}
    bonds = {(bond.attrib["from"], bond.attrib["to"]) for bond in residues["TPO"]}
    assert ("-C", "N") in bonds and ("C", "+N") in bonds
    from openmm.app import ForceField
    forcefield = ForceField(*map(str, parameter_paths()))
    for name in residues:
        assert sum(atom.parameters["charge"] for atom in forcefield._templates[name].atoms) == pytest.approx(-2.)


def test_retained_ca_na_use_element_identity_and_validated_templates(tmp_path):
    raw = tmp_path / "raw.pdb"
    raw.write_text(atom_line(1, "CA", "CA", 261, "Ca", (5., 0., 0.)) +
                   atom_line(2, "NA1", "CA", 262, "Na", (9., 0., 0.)) + "END\n")
    prepared = tmp_path / "prepared.pdb"
    prepared.write_text(atom_line(1, "O", "HOH", 1, "O", (0., 0., 0.)) +
                        atom_line(2, "H1", "HOH", 1, "H", (.957, 0., 0.)) +
                        atom_line(3, "H2", "HOH", 1, "H", (-.24, .927, 0.)) + "END\n")
    result = assemble_retained(raw, prepared, tmp_path / "retained.pdb")
    assert result["all_residues_template_matched"]
    assert result["matched_special_templates"]["A:CA:261"]["charge"] == pytest.approx(2.)
    assert result["matched_special_templates"]["A:NA:262"]["charge"] == pytest.approx(1.)
    assert result["n_particles"] == 5
    assert result["coordination_bonds_added"] is False
    assert json.loads((tmp_path / "retained.model.json").read_text()) == result
    from fmotpiefinal.core.special_chemistry import registered_phosphoresidue_model
    from fmotpiefinal.core.receptor_chemical_graph import audit
    from fmotpiefinal.core.chemical_graph_policy import plan_actions
    output = tmp_path / 'retained.pdb'
    assert registered_phosphoresidue_model(output) == result
    assert plan_actions(audit(output))['quality_status'] == 'PASS_CURRENT_CHEMISTRY_CHECKS'
    model = output.with_suffix('.model.json')
    changed = json.loads(model.read_text())
    changed['matched_special_templates']['A:CA:261']['charge'] = 1.
    model.write_text(json.dumps(changed))
    with pytest.raises(ValueError, match='template charge changed'):
        registered_phosphoresidue_model(output)


def test_unknown_metal_cannot_be_deleted_or_renamed(tmp_path):
    raw, prepared = tmp_path / "raw.pdb", tmp_path / "prepared.pdb"
    raw.write_text(atom_line(1, "ZN", "ZN", 1, "Zn", (0., 0., 0.)) + "END\n")
    prepared.write_text("END\n")
    with pytest.raises(ValueError, match="no validated retained-component"):
        assemble_retained(raw, prepared, tmp_path / "retained.pdb")
    assert not (tmp_path / "retained.pdb").exists()
