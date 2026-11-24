import json

import pytest

from fmotpiefinal.core.structure_input import inspect_structure, normalize_structure


def atom(serial, name="CA", alt=" ", occupancy=1.0, residue="ALA", record="ATOM  ", element="C"):
    return f"{record}{serial:5d} {name:>4}{alt}{residue:>3} A   1       1.000   2.000   3.000{occupancy:6.2f}  0.00          {element:>2}\n"


def test_ordinary_pdb_byte_identical_and_sidecar_preserved(tmp_path):
    source = tmp_path / "source.pdb"
    source.write_text(atom(1) + "TER\nEND\n")
    source.with_suffix(".model.json").write_text('{"pinned": true}')
    result, report = normalize_structure(source, tmp_path / "normalized")
    assert result.read_bytes() == source.read_bytes()
    assert result.with_suffix(".model.json").read_bytes() == source.with_suffix(".model.json").read_bytes()
    assert report["byte_identical"]


def test_altloc_uses_one_residue_consistent_state_not_atomwise(tmp_path):
    source = tmp_path / "source.pdb"
    source.write_text(atom(1, "N", element="N") + atom(2, "CA", "A", .9) + atom(3, "CA", "B", .1)
                      + atom(4, "CB", "A", .3) + atom(5, "CB", "B", .7) + "TER\nEND\n")
    result, report = normalize_structure(source, tmp_path / "normalized")
    assert report["altloc_choices"][0]["selected"] == "A"
    assert report["output_atoms"] == 3
    raw = json.loads((tmp_path / "normalized/raw_component_inventory.json").read_text())
    assert sum(len(row["atoms"]) for row in raw["residues"]) == 5
    assert "TER\n" in result.read_text()
    assert all(line[16] == " " for line in result.read_text().splitlines() if line.startswith("ATOM"))


def test_multi_model_requires_explicit_selection_without_duplicate_atoms(tmp_path):
    source = tmp_path / "source.pdb"
    source.write_text("MODEL        1\n" + atom(1) + "ENDMDL\nMODEL        2\n" + atom(1) + "ENDMDL\nEND\n")
    with pytest.raises(ValueError, match="multiple source models"):
        normalize_structure(source, tmp_path / "ambiguous")
    assert json.loads((tmp_path / "ambiguous/normalization.json").read_text())["status"] == "REVIEW_REQUIRED"
    _, report = normalize_structure(source, tmp_path / "chosen", "2")
    assert report["output_atoms"] == 1
    assert report["excluded_models"] == ["1"]


def test_unknown_component_retained_and_not_claimed_qm_supported(tmp_path):
    source = tmp_path / "source.pdb"
    source.write_text(atom(1, "ZN", residue="ZN", record="HETATM", element="ZN"))
    _, report = normalize_structure(source, tmp_path / "normalized")
    graph = json.loads((tmp_path / "normalized/component_inventory.json").read_text())
    assert graph["residues"][0]["kind"] == "metal"
    assert not report["ranking_eligible"]


def test_mmcif_preserves_all_atoms_and_explicit_chain_mapping(tmp_path):
    gemmi = pytest.importorskip("gemmi")
    pdb = tmp_path / "source.pdb"
    pdb.write_text(atom(1))
    structure = gemmi.read_structure(str(pdb))
    structure[0][0].name = "long_chain"
    source = tmp_path / "source.cif"
    structure.make_mmcif_document().write_file(str(source))
    _, report = normalize_structure(source, tmp_path / "normalized")
    assert report["output_atoms"] == 1
    assert report["chain_mapping"] == {"long_chain": "A"}


def test_incomplete_alternate_is_not_silently_deleted(tmp_path):
    source = tmp_path / "source.pdb"
    source.write_text(atom(1, "CA", "A", .9) + atom(2, "CA", "B", .1) + atom(3, "CB", "B", .1))
    with pytest.raises(ValueError, match="lacks atoms"):
        normalize_structure(source, tmp_path / "normalized")


def test_duplicate_source_tartrate_hydrogen_names_reach_chemical_router(tmp_path):
    source = tmp_path / "source.pdb"
    source.write_text(atom(1, "HO", residue="TLA", record="HETATM", element="H")
                      + atom(2, "HO", residue="TLA", record="HETATM", element="H"))
    result, report = normalize_structure(source, tmp_path / "normalized")
    assert report["byte_identical"]
    assert result.read_bytes() == source.read_bytes()


def test_duplicate_protein_atom_is_not_silently_rebuilt(tmp_path):
    source = tmp_path / "source.pdb"
    source.write_text(atom(1) + atom(2))
    with pytest.raises(ValueError, match="ambiguous atom"):
        normalize_structure(source, tmp_path / "normalized")


def test_unknown_metal_can_be_inspected_without_claiming_support(tmp_path):
    source = tmp_path / "source.pdb"
    source.write_text(atom(1, "ZN", residue="ZN", record="HETATM", element="ZN"))
    report = inspect_structure(source, tmp_path / "inspection")
    assert report["status"] == "PARSED"
    assert report["component_plan"]["status"] == "UNSUPPORTED_CHEMISTRY"
    assert not report["ranking_eligible"]
    assert (tmp_path / "inspection/structure_input/raw_component_inventory.json").exists()


def test_mmcif_chain_remapping_preserves_disulfide_connections(tmp_path):
    gemmi = pytest.importorskip("gemmi")
    structure = gemmi.Structure()
    structure.add_model(gemmi.Model("1"))
    for index, name in enumerate(("long_a", "long_b")):
        chain = gemmi.Chain(name)
        residue = gemmi.Residue()
        residue.name, residue.seqid = "CYS", gemmi.SeqId(1, " ")
        sulfur = gemmi.Atom()
        sulfur.name, sulfur.element = "SG", gemmi.Element("S")
        sulfur.pos = gemmi.Position(index * 2, 0, 0)
        residue.add_atom(sulfur)
        chain.add_residue(residue)
        structure[0].add_chain(chain)
    connection = gemmi.Connection()
    connection.name, connection.type = "disulf1", gemmi.ConnectionType.Disulf
    connection.partner1 = gemmi.AtomAddress("long_a", gemmi.SeqId(1, " "), "CYS", "SG")
    connection.partner2 = gemmi.AtomAddress("long_b", gemmi.SeqId(1, " "), "CYS", "SG")
    structure.connections.append(connection)
    source = tmp_path / "source.cif"
    structure.make_mmcif_document().write_file(str(source))
    result, report = normalize_structure(source, tmp_path / "normalized")
    restored = gemmi.read_structure(str(result))
    assert report["chain_mapping"] == {"long_a": "A", "long_b": "B"}
    assert len(restored.connections) == 1
    assert restored.connections[0].partner1.chain_name == "A"
    assert restored.connections[0].partner2.chain_name == "B"
