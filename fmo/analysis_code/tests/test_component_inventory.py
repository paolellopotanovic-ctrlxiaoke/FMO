import json
from pathlib import Path

import pytest

from fmotpiefinal.core.component_inventory import inventory, enforce_inventory, sha


def write_atoms(path, components):
    lines = []
    for serial, (name, residue, element, record) in enumerate(components, 1):
        lines.append(f"{record:<6}{serial:5d} {name:>4} {residue:>3} A{serial:4d}    {serial:8.3f}{0.:8.3f}{0.:8.3f}  1.00  0.00          {element:>2}")
    path.write_text("\n".join(lines) + "\nEND\n")
    return path


def test_raw_inventory_uses_element_not_residue_name(tmp_path):
    source = write_atoms(tmp_path / "raw.pdb", [("CA", "CA", "Ca", "HETATM"), ("NA1", "CA", "Na", "HETATM")])
    result = inventory(source)
    assert [row["elements"] for row in result["residues"]] == [{"CA": 1}, {"NA": 1}]
    output = tmp_path / "audit.json"
    assert enforce_inventory(source, output)["status"] == "RAW_COMPONENT_CONTRACT_PASS"
    assert json.loads(output.read_text())["tasks"] == []


def test_modified_residue_detected_even_when_atom_record(tmp_path):
    source = write_atoms(tmp_path / "raw.pdb", [("P", "TPO", "P", "ATOM")])
    with pytest.raises(ValueError, match="unresolved chemistry"):
        enforce_inventory(source, tmp_path / "audit.json")
    assert inventory(source)["residues"][0]["kind"] == "nonstandard_component"


def test_complete_phosphoresidue_uses_registered_route_without_deletion(tmp_path):
    import xml.etree.ElementTree as ET
    from fmotpiefinal.core.special_chemistry import parameter_paths
    root = ET.parse(parameter_paths()[1]).getroot()
    template = next(row for row in root.find("Residues") if row.attrib["name"] == "TPO")
    elements = {row.attrib["name"]: row.attrib["element"].upper()
                for path in parameter_paths() for row in ET.parse(path).getroot().find("AtomTypes")}
    rows = [f"ATOM  {index:5d} {atom.attrib['name']:>4} TPO A 160    {float(index):8.3f}{0.:8.3f}{0.:8.3f}  1.00  0.00          {elements[atom.attrib['type']]:>2}"
            for index, atom in enumerate(template.findall("Atom"), 1)]
    source = tmp_path / "raw.pdb"
    source.write_text("\n".join(rows) + "\nEND\n")
    original = source.read_bytes()
    report = enforce_inventory(source, tmp_path / "audit.json")
    assert report["registered_mm_route"] == "PHOSPHORESIDUE_ALKALI_ALKALINE_EARTH_MM"
    assert source.read_bytes() == original and report["tasks"] == []


def test_clean_standard_input_does_not_change_coordinates(tmp_path):
    source = write_atoms(tmp_path / "raw.pdb", [("CA", "ALA", "C", "ATOM"), ("O", "HOH", "O", "HETATM")])
    before = source.read_bytes()
    report = enforce_inventory(source, tmp_path / "audit.json")
    assert report["status"] == "RAW_COMPONENT_CONTRACT_PASS"
    assert source.read_bytes() == before


def test_only_explicit_pinned_exclusion_allowed(tmp_path):
    source = write_atoms(tmp_path / "raw.pdb", [("CL", "CL", "Cl", "HETATM")])
    evidence = tmp_path / "evidence.json"
    evidence.write_text('{"structural_reason":"declared crystallization counterion"}')
    decisions = tmp_path / "decisions.json"
    data = {"source_sha256": sha(source), "experimental_labels_used": False, "components": [
        {"id": "A:CL:1", "action": "EXCLUDE", "reason": "declared nonbinding crystallization counterion", "evidence_path": str(evidence), "evidence_sha256": sha(evidence)}]}
    decisions.write_text(json.dumps(data))
    assert enforce_inventory(source, tmp_path / "audit.json", decisions)["status"] == "RAW_COMPONENT_CONTRACT_PASS"
    evidence.write_text("changed")
    with pytest.raises(ValueError, match="evidence changed"):
        enforce_inventory(source, tmp_path / "other.json", decisions)


def test_generic_path_cannot_approve_metal_deletion(tmp_path):
    source = write_atoms(tmp_path / "raw.pdb", [("ZN", "ZN", "Zn", "HETATM")])
    decisions = tmp_path / "decisions.json"
    decisions.write_text(json.dumps({"source_sha256": sha(source), "experimental_labels_used": False,
        "components": [{"id": "A:ZN:1", "action": "EXCLUDE", "reason": "not allowed to silently discard metal"}]}))
    with pytest.raises(ValueError, match="unresolved chemistry"):
        enforce_inventory(source, tmp_path / "audit.json", decisions)


def test_wrong_source_decisions_rejected(tmp_path):
    source = write_atoms(tmp_path / "raw.pdb", [("CA", "ALA", "C", "ATOM")])
    decisions = tmp_path / "decisions.json"
    decisions.write_text(json.dumps({"source_sha256": "0" * 64, "experimental_labels_used": False, "components": []}))
    with pytest.raises(ValueError, match="match raw source"):
        enforce_inventory(source, tmp_path / "audit.json", decisions)


def test_polymer_modified_residue_cannot_use_heteroatom_exclusion(tmp_path):
    source = write_atoms(tmp_path / "raw.pdb", [("P", "TPO", "P", "ATOM")])
    decisions = tmp_path / "decisions.json"
    decisions.write_text(json.dumps({"source_sha256": sha(source), "experimental_labels_used": False,
        "components": [{"id": "A:TPO:1", "action": "EXCLUDE", "reason": "cannot delete polymer component to pass"}]}))
    with pytest.raises(ValueError, match="unresolved chemistry"):
        enforce_inventory(source, tmp_path / "audit.json", decisions)


def test_intake_checks_raw_components_before_chain_selection(tmp_path):
    from fmotpiefinal.pipeline.target_intake import audit
    source = write_atoms(tmp_path / "raw.pdb", [("ZN", "ZN", "Zn", "HETATM")])
    ligands = tmp_path / "ligands.csv"
    ligands.write_text("ligand_id,smiles\none,CCO\n")
    poses = tmp_path / "poses.sdf"
    poses.write_text("not parsed before raw-component contract\n")
    with pytest.raises(ValueError, match="unresolved chemistry"):
        audit("new", ligands, tmp_path / "output", source, "auto", poses)
    assert (tmp_path / "output/raw_component_inventory.json").exists()
    assert not (tmp_path / "output/structural_cohort.csv").exists()
