import json

import pytest

from fmotpiefinal.core.chemical_preflight import audit_chemistry, write_component_plan
from fmotpiefinal.core.component_inventory import enforce_inventory, inventory


def source_file(path, rows):
    lines = []
    for serial, (name, residue, element, point) in enumerate(rows, 1):
        lines.append(f"HETATM{serial:5d} {name:>4} {residue:>3} A{serial:4d}    "
                     f"{point[0]:8.3f}{point[1]:8.3f}{point[2]:8.3f}  1.00  0.00          {element:>2}")
    path.write_text("\n".join(lines) + "\nEND\n")
    return path


def test_plan_is_saved_before_intake_rejects_metal(tmp_path):
    source = source_file(tmp_path / "raw.pdb", [("ZN", "ZN", "Zn", (0., 0., 0.))])
    with pytest.raises(ValueError, match="unresolved chemistry"):
        enforce_inventory(source, tmp_path / "audit/raw_component_inventory.json")
    report = json.loads((tmp_path / "audit/component_model_plan.json").read_text())
    assert report["status"] == "UNSUPPORTED_CHEMISTRY"
    assert not report["automatic_qm_approval"]
    assert report["components"][0]["elements"] == {"ZN": 1}


def test_supported_mm_route_never_approves_qm(tmp_path):
    source = source_file(tmp_path / "raw.pdb", [("CL", "CL", "Cl", (0., 0., 0.))])
    before = source.read_bytes()
    report = write_component_plan(inventory(source), tmp_path / "plan.json")
    assert report["available_mm_route"] == "MONATOMIC_HALIDE_MM"
    assert report["status"] == "VALIDATED_MM_PREPARATION_REQUIRED"
    assert not report["ranking_eligible"]
    assert source.read_bytes() == before


def test_plan_cannot_silently_change_existing_registration(tmp_path):
    source = source_file(tmp_path / "raw.pdb", [("CA", "ALA", "C", (0., 0., 0.))])
    graph = inventory(source)
    write_component_plan(graph, tmp_path / "plan.json")
    graph["source_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="changed"):
        write_component_plan(graph, tmp_path / "plan.json")


def test_packaged_graph_audit_preserves_unsupported_chemistry(tmp_path):
    source = source_file(tmp_path / "raw.pdb", [("ZN", "ZN", "Zn", (0., 0., 0.))])
    before = source.read_bytes()
    report = audit_chemistry(source, tmp_path / "audit")
    assert report["status"] == "BLOCKED_UNSUPPORTED_CHEMISTRY"
    assert not report["coordinates_modified"]
    assert (tmp_path / "audit/receptor_chemical_graph.json").is_file()
    assert source.read_bytes() == before


def test_missing_cap_atoms_are_not_treated_as_supported(tmp_path):
    source = source_file(tmp_path / "raw.pdb", [("C", "ACE", "C", (0., 0., 0.))])
    report = write_component_plan(inventory(source), tmp_path / "plan.json")
    assert report["status"] == "UNSUPPORTED_CHEMISTRY"
    assert "complete atom" in report["blocking_reason"]


def test_unknown_component_without_ccd_returns_model_requirement_not_key_error(tmp_path):
    source = source_file(tmp_path / "raw.pdb", [("C1", "TLA", "C", (0., 0., 0.))])
    report = audit_chemistry(source, tmp_path / "audit")
    assert report["status"] == "BLOCKED_UNSUPPORTED_CHEMISTRY"
    graph = json.loads((tmp_path / "audit/receptor_chemical_graph.json").read_text())
    assert graph["residues"][0]["ccd"]["status"] == "CCD_REQUIRED"
    actions = json.loads((tmp_path / "audit/chemical_action_plan.json").read_text())
    assert any(task["code"] == "NONSTANDARD_COMPONENT_REQUIRES_MODEL" for task in actions["tasks"])


def test_ccd_roles_are_assigned_to_current_atom_not_previous_residue(tmp_path, monkeypatch):
    from fmotpiefinal.core import receptor_chemical_graph
    source = source_file(tmp_path / "raw.pdb", [("C1", "TLA", "C", (0., 0., 0.))])
    cache = tmp_path / "ccd"
    cache.mkdir()
    (cache / "TLA.cif").write_text("pinned dictionary placeholder for mock\n")
    monkeypatch.setattr(receptor_chemical_graph, "ccd_component", lambda *args: {
        "status": "CCD_MODEL_ANNOTATED", "atoms": {"C1": {"feature_families": ["Acceptor"]}}})
    graph = receptor_chemical_graph.audit(source, ccd_cache=cache)
    role = graph["atoms"]["A:TLA:1/C1"]["chemical_role"]
    assert role["acceptor"] and not role["donor"]
    assert role["role_uncertain"]


def test_failed_model_validation_does_not_publish_partial_receptor(tmp_path, monkeypatch):
    from fmotpiefinal.core import counterion_preparation
    from fmotpiefinal.core.special_dispatch import prepare_components
    source = source_file(tmp_path / "raw.pdb", [("CL", "CL", "Cl", (0., 0., 0.))])
    prepared = tmp_path / "prepared.pdb"
    prepared.write_text("unchanged prepared receptor\n")

    def failed_validation(raw, original, output):
        output.write_text("invalid candidate\n")
        raise ValueError("forcefield validation failed")

    monkeypatch.setattr(counterion_preparation, "assemble_counterions", failed_validation)
    output = tmp_path / "retained.pdb"
    with pytest.raises(ValueError, match="forcefield validation"):
        prepare_components(source, prepared, output)
    assert not output.exists()
    assert not output.with_suffix(".dispatch.json").exists()
    assert prepared.read_text() == "unchanged prepared receptor\n"
