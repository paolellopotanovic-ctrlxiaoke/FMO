import pytest

from fmotpiefinal.core.component_inventory import inventory
from fmotpiefinal.core.special_dispatch import choose_route, prepare_components


def graph(*components):
    return {"n_models": 1, "residues": [{"id": f"A:{name}:{number}", "kind": kind,
        "elements": {element: 1}, "atoms": [{"element": element}]} for number, (name, kind, element) in enumerate(components)]}


def test_dispatch_uses_chemistry_without_target_identity():
    assert choose_route(graph(("CA", "metal", "NA"), ("TPO", "nonstandard_component", "P"))) == "PHOSPHORESIDUE_ALKALI_ALKALINE_EARTH_MM"
    assert choose_route(graph(("CL", "ion", "CL"))) == "MONATOMIC_HALIDE_MM"
    assert choose_route(graph(("ALA", "standard_residue", "C"))) == "NO_SPECIAL_COMPONENTS"


def test_unknown_covalent_or_transition_metal_models_do_not_silently_fallback():
    for component in (("ZN", "metal", "ZN"), ("TLA", "nonstandard_component", "C")):
        with pytest.raises(ValueError, match="do not delete or mutate"):
            choose_route(graph(component))


def test_untriggered_prepared_coordinates_are_byte_identical(tmp_path):
    raw, prepared, output = tmp_path / "raw.pdb", tmp_path / "prepared.pdb", tmp_path / "output.pdb"
    raw.write_text("ATOM      1   CA ALA A   1       1.000   2.000   3.000  1.00  0.00           C\nEND\n")
    prepared.write_text("registered standard preparation stays byte-identical\n")
    assert inventory(raw)["residues"][0]["kind"] == "standard_residue"
    result = prepare_components(raw, prepared, output)
    assert output.read_bytes() == prepared.read_bytes()
    assert result["dispatch"]["route"] == "NO_SPECIAL_COMPONENTS"
    assert result["forcefield_system_created"] is False
