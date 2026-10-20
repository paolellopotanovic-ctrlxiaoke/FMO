import pytest

from fmotpiefinal.core.terminal_caps import cap_identity


@pytest.mark.parametrize("identifier,expected", [("298A", " 298A"), ("298", " 298 "),
    ("-1B", "  -1B"), ("9999", "9999 "), ("0", "   0 ")])
def test_cap_residue_fields_preserve_insertion_codes(identifier, expected):
    from fmotpiefinal.core.terminal_caps import pdb_residue_fields
    assert pdb_residue_fields(identifier) == expected


@pytest.mark.parametrize("identifier", ["298AB", "A298", "10000", "-1000", "", "1.5"])
def test_cap_residue_fields_reject_ambiguous_or_overwidth_ids(identifier):
    from fmotpiefinal.core.terminal_caps import pdb_residue_fields
    with pytest.raises(ValueError):
        pdb_residue_fields(identifier)


def test_structure_reader_preserves_residue_insertion_identity(tmp_path):
    from fmotpiefinal.core.structure import read_pdb, residue_key

    source = tmp_path / 'insertion.pdb'
    source.write_text('ATOM      1  CA  ALA A 298A      1.000   2.000   3.000  1.00 20.00           C\n')
    atom = read_pdb(source).atoms[0]
    assert atom.residue_number == '298A'
    assert residue_key(atom) == 'A:ALA:298A'


def methylamide():
    return {"id": "A:NMA:298A", "atoms": [{"name": name, "element": element} for name, element in
        (("N", "N"), ("CA", "C"), ("H", "H"), ("1HA", "H"), ("2HA", "H"), ("3HA", "H"))]}


def test_methylamide_alias_preserves_full_chemical_inventory():
    target, aliases = cap_identity(methylamide())
    assert target == "NME" and aliases["CA"] == "CH3"
    assert len(aliases) == 6


def test_missing_cap_atom_is_not_an_alias():
    record = methylamide()
    record["atoms"].pop()
    with pytest.raises(ValueError, match="complete atom"):
        cap_identity(record)


def test_unknown_nonstandard_cannot_be_renamed_to_cap():
    record = methylamide()
    record["id"] = "A:TLA:502"
    with pytest.raises(ValueError, match="unsupported"):
        cap_identity(record)


def test_cap_dispatch_is_chemical_not_target_specific():
    from fmotpiefinal.core.special_dispatch import choose_route
    record = methylamide()
    record["kind"] = "nonstandard_component"
    assert choose_route({"n_models": 1, "residues": [record]}) == "SOURCE_TERMINAL_CAP_MM"


def unknown_methylamide():
    return {"id": "_:UNK:1177", "kind": "nonstandard_component", "atoms": [
        {"name": name, "element": element, "xyz": xyz} for name, element, xyz in [
            ("C1", "C", (-29.865, 33.940, -15.759)), ("N1", "N", (-28.811, 32.886, -15.926)),
            ("H1", "H", (-28.273, 32.840, -16.779)), ("H2", "H", (-30.336, 33.833, -14.782)),
            ("H3", "H", (-29.408, 34.927, -15.836)), ("H4", "H", (-30.618, 33.828, -16.539))]]}


def test_unknown_cap_requires_atom_identical_complete_geometry():
    from fmotpiefinal.core.special_dispatch import choose_route
    record = unknown_methylamide()
    target, aliases = cap_identity(record)
    assert target == "NME" and aliases["C1"] == "CH3" and aliases["H1"] == "H"
    assert len(set(aliases.values())) == 6
    assert choose_route({"n_models": 1, "residues": [record]}) == "SOURCE_TERMINAL_CAP_MM"


def test_same_unknown_atom_counts_with_wrong_geometry_are_rejected():
    record = unknown_methylamide()
    record["atoms"][-1]["xyz"] = (0., 0., 0.)
    with pytest.raises(ValueError, match="connectivity|hydrogen"):
        cap_identity(record)


def test_unknown_acetyl_cap_graph():
    record = {"id": "_:UNK:890", "atoms": [{"name": name, "element": element, "xyz": xyz}
        for name, element, xyz in [("C1", "C", (-11.428, 14.683, -48.309)),
            ("O1", "O", (-12.334, 15.445, -48.549)), ("C2", "C", (-11.713, 13.288, -47.837)),
            ("H1", "H", (-10.774, 12.759, -47.676)), ("H2", "H", (-12.273, 13.327, -46.903)),
            ("H3", "H", (-12.300, 12.762, -48.590))]]}
    target, aliases = cap_identity(record)
    assert target == "ACE" and aliases["C1"] == "C" and aliases["C2"] == "CH3"
