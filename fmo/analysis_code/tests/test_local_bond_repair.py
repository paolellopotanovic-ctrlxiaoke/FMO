from fmotpiefinal.core.local_bond_repair import identify_amide_nitrogens


def task(first, second, code="PROTEIN_BOND_LENGTH_SCREEN"):
    return {"code": code, "evidence": {"source": first, "target": second}}


def test_graph_selects_only_same_residue_n_ca_or_cross_residue_c_n():
    tasks = [task("A:VAL:121/C", "A:GLN:122/N", "INVALID_COVALENT_LINK_GEOMETRY"),
             task("A:GLN:122/N", "A:GLN:122/CA"), task("A:GLN:122/CA", "A:GLN:122/N"),
             task("A:GLN:122/CD", "A:GLN:122/NE2"), task("A:GLN:122/C", "A:GLN:122/O"),
             task("A:GLN:122/N", "A:GLN:123/CA"), task("A:GLN:122/N", "A:GLN:122/CA", "NOT_A_BOND_FLAG")]
    assert identify_amide_nitrogens(tasks) == ["A:GLN:122/N"]


def test_no_trigger_does_not_select_a_receptor():
    assert identify_amide_nitrogens([]) == []
