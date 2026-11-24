from fmotpiefinal.core.sidechain_network import extend_one_contact_shell
from fmotpiefinal.core.structure import Atom, Protein


def test_single_shell_is_nonrecursive_and_respects_protection():
    atoms = tuple(Atom(index, "C", (float(index*3), 0., 0.), "A", "LYS", str(index+1), "ATOM", "CB")
                  for index in range(4))
    nodes = [{"id": f"A:LYS:{index+1}", "residue": ["A", str(index+1), "LYS"],
              "protection": [], "action": "KEEP_FIXED_NO_TRIGGER", "evidence": []} for index in range(4)]
    nodes[3]["protection"] = ["METAL_COORDINATION_CANDIDATE"]
    plan = {"selected_residues": [["A", "1", "LYS"]], "nodes": nodes,
        "edges": [{"kind": "sidechain_ligand_contact", "receptor": node["id"]} for node in nodes],
        "policy": {"contact_radius_angstrom": 4.0}, "plan_sha256": "seed"}
    result = extend_one_contact_shell(Protein(atoms), plan)
    assert result["selected_residues"] == [["A", "1", "LYS"], ["A", "2", "LYS"]]
    assert plan["selected_residues"] == [["A", "1", "LYS"]]
    assert result["network_extension"]["recursive_expansion"] is False


def test_empty_seed_does_not_expand_nearby_nodes():
    plan = {"selected_residues": [], "nodes": [], "edges": [],
            "policy": {"contact_radius_angstrom": 4.0}, "plan_sha256": "seed"}
    assert extend_one_contact_shell(Protein(()), plan)["selected_residues"] == []


def test_two_shells_expand_only_to_second_neighbor_and_not_protected_nodes():
    atoms = tuple(Atom(index, "C", (float(index*3), 0., 0.), "A", "LYS", str(index+1), "ATOM", "CB")
                  for index in range(4))
    nodes = [{"id": f"A:LYS:{index+1}", "residue": ["A", str(index+1), "LYS"],
              "protection": [], "action": "KEEP_FIXED_NO_TRIGGER", "evidence": []} for index in range(4)]
    nodes[3]["protection"] = ["COVALENT_LINK"]
    plan = {"selected_residues": [["A", "1", "LYS"]], "nodes": nodes,
        "edges": [{"kind": "sidechain_ligand_contact", "receptor": node["id"]} for node in nodes],
        "policy": {"contact_radius_angstrom": 4.0}, "plan_sha256": "seed"}
    result = extend_one_contact_shell(Protein(atoms), plan, n_hops=2)
    assert result["selected_residues"] == [["A", "1", "LYS"], ["A", "2", "LYS"], ["A", "3", "LYS"]]
    assert result["network_extension"]["n_hops"] == 2
    assert not result["network_extension"]["unbounded_expansion"]


def test_onehop_complete_evidence_retains_untriggered_cohorts():
    import json
    from pathlib import Path
    report = json.loads((Path(__file__).parents[1] / "docs/GRAPH_SIDECHAIN_ONEHOP_RESULTS.zh-CN.json").read_text())
    assert report["final_status"]["status"] == "COMPLETE"
    assert report["total_new_qm"] == 10
    assert report["total_unchanged_reused"] == 260
    assert len(report["groups"]["CDK8"]["selected_residues"]) == 5
    for target, group in report["groups"].items():
        assert not group["failures"]
        if target != "CDK8":
            assert group["graph_guided"] == group["baseline"]
