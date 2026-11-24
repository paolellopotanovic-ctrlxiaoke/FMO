import json

import pytest

from fmotpiefinal.core import graph_preparation


def test_untriggered_graph_is_byte_identical(tmp_path, monkeypatch):
    source = tmp_path / "source.pdb"
    source.write_text("registered receptor bytes\n")
    monkeypatch.setattr(graph_preparation, "audit", lambda path: {})
    monkeypatch.setattr(graph_preparation, "plan_actions", lambda graph: {
        "quality_status": "PASS_CURRENT_CHEMISTRY_CHECKS", "tasks": []})
    candidate, report = graph_preparation.prepare_graph_receptor(source, tmp_path / "graph")
    assert candidate.read_bytes() == source.read_bytes()
    assert report["status"] == "UNCHANGED_VERIFIED"
    assert not report["coordinates_modified"]


def test_carbonyl_only_repair_is_reported_as_modified_not_unchanged(tmp_path, monkeypatch):
    source = tmp_path / 'source.pdb'
    source.write_text('initial receptor\n')
    task = {'code': 'PROTEIN_BOND_LENGTH_SCREEN', 'evidence': {'source': 'A:GLY:1/C', 'target': 'A:GLY:1/O'}}
    monkeypatch.setattr(graph_preparation, 'audit', lambda path: {'repaired': path.read_text() == 'repaired receptor\n'})
    monkeypatch.setattr(graph_preparation, 'plan_actions', lambda graph: {
        'quality_status': 'PASS_CURRENT_CHEMISTRY_CHECKS' if graph['repaired'] else 'BLOCKED',
        'tasks': [] if graph['repaired'] else [task]})
    monkeypatch.setattr(graph_preparation, 'hydrogen_clashes', lambda path: [])
    def repair(source, tasks, output):
        output.write_text('repaired receptor\n')
        return {'selected_oxygens': ['A:GLY:1/O']}
    monkeypatch.setattr(graph_preparation, 'repair_carbonyl_oxygens', repair)
    candidate, report = graph_preparation.prepare_graph_receptor(source, tmp_path / 'graph')
    assert report['status'] == 'REPAIRED_AND_VERIFIED' and report['coordinates_modified']
    assert report['carbonyl_repair']['selected_oxygens'] == ['A:GLY:1/O']
    assert source.read_text() == 'initial receptor\n'


def test_unsupported_graph_tasks_stop_before_scoring(tmp_path, monkeypatch):
    source = tmp_path / "source.pdb"
    source.write_text("registered receptor bytes\n")
    monkeypatch.setattr(graph_preparation, "audit", lambda path: {})
    monkeypatch.setattr(graph_preparation, "plan_actions", lambda graph: {
        "quality_status": "BLOCKED", "tasks": [{"code": "METAL_MODEL_CONTRACT_REQUIRED"}]})
    with pytest.raises(ValueError, match="unresolved"):
        graph_preparation.prepare_graph_receptor(source, tmp_path / "graph")
    report = json.loads((tmp_path / "graph/graph_preparation.json").read_text())
    assert report["status"] == "REVIEW_REQUIRED"


def test_repair_must_pass_independent_second_graph(tmp_path, monkeypatch):
    source = tmp_path / "source.pdb"
    source.write_text("registered receptor bytes\n")
    task = {"code": "INVALID_COVALENT_LINK_GEOMETRY", "evidence": {"source": "A:VAL:1/C", "target": "A:GLN:2/N"}}
    monkeypatch.setattr(graph_preparation, "audit", lambda path: {})
    monkeypatch.setattr(graph_preparation, "plan_actions", lambda graph: {"quality_status": "BLOCKED", "tasks": [task]})
    monkeypatch.setattr(graph_preparation, "repair_amide_nitrogens", lambda raw, tasks, output: output.write_text("attempted repair\n"))
    with pytest.raises(ValueError, match="unresolved"):
        graph_preparation.prepare_graph_receptor(source, tmp_path / "graph")
    assert source.read_text() == "registered receptor bytes\n"
