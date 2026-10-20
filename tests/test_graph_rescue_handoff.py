import json

from fmotpiefinal.pipeline.cohort_rescue import preparation_stages_verified
from fmotpiefinal.pipeline.target_intake import digest


def test_graph_preparation_requires_report_hash_and_matching_receptor(tmp_path):
    prepared = tmp_path / 'receptor.pdb'
    prepared.write_text('registered receptor')
    report = tmp_path / 'graph_preparation/graph_preparation.json'
    report.parent.mkdir()
    report.write_text(json.dumps({'status': 'UNCHANGED_VERIFIED', 'output_sha256': digest(prepared)}))
    status = {'stage_completed': ['pdbfixer+geometry', 'chemical_graph+repair+reaudit', 'prepare', 'minimize'],
              'graph_preparation_sha256': digest(report)}
    assert preparation_stages_verified(tmp_path, status, prepared)
    prepared.write_text('changed receptor')
    assert not preparation_stages_verified(tmp_path, status, prepared)
    prepared.write_text('registered receptor')
    report.write_text(report.read_text() + '\n')
    assert not preparation_stages_verified(tmp_path, status, prepared)


def test_legacy_stages_cannot_hide_graph_stage_or_skip_minimization(tmp_path):
    prepared = tmp_path / 'receptor.pdb'
    prepared.write_text('receptor')
    status = {'stage_completed': ['pdbfixer+geometry', 'prepare', 'minimize']}
    assert preparation_stages_verified(tmp_path, status, prepared)
    status['graph_preparation_sha256'] = 'unverified'
    assert not preparation_stages_verified(tmp_path, status, prepared)
    status = {'stage_completed': ['pdbfixer+geometry', 'chemical_graph+repair+reaudit', 'prepare']}
    assert not preparation_stages_verified(tmp_path, status, prepared)
