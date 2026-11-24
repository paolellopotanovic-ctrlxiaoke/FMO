import json

import pytest

from fmotpiefinal.core.component_inventory import sha
from fmotpiefinal.core.water_hydrogen_clash import (
    hydrogen_clashes, repair_water_hydrogen_clashes, water_geometry_valid,
)
from fmotpiefinal.core.water_hydrogen_transfer import waters


def atom(serial, name, residue, number, position, element):
    return (f'HETATM{serial:5d} {name:4s} {residue:3s} A{number:4d}    '
            f'{position[0]:8.3f}{position[1]:8.3f}{position[2]:8.3f}'
            f'  1.00  0.00          {element:>2s}  \n')


def receptor(safe=False):
    return ''.join([
        atom(1, 'O', 'HOH', 1, (0, 0, 0), 'O'),
        atom(2, 'H1', 'HOH', 1, (-0.957 if safe else 0.957, 0, 0), 'H'),
        atom(3, 'H2', 'HOH', 1, (0.240 if safe else -0.240, 0.927, 0), 'H'),
        atom(4, 'HE21', 'GLN', 2, (1.8, 0, 0), 'H'),
        atom(5, 'O', 'HOH', 3, (10, 0, 0), 'O'),
        atom(6, 'H1', 'HOH', 3, (10.957, 0, 0), 'H'),
        atom(7, 'H2', 'HOH', 3, (9.760, 0.927, 0), 'H'),
    ])


def test_clash_detection_selects_water_and_not_intrawater_hydrogens(tmp_path):
    source = tmp_path / 'source.pdb'
    source.write_text(receptor())
    clashes = hydrogen_clashes(source)
    assert len(clashes) == 1
    assert clashes[0]['water_keys'] == [['A', '   1 ', 'HOH']]
    assert clashes[0]['distance_angstrom'] == pytest.approx(0.843)


def test_only_selected_water_hydrogens_change_and_metadata_is_rebound(tmp_path):
    source, original, output = [tmp_path / name for name in ('source.pdb', 'original.pdb', 'output.pdb')]
    source.write_text(receptor())
    original.write_text(receptor(safe=True))
    source.with_suffix('.model.json').write_text(json.dumps({'prepared_sha256': sha(source), 'ions': ['Ca']}))
    report = repair_water_hydrogen_clashes(source, original, output)
    assert report['n_changed_water_hydrogens'] == 2
    assert report['after'] == []
    assert not report['mm_fallback_used']
    assert output.read_text() == original.read_text()
    assert source.read_text() == receptor()
    model = json.loads(output.with_suffix('.model.json').read_text())
    assert model['prepared_sha256'] == sha(output)
    assert model['coordinate_repair_parent_sha256'] == sha(source)
    assert model['ions'] == ['Ca']


def test_untriggered_structure_is_byte_identical(tmp_path):
    source, output = tmp_path / 'source.pdb', tmp_path / 'output.pdb'
    source.write_text(receptor(safe=True))
    report = repair_water_hydrogen_clashes(source, None, output)
    assert report['n_changed_water_hydrogens'] == 0
    assert report['selected_water_keys'] == []
    assert source.read_bytes() == output.read_bytes()


def test_wrong_source_oxygen_cannot_transfer_hydrogens(tmp_path):
    source, original = tmp_path / 'source.pdb', tmp_path / 'original.pdb'
    source.write_text(receptor())
    original.write_text(receptor(safe=True).replace('   0.000   0.000   0.000', '   0.100   0.000   0.000', 1))
    with pytest.raises(ValueError, match='oxygen mapping'):
        repair_water_hydrogen_clashes(source, original, tmp_path / 'output.pdb')


def test_changed_registration_cannot_be_silently_reused(tmp_path):
    source, original = tmp_path / 'source.pdb', tmp_path / 'original.pdb'
    source.write_text(receptor())
    original.write_text(receptor(safe=True))
    source.with_suffix('.model.json').write_text(json.dumps({'prepared_sha256': 'wrong'}))
    with pytest.raises(ValueError, match='parent receptor changed'):
        repair_water_hydrogen_clashes(source, original, tmp_path / 'output.pdb')


def test_fresh_output_required(tmp_path):
    source = tmp_path / 'source.pdb'
    source.write_text(receptor())
    with pytest.raises(FileExistsError):
        repair_water_hydrogen_clashes(source, None, source)


def test_water_only_mm_fallback_fixes_clash_without_moving_oxygen(tmp_path):
    source, output = tmp_path / 'source.pdb', tmp_path / 'output.pdb'
    text = ''.join([
        atom(1, 'O', 'HOH', 1, (0, 0, 0), 'O'),
        atom(2, 'H1', 'HOH', 1, (0.957, 0, 0), 'H'),
        atom(3, 'H2', 'HOH', 1, (-0.240, 0.927, 0), 'H'),
        atom(4, 'O', 'HOH', 2, (2.6, 0, 0), 'O'),
        atom(5, 'H1', 'HOH', 2, (1.643, 0, 0), 'H'),
        atom(6, 'H2', 'HOH', 2, (2.840, -0.927, 0), 'H'),
    ])
    source.write_text(text)
    report = repair_water_hydrogen_clashes(source, None, output)
    assert report['mm_fallback_used']
    assert report['after'] == []
    assert report['all_other_atoms_byte_identical']
    assert all(water_geometry_valid(entry) for entry in waters(output).values())
    assert source.read_text().splitlines()[0] == output.read_text().splitlines()[0]
    assert source.read_text().splitlines()[3] == output.read_text().splitlines()[3]


def test_graph_preparation_executes_selected_water_repair_then_reaudits(tmp_path, monkeypatch):
    from fmotpiefinal.core import graph_preparation

    source, original = tmp_path / 'source.pdb', tmp_path / 'original.pdb'
    source.write_text(receptor())
    original.write_text(receptor(safe=True))
    audited = []

    def graph_audit(path):
        audited.append(str(path))
        return {'contacts': hydrogen_clashes(path)}

    monkeypatch.setattr(graph_preparation, 'audit', graph_audit)
    monkeypatch.setattr(graph_preparation, 'plan_actions', lambda graph: {
        'quality_status': 'REVIEW_REQUIRED' if graph['contacts'] else 'PASS_CURRENT_CHEMISTRY_CHECKS',
        'tasks': [{'code': 'WATER_HYDROGEN_STERIC_CLASH'}] if graph['contacts'] else []})
    candidate, report = graph_preparation.prepare_graph_receptor(source, tmp_path / 'graph', source_water=original)
    assert len(audited) == 2
    assert report['status'] == 'REPAIRED_AND_VERIFIED'
    assert report['water_hydrogen_repair']['n_changed_water_hydrogens'] == 2
    assert candidate.read_text() == receptor(safe=True)
    assert source.read_text() == receptor()
