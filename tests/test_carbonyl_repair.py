import pytest

from fmotpiefinal.core.carbonyl_repair import selected_carbonyl_oxygens, repair_carbonyl_oxygens


def task(source, target, code='PROTEIN_BOND_LENGTH_SCREEN'):
    return {'code': code, 'evidence': {'source': source, 'target': target}}


def test_carbonyl_selection_requires_same_residue_backbone_C_O_names():
    tasks = [task('A:SER:87/C', 'A:SER:87/O'), task('A:GLY:111/O', 'A:GLY:111/C'),
             task('A:GLN:20/CD', 'A:GLN:20/OE1'), task('A:ALA:2/C', 'A:ALA:3/O'),
             task('A:ALA:2/N', 'A:ALA:2/CA'), task('A:ALA:2/C', 'A:ALA:2/O', 'OTHER')]
    assert selected_carbonyl_oxygens(tasks) == ['A:GLY:111/O', 'A:SER:87/O']


def test_carbonyl_repair_requires_trigger_and_preserves_source(tmp_path):
    source, output = tmp_path / 'source.pdb', tmp_path / 'output.pdb'
    source.write_text('unchanged\n')
    with pytest.raises(ValueError, match='graph-selected'):
        repair_carbonyl_oxygens(source, [], output)
    assert source.read_text() == 'unchanged\n' and not output.exists()


def test_carbonyl_special_component_combination_still_requires_joint_validation(tmp_path):
    source, output = tmp_path / 'source.pdb', tmp_path / 'output.pdb'
    source.write_text('unchanged\n')
    source.with_suffix('.components.json').write_text('{}')
    with pytest.raises(ValueError, match='joint validation'):
        repair_carbonyl_oxygens(source, [task('A:SER:87/C', 'A:SER:87/O')], output)
    assert not output.exists()
