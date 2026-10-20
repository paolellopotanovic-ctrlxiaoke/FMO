import pytest
import json
from pathlib import Path

from fmotpiefinal.core.scf_stall import MonomerSCFBudget


def cycle(budget, iteration, fragment=4, failed=True):
    budget.consume(f'RHF monomer {fragment} corr= NONE')
    if failed:
        budget.consume('SCF IS UNCONVERGED')
    return budget.consume(f'FMO iter {iteration}.1, 0 fragment(s) converged, maxE 1.0')


def test_three_consecutive_failed_cycles_require_retry_not_success():
    budget = MonomerSCFBudget()
    assert cycle(budget, 1) is None
    assert cycle(budget, 2) is None
    evidence = cycle(budget, 3)
    assert evidence['fragments'] == [4]
    assert evidence['requires_numerical_retry']
    assert not evidence['ranking_eligible']


def test_recovered_early_failure_is_not_aborted():
    budget = MonomerSCFBudget()
    assert cycle(budget, 1) is None
    assert cycle(budget, 2) is None
    assert cycle(budget, 3, failed=False) is None
    assert cycle(budget, 4) is None
    assert cycle(budget, 5) is None


def test_different_fragments_do_not_accumulate():
    budget = MonomerSCFBudget()
    assert cycle(budget, 1, fragment=1) is None
    assert cycle(budget, 2, fragment=2) is None
    assert cycle(budget, 3, fragment=1) is None


def test_dimer_failures_are_not_attributed_to_last_monomer():
    budget = MonomerSCFBudget()
    budget.consume('RHF monomer 4 corr= NONE')
    budget.consume('RHF dimer 4 2 corr= NONE')
    budget.consume('SCF IS UNCONVERGED')
    assert not budget.failed


def test_unreasonably_short_budget_rejected():
    with pytest.raises(ValueError):
        MonomerSCFBudget(1)


def test_persistent_outer_energy_cycle_is_bounded_without_accepting_result():
    budget = MonomerSCFBudget()
    for iteration in range(1, 25):
        value = [1.09e-4, 4.7e-7, 1.094e-4][(iteration - 1) % 3]
        evidence = budget.consume(f'FMO iter {iteration}.1, 30 fragment(s) converged, maxE {value}')
        if iteration < 24:
            assert evidence is None
    assert evidence['reason'] == 'OUTER_FMO_RESIDUAL_CYCLE_BUDGET'
    assert evidence['period'] == 3
    assert not evidence['ranking_eligible']


def test_decreasing_outer_residual_is_not_a_stalled_cycle():
    budget = MonomerSCFBudget()
    for iteration in range(1, 40):
        assert budget.consume(f'FMO iter {iteration}.1, 30 fragment(s) converged, maxE {0.9 ** iteration}') is None


def test_iteration_restart_does_not_join_residuals_from_different_attempts():
    budget = MonomerSCFBudget()
    for attempt in range(2):
        for iteration in range(1, 13):
            residual = 1e-4 if iteration % 2 else 2e-4
            assert budget.consume(f'FMO iter {iteration}.1, 30 fragment(s) converged, maxE {residual}') is None


@pytest.mark.skipif(not Path('/bin/tcsh').exists(), reason='GAMESS launcher requires tcsh')
def test_runner_stops_stalled_process_and_persists_retry_evidence(tmp_path):
    from fmotpiefinal.core.runner import run_gamess

    installation = tmp_path / 'gamess'
    installation.mkdir()
    (installation / 'gamess.00.x').write_bytes(b'fixture')
    script = ['#!/bin/tcsh']
    for iteration in range(1, 4):
        script.extend(['echo "RHF monomer 4 corr= NONE"', 'echo "SCF IS UNCONVERGED"',
                       f'echo "FMO iter {iteration}.1, 0 fragment(s) converged, maxE 1.0"'])
    script.append('sleep 20')
    (installation / 'rungms').write_text('\n'.join(script) + '\n')
    source = tmp_path / 'input.inp'
    source.write_text('fixture')
    output = tmp_path / 'run'
    returned = run_gamess(source, output, gamess_dir=installation, library_dir=tmp_path,
                          ligand_fragment=5, timeout_seconds=10)
    result = json.loads((output / 'result.json').read_text())
    assert returned == result
    assert result['status'] == 'FAILED_SCF_CONVERGENCE'
    assert 'numerical retry required' in result['error']
    assert result['numerical_budget']['fragments'] == [4]
    assert result['wall_seconds'] < 10
