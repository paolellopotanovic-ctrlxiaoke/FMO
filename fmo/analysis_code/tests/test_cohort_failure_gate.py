import json
from threading import Lock

import pytest

from fmotpiefinal.core.cohort_failure_gate import first_failed_monomer, repeated_scf_failure
from fmotpiefinal.pipeline.target_stages import collect_cohort_scores


DECK = 'FRGNAM(1)=AHIS57, AGLN61, ligand INDAT(1)=0 1 0 2 0 3 0'


def test_failed_monomer_uses_deck_inventory(tmp_path):
    log = tmp_path / 'gamess.log'
    log.write_text(' RHF      monomer      2 corr= NONE\nSCF IS UNCONVERGED\n')
    assert first_failed_monomer(DECK, log) == 'AGLN61'
    log.write_text(' RHF      dimer      1 2 corr= NONE\nSCF IS UNCONVERGED\n')
    assert first_failed_monomer(DECK, log) is None


def test_gate_requires_two_actual_scf_failures_with_same_fragment(tmp_path):
    for name in ['first', 'second']:
        folder = tmp_path / 'qm' / name
        (folder / 'run').mkdir(parents=True)
        (folder / 'input.inp').write_text(DECK)
        (folder / 'run/result.json').write_text(json.dumps({'status': 'FAILED_SCF_CONVERGENCE'}))
        (folder / 'run/gamess.log').write_text('RHF monomer 2 corr= NONE\nSCF IS UNCONVERGED\n')
    assert repeated_scf_failure(tmp_path, ['first']) is None
    assert repeated_scf_failure(tmp_path, ['first', 'second'])['fragment'] == 'AGLN61'
    (tmp_path / 'qm/second/run/result.json').write_text(json.dumps({'status': 'SUCCEEDED'}))
    assert repeated_scf_failure(tmp_path, ['first', 'second']) is None


@pytest.mark.parametrize('workers', [1, 6])
def test_repeated_failure_defers_unattempted_cases_without_partial_success(workers):
    attempted = []
    lock = Lock()

    def score(name, ligand):
        with lock:
            attempted.append(name)
        raise RuntimeError('SCF failed')

    jobs = [(index, str(index), None) for index in range(27)]
    with pytest.raises(RuntimeError, match='deferred to whole-cohort rescue'):
        collect_cohort_scores(jobs, workers, 27, score, failure_gate=lambda names: len(names) >= 2)
    assert len(attempted) == 2


def test_successful_probe_releases_entire_cohort():
    jobs = [(index, str(index), None) for index in range(27)]
    rows = collect_cohort_scores(jobs, 6, 27, lambda name, ligand: {'ligand': name, 'tpie': -1},
                                 failure_gate=lambda names: False)
    assert len(rows) == 27
