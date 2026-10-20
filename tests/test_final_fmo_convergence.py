import pytest
import json
from pathlib import Path

from fmotpiefinal.core.fmo_convergence import final_fmo_convergence
from fmotpiefinal.core.runner import classify_gamess_log


DECK = (" $CONTRL RUNTYP=ENERGY SCFTYP=RHF $END\n"
        " $SCF CONV=1.0E-8 $END\n $DFTB SCC=.TRUE. $END\n"
        " $FMOPRP CONV=1.0E-7 NGUESS=2 $END\n"
        " $PCM IFMO=-1 $END\n $FMO NFRAG=2 NBODY=2 $END\n")


def monomer(fragment, failed=False, tolerance="2.00E-08"):
    return (f" RHF monomer {fragment} corr= NONE\n DENSITY MATRIX CONV= {tolerance}\n" +
            (" SCF IS UNCONVERGED, TOO MANY ITERATIONS\n FINAL RDFTB ENERGY IS 0.0 AFTER 200 ITERATIONS\n"
             if failed else " DENSITY CONVERGED\n FINAL RDFTB ENERGY IS -10.0 AFTER 4 ITERATIONS\n"))


def example():
    return (monomer(1, failed=True, tolerance="2.00E-05") + monomer(2) +
            " FMO iter 1.1, 0 fragment(s) converged, maxD 0.1 maxE 1.0\n" +
            monomer(1) + monomer(2) +
            " FMO iter 2.1, 2 fragment(s) converged, maxD 0.0000001 maxE 0.00000001\n"
            " ------------Monomer SCF for layer 1 converged in 2 iterations!------------\n"
            " RHF dimer 2 1 corr= NONE\n DENSITY MATRIX CONV= 2.00E-08\n"
            " DENSITY CONVERGED\n FINAL RDFTB ENERGY IS -21.0 AFTER 4 ITERATIONS\n"
            " Layer 1 converged in 2 iterations (max200, conv .10E-06).\n"
            " EXECUTION OF GAMESS TERMINATED NORMALLY\n")


def audit(text, deck=DECK):
    return final_fmo_convergence(deck, text.splitlines())


def test_recovered_early_monomer_requires_complete_strict_final_state():
    result = audit(example())
    assert result["n_strict_final_monomers"] == 2
    assert result["n_strict_quantum_dimers"] == 1
    assert result["early_monomer_failures"][0]["fragment"] == 1


def test_nonconverged_initial_residual_overflow_does_not_invalidate_strict_final_state():
    text = example().replace("maxE 1.0", "maxE************")
    result = audit(text)
    assert result['nonconverged_residual_overflow_iterations'] == [1]
    assert result['n_strict_final_monomers'] == 2
    assert result['n_strict_quantum_dimers'] == 1


@pytest.mark.parametrize('replacement', ['maxE************', 'maxE -0.001'])
def test_final_residual_must_be_numeric_and_within_absolute_tolerance(replacement):
    with pytest.raises(ValueError):
        audit(example().replace('maxE 0.00000001', replacement))


def test_overflow_does_not_mask_bad_final_monomer_or_dimer():
    text = example().replace('maxE 1.0', 'maxE************')
    with pytest.raises(ValueError):
        audit(text.replace(monomer(1), monomer(1, tolerance='2.0E-5')))
    with pytest.raises(ValueError):
        audit(text.replace('FINAL RDFTB ENERGY IS -21.0', 'SCF IS UNCONVERGED\n FINAL RDFTB ENERGY IS -21.0'))


def test_residual_overflow_cannot_certify_partial_convergence():
    with pytest.raises(ValueError):
        audit(example().replace('0 fragment(s) converged, maxD 0.1 maxE 1.0',
                               '1 fragment(s) converged, maxD 0.1 maxE************'))


def test_runner_needs_registered_deck_and_explicit_strict_final_proof():
    assert classify_gamess_log(0, example()) == "FAILED_SCF_CONVERGENCE"
    assert classify_gamess_log(0, example(), deck=DECK) == "SUCCEEDED"
    assert classify_gamess_log(1, example(), deck=DECK) == "FAILED_SCF_CONVERGENCE"


@pytest.mark.skipif(not Path('/bin/tcsh').exists(), reason='GAMESS launcher requires tcsh')
def test_normal_exit_without_explicit_scf_failure_still_requires_strict_final_proof(tmp_path):
    from fmotpiefinal.core.runner import run_gamess
    installation = tmp_path / 'gamess'
    installation.mkdir()
    (installation / 'gamess.00.x').write_bytes(b'fixture')
    text = example().replace(monomer(1, failed=True, tolerance='2.00E-05'), monomer(1))
    text = text.replace('2.00E-08', '9.17E-08')
    assert 'SCF IS UNCONVERGED' not in text
    assert classify_gamess_log(0, text, deck=DECK) == 'SUCCEEDED'
    (installation / 'rungms').write_text('#!/bin/tcsh\ncat << LOGEND\n' + text + 'LOGEND\n')
    source = tmp_path / 'input.inp'
    source.write_text(DECK)
    result = run_gamess(source, tmp_path / 'run', gamess_dir=installation,
                        library_dir=tmp_path, ligand_fragment=2, timeout_seconds=10)
    assert result['status'] == 'FAILED_SCF_CONVERGENCE'
    assert 'numerical retry required' in result['error']
    assert json.loads((tmp_path / 'run/result.json').read_text()) == result


def test_free_monomer_guess_needs_actual_log_marker_and_same_final_proof():
    deck = DECK.replace('NGUESS=2', 'NGUESS=3')
    text = 'Using free molecule guess for monomers.\n' + example()
    assert audit(text, deck)['free_monomer_initial_guess_verified']
    with pytest.raises(ValueError):
        audit(text.replace('2 fragment(s) converged', '1 fragment(s) converged'), deck)


@pytest.mark.parametrize("text", [
    example().replace(monomer(1), monomer(1, failed=True)),
    example().replace(monomer(1), monomer(1, tolerance="2.00E-05")),
    example().replace("2 fragment(s) converged", "1 fragment(s) converged"),
    example().replace("maxE 0.00000001", "maxE 0.000001"),
    example().replace(" RHF dimer 2 1 corr= NONE", " SCF IS UNCONVERGED\n RHF dimer 2 1 corr= NONE"),
    example().replace(" FINAL RDFTB ENERGY IS -21.0", " SCF IS UNCONVERGED\n FINAL RDFTB ENERGY IS -21.0"),
    example().replace(" EXECUTION OF GAMESS TERMINATED NORMALLY", ""),
    example().replace(" Layer 1 converged in 2 iterations (max200, conv .10E-06).", ""),
    example().replace(monomer(2), ""),
    example().replace(" FINAL RDFTB ENERGY IS -21.0", " FINAL RDFTB ENERGY IS 0.0"),
    example() + " TERMINATED -ABNORMALLY-\n",
    "SCF IS UNCONVERGED\n" + example(),
])
def test_incomplete_unconverged_or_ambiguous_states_rejected(text):
    with pytest.raises(ValueError):
        audit(text)
    assert classify_gamess_log(0, text, deck=DECK) != "SUCCEEDED"


@pytest.mark.parametrize("deck", [DECK.replace("NBODY=2", "NBODY=1"),
                                  DECK.replace("NGUESS=2", "NGUESS=3"),
                                  DECK.replace("CONV=1.0E-8", "CONV=1.0E-5"),
                                  DECK.replace("NFRAG=2", "NFRAG=2 NLAYER=2"),
                                  DECK.replace("IFMO=-1", "IFMO=3")])
def test_other_protocols_do_not_gain_automatic_exception(deck):
    with pytest.raises(ValueError):
        audit(example(), deck)
