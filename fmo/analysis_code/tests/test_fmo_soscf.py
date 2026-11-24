import pytest

from fmotpiefinal.core.fmo_soscf import enable_fmo_soscf, verify_fmo_soscf_log
from fmotpiefinal.pipeline.target_stages import collect_cohort_scores


DECK = (" $SCF CONV=1.0E-8 DIIS=.FALSE. SOSCF=.TRUE. $END\n"
        " $DFTB SCC=.TRUE. DFTB3=.TRUE. ITYPMX=-1 $END\n"
        " $FMOPRP MAXIT=200 NGUESS=2 MCONV(2)=530 MCONV(4)=530 $END\n")


def test_only_converger_bits_change():
    candidate = enable_fmo_soscf(DECK)
    assert candidate.replace("MCONV(2)=576 MCONV(4)=576", "MCONV(2)=530 MCONV(4)=530") == DECK
    assert 576 == 512 + 64


@pytest.mark.parametrize("deck", [DECK.replace("ITYPMX=-1", "ITYPMX=2"),
                                 DECK.replace("MCONV(4)=530", "MCONV(4)=534"),
                                 DECK + DECK])
def test_conflicting_controls_rejected(deck):
    with pytest.raises(ValueError):
        enable_fmo_soscf(deck)


@pytest.mark.parametrize("text", ["", "SOSCF IN EFFECT", "USE GAMESS STANDARD SUBROUTINES FOR SCF",
                                 "USE GAMESS STANDARD SUBROUTINES FOR SCF\nSOSCF IN EFFECT\nDIIS CHARGE UPDATE"])
def test_log_must_prove_actual_solver(tmp_path, text):
    path = tmp_path / "gamess.log"
    path.write_text(text)
    with pytest.raises(ValueError):
        verify_fmo_soscf_log(path)


def test_actual_solver_log(tmp_path):
    path = tmp_path / "gamess.log"
    path.write_text("USE GAMESS STANDARD SUBROUTINES FOR SCF\nSOSCF IN EFFECT\n")
    verify_fmo_soscf_log(path)


@pytest.mark.parametrize("workers", [1, 3])
def test_failure_does_not_skip_rest_of_cohort(workers):
    attempted = []

    def scoring(name, ligand):
        attempted.append(name)
        if name in {"0", "1", "2"}:
            raise RuntimeError("SCF failure")
        return {"ligand": name, "tpie": -10.0}

    jobs = [(index, str(index), None) for index in range(9)]
    with pytest.raises(RuntimeError, match="all 9 requested cases attempted"):
        collect_cohort_scores(jobs, workers, 9, scoring)
    assert sorted(attempted) == [str(index) for index in range(9)]


@pytest.mark.parametrize("workers", [1, 3])
def test_successful_cohort_unchanged(workers):
    jobs = [(index, str(index), None) for index in range(5)]
    rows = collect_cohort_scores(jobs, workers, 5, lambda name, ligand: {"ligand": name, "tpie": -2.0})
    assert sorted(row["ligand"] for row in rows) == [str(index) for index in range(5)]
