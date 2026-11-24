import pytest

from fmotpiefinal.core.scc_convergence import fragment_failure_graph, mixer_evidence, scc_diis_deck


DECK = " $CONTRL MAXIT=200 $END\n $DFTB SCC=.TRUE. DFTB3=.TRUE. $END\n $FMO NFRAG=2 FRGNAM(1)=AASP1,LIG INDAT(1)=0 $END\n"


def test_only_actual_scc_mixer_changes():
    candidate = scc_diis_deck(DECK)
    assert candidate.replace("ITYPMX=2 ", "") == DECK
    with pytest.raises(ValueError):
        scc_diis_deck(candidate)


def test_ordinary_scf_diis_not_evidence():
    assert not mixer_evidence("DIIS=.TRUE. BROYDEN'S CHARGE MIXING")["verified_scc_diis"]
    assert mixer_evidence("DIIS CHARGE UPDATE")["verified_scc_diis"]
    assert not mixer_evidence("DIIS CHARGE UPDATE UNKNOWN ITYPMX")["verified_scc_diis"]


def test_fragment_headers_are_context_not_causality():
    graph = fragment_failure_graph(DECK, " RHF monomer 1\n SCF IS UNCONVERGED\n RHF dimer 2 1\n SCF IS UNCONVERGED")
    assert graph["unconverged_events"][0]["last_explicit_header"]["names"] == ["AASP1"]
    assert graph["unconverged_events"][1]["last_explicit_header"]["names"] == ["LIG", "AASP1"]
    assert graph["header_context_not_unique_causal_attribution"]


def test_invalid_fragment_rejected():
    with pytest.raises(ValueError):
        fragment_failure_graph(DECK, " RHF monomer 3")
