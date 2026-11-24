import pytest

from report_local_strain_transfer import ARMS, aggregate_crystals, score_metrics
from run_local_strain_transfer import optimization_recovery


def test_recovery_preserves_model_coordinates_and_strict_tolerance():
    deck = " $STATPT NSTEP=200 OPTTOL=0.000100 $END\n $PCM IDISP=0 $END\n $FMOXYZ\n H 1 0 0 0\n $END"
    for profile in ("original_200", "rfo_500", "qa_500"):
        recovered = optimization_recovery(deck, profile)
        assert "OPTTOL=0.000100" in recovered
        assert recovered.split(" $PCM", 1)[1] == deck.split(" $PCM", 1)[1]
    with pytest.raises(ValueError, match="unregistered"):
        optimization_recovery(deck, "choose_best_energy")


def test_openbind_uses_mean_not_best_crystal():
    compounds = [{"compound_index": "first", "fragalysis_codes": "first;second", "experimental_pKD": "7"}]
    cases = {"first": dict.fromkeys(ARMS, -10.), "second": dict.fromkeys(ARMS, -30.)}
    result = aggregate_crystals(compounds, cases)
    assert result["first"][ARMS[-1]] == -20.
    assert result["first"]["n_crystals"] == 2


def test_openbind_requires_complete_nonduplicated_mapping():
    compound = {"compound_index": "first", "fragalysis_codes": "first;missing", "experimental_pKD": "7"}
    with pytest.raises(ValueError, match="missing"):
        aggregate_crystals([compound], {"first": dict.fromkeys(ARMS, -10.)})
    compound["fragalysis_codes"] = "first;first"
    with pytest.raises(ValueError, match="duplicate"):
        aggregate_crystals([compound], {"first": dict.fromkeys(ARMS, -10.)})


def test_soft_compression_preserves_rank_and_does_not_reverse_score_sign():
    truth = {str(index): float(index) for index in range(5)}
    records = {name: dict.fromkeys(ARMS, -value) for name, value in truth.items()}
    result = score_metrics(truth, records)
    for arm in ARMS:
        assert result[arm]["raw"]["spearman"] == pytest.approx(1.)
        assert result[arm]["soft"]["spearman"] == pytest.approx(1.)
        assert result[arm]["soft"]["kendall"] == result[arm]["raw"]["kendall"]
    with pytest.raises(ValueError, match="full cohort"):
        score_metrics({"missing": 1.}, records)
