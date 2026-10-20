import pytest

from fmotpiefinal.preparation_sensitivity import parse_propka_summary


def test_propka_parse_summary_not_determinants():
    text = ("ASP  93 A   4.52* 100 % more determinants\n"
            "SUMMARY OF THIS PREDICTION\nGroup pKa model-pKa\n"
            "   ASP  93 A   4.52   3.80\n   ASP 289 A 7.51 3.80\n"
            "   HIS  45 B 6.75 6.50\nFree energy of unfolding\n")
    assert parse_propka_summary(text) == {"A:ASP:93": 4.52, "A:ASP:289": 7.51, "B:HIS:45": 6.75}


@pytest.mark.parametrize("text", ["no summary", "SUMMARY OF THIS PREDICTION\n",
    "SUMMARY OF THIS PREDICTION\nASP 93 A 4.52 3.80\nASP 93 A 4.52 3.80\n",
    "SUMMARY OF THIS PREDICTION\nASP 93 A 4.52 3.80\nSUMMARY OF THIS PREDICTION"])
def test_propka_ambiguity_fails_closed(text):
    with pytest.raises(ValueError):
        parse_propka_summary(text)


def test_audit_has_no_access_to_affinity_properties():
    import inspect
    from fmotpiefinal.preparation_sensitivity import audit_preparation
    source = inspect.getsource(audit_preparation)
    assert "r_exp_dg" not in source
    assert "experimental_minus_dG" not in source
    assert "fep_rows" not in source
