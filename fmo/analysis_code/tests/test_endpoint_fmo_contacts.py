from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from analyze_endpoint_fmo_contacts import descriptive_correlation, fragment_label, rotamer_fractions


def test_contact_residue_matches_fragment_label():
    assert fragment_label("ASP:A  93") == "AASP93"
    with pytest.raises(ValueError, match="ambiguous"):
        fragment_label("ASP:93")


def test_constant_contacts_do_not_create_nan_correlations():
    assert descriptive_correlation([1, 1, 1], [2, 3, 4]) is None
    assert descriptive_correlation([1, 2, 3], [3, 2, 1]) == pytest.approx(-1)


def test_rotamers_use_periodic_angles():
    fractions = rotamer_fractions([-179, 179, -60, 60])
    assert fractions == {"gauche_minus": .25, "gauche_plus": .25, "trans": .5}
