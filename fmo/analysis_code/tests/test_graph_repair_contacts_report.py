import pytest

from fmotpiefinal.core.parser import PairInteraction
from report_graph_repair_contacts import contact_map, ordering_errors


def test_contact_inventory_and_energy_closure():
    records = [PairInteraction(3, 1, -2.), PairInteraction(2, 3, 4.)]
    contacts = contact_map(records, {1: "ASP1", 2: "HOH2", 3: "LIG"}, 3)
    assert contacts == {"ASP1": -2., "HOH2": 4.}
    assert sum(contacts.values()) == 2.


def test_incomplete_or_duplicate_contact_inventory_is_rejected():
    with pytest.raises(ValueError, match="incomplete"):
        contact_map([PairInteraction(3, 1, -2.)], {1: "ASP1", 2: "HOH2"}, 3)
    with pytest.raises(ValueError, match="nonunique"):
        contact_map([PairInteraction(3, 1, -2.), PairInteraction(3, 2, 4.)], {1: "ASP1", 2: "ASP1"}, 3)


def test_ordering_uses_lower_energy_as_stronger_and_excludes_truth_ties():
    rows = ordering_errors(["strong", "weak"], {"strong": 9., "weak": 6.},
                           {"strong": -1., "weak": -2.}, {"strong": -3., "weak": -2.})
    assert rows[0]["category"] == "corrected"
    assert not ordering_errors(["first", "second"], {"first": 1., "second": 1.},
                               {"first": 1., "second": 2.}, {"first": 1., "second": 2.})


def test_still_wrong_does_not_mislabel_a_regression():
    rows = ordering_errors(["strong", "weak"], {"strong": 9., "weak": 6.},
                           {"strong": -3., "weak": -2.}, {"strong": -1., "weak": -2.})
    assert rows[0]["category"] == "regressed"
