from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from audit_endpoint_water_coverage import coverage_counts


def test_nearby_water_is_not_assumed_to_be_in_qm():
    result = coverage_counts([1, 2, 3], [1, 4])
    assert result["n_nearby_represented_water_oxygens"] == 1
    assert result["n_nearby_omitted_water_oxygens"] == 2
    assert result["nearby_omitted_saved_atom_indices"] == [2, 3]


def test_no_polar_contact_is_zero_not_a_missing_water():
    result = coverage_counts([], [1, 4])
    assert result["n_nearby_water_oxygens"] == 0
    assert result["n_nearby_omitted_water_oxygens"] == 0
