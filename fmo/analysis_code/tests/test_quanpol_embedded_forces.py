import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from build_quanpol_force_export import instrument
from validate_quanpol_embedded_forces import BOHR_ANGSTROM, HARTREE_KCAL, parse_embedded, replace_coordinates


def example():
    return """QMMM_FORCE_EXPORT_BEGIN 1 2 0.0
QMG 1 1.0 0.0 0.0
MMG 1 1 1.0 0.0 0.0
MMG 2 0 -1.0 0.0 0.0
QMMM_FORCE_EXPORT_END
POTENTIAL ENERGY = -12.345 KCAL/MOL
EXECUTION OF GAMESS TERMINATED NORMALLY
"""


def test_gradient_units_no_duplicate_count_and_final_energy():
    result = parse_embedded(example())
    assert result["energy_kcal_mol"] == -12.345
    assert result["gradient_kcal_mol_angstrom"].shape == (2, 3)
    assert result["gradient_kcal_mol_angstrom"][0, 0] == pytest.approx(HARTREE_KCAL / BOHR_ANGSTROM)
    np.testing.assert_allclose(result["gradient_kcal_mol_angstrom"].sum(0), 0)


@pytest.mark.parametrize("change", [lambda text: text.replace("MMG 2", "MMG 3"),
                                       lambda text: text.replace("MMG 1 1 1.0", "MMG 1 1 2.0"),
                                       lambda text: text + "SCF IS UNCONVERGED", lambda text: text + text,
                                       lambda text: text.replace("POTENTIAL ENERGY", "TOTAL ENERGY")])
def test_gradient_fail_closed(change):
    with pytest.raises(ValueError):
        parse_embedded(change(example()))


def test_qm_atom_perturbation_updates_both_coordinate_copies():
    text = " $DATA\ntitle\nC1\n C 6 1 2 3\n $END\n $FFDATA\nCOORDINATES\n C 6 1 2 3\n O 8 4 5 6\nSTOP\n $END\n"
    result = replace_coordinates(text, [[1.1, 2, 3], [4, 5, 6]], 1)
    assert result.count("1.100000000000") == 2
    with pytest.raises(ValueError):
        replace_coordinates(text, [[1, 2, 3]], 1)


def test_instrumentation_is_output_only_before_force_masking():
    source = "      SUBROUTINE FFOPT2(\nC\nC        - ZERO OFF SOME FORCES -\n"
    result = instrument(source)
    assert result.index("QMMM_FORCE_EXPORT_BEGIN") < result.index("ZERO OFF SOME FORCES")
    assert "CALL " not in result
    with pytest.raises(ValueError):
        instrument(result)
