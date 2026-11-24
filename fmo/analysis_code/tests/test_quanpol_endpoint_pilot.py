import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

from build_quanpol_endpoint_pilot import whole_components
from report_quanpol_endpoint_pilot import parse


def test_whole_components_preserves_disconnected_water_molecules():
    positions = np.array([[4.8, 0., 0.], [-4.3, 0., 0.],
                          [-4.8, 2., 0.], [4.3, 2., 0.]])
    expected = np.array([[4.8, 0., 0.], [5.7, 0., 0.],
                         [-4.8, 2., 0.], [-5.7, 2., 0.]])
    original = positions.copy()
    actual = whole_components(positions, [(0, 1), (2, 3)], np.eye(3) * 10.)
    np.testing.assert_allclose(actual, expected)
    np.testing.assert_array_equal(positions, original)


def test_whole_protein_component_is_not_folded_at_every_residue():
    positions = np.array([[1., 0., 0.], [4., 0., 0.], [7., 0., 0.], [0., 0., 0.]])
    actual = whole_components(positions, [(0, 1), (1, 2), (2, 3)], np.eye(3) * 10.)
    np.testing.assert_allclose(actual[:, 0], [1., 4., 7., 10.])


def test_isolated_ions_keep_coordinates():
    positions = np.array([[1., 2., 3.], [-4., -3., 0.]])
    np.testing.assert_array_equal(whole_components(positions, [], np.eye(3) * 10.), positions)


def test_parse_qm_center_energy_uses_kcal_record_not_hartree(tmp_path):
    log = tmp_path / "gamess.out"
    log.write_text(
        " QM CENTER        ENERGY = 0.0004235080\n"
        " QM CENTER                     ENERGY = 0.2657552496 KCAL/MOL\n"
        " EXECUTION OF GAMESS TERMINATED NORMALLY\n"
    )

    result = parse(log)

    assert result["qm_center_kcal_mol"] == 0.2657552496


def test_parse_rejects_ewald_abort_even_without_process_exit_code(tmp_path):
    log = tmp_path / "ewald.out"
    log.write_text(
        "ERROR: CURRENTLY NO EWALD CALCULATION IS ALLOWED FOR QMMM SYSTEMS.\n"
        "EXECUTION OF GAMESS TERMINATED -ABNORMALLY-\n"
    )
    assert parse(log)["status"] == "FAILED"
