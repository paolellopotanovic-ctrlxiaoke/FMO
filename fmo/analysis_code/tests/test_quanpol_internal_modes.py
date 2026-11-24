import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

from diagnose_quanpol_internal_modes import minimize_hydrogens, rebind_parameters, replace_vacuum_coordinates, summarize


def test_parameter_rebinding_validates_species_and_records_content(tmp_path):
    (tmp_path / "C-H.skf").write_text("test parameter")
    updated, hashes = rebind_parameters(' $DFTBSK\n C H "/tmp/stale/sk/C-H.skf"\n $END\n', tmp_path)
    assert str(tmp_path / "C-H.skf") in updated
    assert len(hashes["C-H.skf"]) == 64
    with pytest.raises(ValueError, match="mismatch"):
        rebind_parameters(' C C "/tmp/stale/C-H.skf"\n', tmp_path)


def test_coordinate_replacement_preserves_qm_definition():
    text = " $CONTRL ICHARG=0 $END\n $DATA\nname\nC1\n C 6.0 0 0 0\n H 1.0 1 0 0\n $END\n"
    updated = replace_vacuum_coordinates(text, np.array([[0., 0., 0.], [1.1, .2, 0.]]))
    assert "ICHARG=0" in updated
    assert "H 1.0 1.1000000000 0.2000000000 0.0000000000" in updated
    with pytest.raises(ValueError, match="inventory"):
        replace_vacuum_coordinates(text, np.zeros((1, 3)))
    with pytest.raises(ValueError, match="isolated"):
        replace_vacuum_coordinates(text + " $FFDATA\n $END\n", np.zeros((2, 3)))


def test_hydrogen_minimization_ignores_diagnostic_duplicate_force():
    import openmm
    from openmm import unit

    system = openmm.System()
    system.addParticle(12)
    system.addParticle(1)
    bond = openmm.HarmonicBondForce()
    bond.addBond(0, 1, .1, 1000)
    system.addForce(bond)
    duplicate = openmm.CustomBondForce("100000*(r-0.5)^2")
    duplicate.addBond(0, 1, [])
    duplicate.setForceGroup(31)
    system.addForce(duplicate)
    integrator = openmm.VerletIntegrator(.001 * unit.picosecond)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName("Reference"))
    original = np.array([[0., 0., 0.], [1.5, 0., 0.]])
    relaxed, metadata = minimize_hydrogens(context, {0: "bond"}, original, [1])
    np.testing.assert_array_equal(relaxed[0], original[0])
    assert relaxed[1, 0] == pytest.approx(1.0, abs=1e-3)
    assert metadata["mm_energy_change_kcal_mol"] < 0
    assert original[1, 0] == 1.5


def test_summary_keeps_all_frames_and_rejects_missing_states():
    rows = []
    for repeat in (0, 1):
        for leg in ("bound", "solution"):
            for endpoint in (0, 10):
                value = endpoint * (2 if leg == "bound" else 1) + repeat
                rows.append({"repeat": repeat, "iteration": 400, "leg": leg, "endpoint": endpoint,
                             "original": {"internal_correction_kcal_mol": value},
                             "hydrogen_relaxed": {"internal_correction_kcal_mol": value * .5},
                             "minimization": {"heavy_atoms_displacement_angstrom": 0,
                                              "maximum_hydrogen_displacement_angstrom": .1}})
    report = summarize(rows)
    assert report["cycle_statistics"]["original"]["mean"] == 10
    assert report["cycle_statistics"]["hydrogen_relaxed"]["mean"] == 5
    with pytest.raises(ValueError, match="incomplete"):
        summarize(rows[:-1])
    with pytest.raises(ValueError, match="unique"):
        summarize(rows + rows[:1])
