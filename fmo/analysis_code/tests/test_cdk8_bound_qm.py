"""Frozen-receptor geometry controls preserve model and atom identity."""
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from run_cdk8_bound_qm import converged_coordinates, optimization_deck


DECK = """ $CONTRL RUNTYP=ENERGY $END
 $FMO NFRAG=2 NBODY=2 NACUT=0 $END
 $PCM IDISP=0 ICAV=1 $END
 $FMOXYZ
 C 6 0 0 0
 H 1 1 0 0
 O 8 2 0 0
 $END
"""


def test_only_ligand_active_and_no_global_projection():
    result = optimization_deck(DECK, 2, 100)
    assert "IACTAT(1)=2,-3" in result
    assert "MODGRD=96 MODESP=0" in result
    assert "NZVAR=0" in result
    assert "$PCM IDISP=0 ICAV=1" in result
    assert all(len(line) <= 80 for line in result.splitlines())


def test_reject_selective_energy_and_invalid_budget():
    with pytest.raises(ValueError):
        optimization_deck(DECK.replace("NACUT=0", "NACUT=0 MODMOL=1"), 2, 100)
    with pytest.raises(ValueError):
        optimization_deck(DECK, 2, 0)


def test_frozen_receptor_is_checked_independently():
    text = """***** EQUILIBRIUM GEOMETRY LOCATED *****
 COORDINATES OF ALL ATOMS ARE (ANGS)
 ATOM CHARGE X Y Z
 ------------------
 C 6.0 0 0 0
 H 1.0 1.1 0 0
 O 8.0 2.1 0 0
 EXECUTION OF GAMESS TERMINATED NORMALLY
"""
    assert len(converged_coordinates(text, DECK, 2)) == 3
    with pytest.raises(ValueError, match="frozen receptor moved"):
        converged_coordinates(text.replace("C 6.0 0", "C 6.0 0.1"), DECK, 2)
    with pytest.raises(ValueError, match="did not converge"):
        converged_coordinates(text.replace("EQUILIBRIUM", "NOT EQUILIBRIUM"), DECK, 2)


def test_matching_free_reference_has_positive_strain_cost():
    from evaluate_cdk8_bound_qm import geometry_scores
    result = geometry_scores(-100, -70, -20, -21, -8)
    assert result["qm_strain_kcal_mol"] > 0
    assert result["qm_free_three_state"] > result["qm_bound_three_state"]
    assert result["qm_bound_strain_hybrid"] > result["qm_bound_tpie"]


def test_pieda_recovery_changes_only_output_control():
    from recover_cdk8_bound_pieda import pieda_deck
    source = DECK + " $FMOPRP IPIEDA=0 $END\n"
    result = pieda_deck(source)
    assert result.replace("IPIEDA=1", "IPIEDA=0") == source
    with pytest.raises(ValueError):
        pieda_deck(source.replace("NACUT=0", "NACUT=0 MODMOL=1"))


def test_fragment_membership_is_explicit_not_residue_assumption():
    from analyze_cdk8_geometry_evidence import fragments
    source = DECK.replace("NACUT=0", """NACUT=0
 ICHARG(1)=0,0
 FRGNAM(1)=ALA100,LIG
 INDAT(1)=0 1 -2 0 3 0""")
    names, membership = fragments(source)
    assert names == ["ALA100", "LIG"]
    assert membership == [[1, 2], [3]]


def test_bond_distance_audit_does_not_hide_disconnect():
    from analyze_cdk8_geometry_evidence import bond_audit
    from rdkit import Chem
    import numpy as np
    molecule = Chem.MolFromSmiles("CO")
    assert bond_audit(molecule, np.array([[0,0,0], [3,0,0]]))["n_flagged"] == 1


def test_amide_planarity_has_explicit_atom_mapping():
    from analyze_cdk8_geometry_evidence import amide_torsions
    from rdkit import Chem
    import numpy as np
    molecule = Chem.MolFromSmiles("CC(=O)NC")
    positions = np.array([[-1,1,0], [0,0,0], [0,1,0], [1,0,0], [1,-1,0]])
    records = amide_torsions(molecule, positions)
    assert len(records) == 1
    assert records[0]["departure_from_planarity_deg"] == pytest.approx(0)
    assert records[0]["atoms_zero_based"] == [2, 1, 3, 4]


def test_first_gradient_reads_requested_atom_not_frozen_atom(tmp_path):
    from check_cdk8_qm_gradient import first_gradient
    path = tmp_path / "gamess.log"
    path.write_text("""GRADIENT (HARTREE/BOHR)
 ATOM ZNUC DE/DX DE/DY DE/DZ
 1 N 7.0 1.0 2.0 3.0
 2 O 8.0 0.01 -0.02 0.03
""")
    assert first_gradient(path, 2) == [0.01, -0.02, 0.03]
    with pytest.raises(ValueError):
        first_gradient(path, 3)


def test_crystal_frame_alignment_preserves_proper_rotation():
    from audit_cdk8_crystal_frame import align
    import numpy as np
    reference = np.array([[0,0,0], [1,0,0], [0,1,0], [0,0,1]], dtype=float)
    target = reference @ np.array([[0,1,0], [-1,0,0], [0,0,1]]) + 10
    rotation, first, second = align(reference, target)
    assert np.linalg.det(rotation) == pytest.approx(1)
    assert (reference-first) @ rotation + second == pytest.approx(target)


def test_followup_retains_complete_negative_input_control():
    import json
    root = Path(__file__).resolve().parents[1]
    evidence = json.loads((root / "docs/CDK8_FOLLOWUP_ANALYSIS.json").read_text())
    assert evidence["n"] == 10
    assert evidence["production_changed"] is False
    assert evidence["completed_controls"]["input_pose"]["status"] == "COMPLETE"
    assert len(evidence["completed_controls"]["input_pose"]["rows"]) == 10
    assert evidence["metrics"]["input_pose_tpie"]["spearman"] == pytest.approx(-0.0545454545)
    assert len(evidence["completed_controls"]["gradient"]["rows"]) == 4


def test_full_qm_geometry_evidence_retains_all_cases_and_negative_result():
    import json
    root = Path(__file__).resolve().parents[1]
    evidence = json.loads((root / "docs/evidence/cdk8/geometry_evidence.json").read_text())
    assert evidence["n"] == len(evidence["rows"]) == 10
    assert len(evidence["pairs"]) == 45
    assert evidence["production_changed"] is False
    assert evidence["metrics"]["qm_bound_tpie"]["spearman"] == pytest.approx(0.1636363636)
    assert max(row["max_frozen_receptor_displacement"] for row in evidence["geometry"].values()) < 2e-8
    assert sum(row["bond_audit"]["n_flagged"] for row in evidence["geometry"].values()) == 0
