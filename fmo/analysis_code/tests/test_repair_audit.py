import json

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from fmotpiefinal.core.conjugated_amine_geometry import conjugated_nh2_groups
from fmotpiefinal.core.repair_audit import (
    apply_experimental_nh2_repair,
    audit_local_strain_geometry,
)
from fmotpiefinal.pipeline.repair_audit import audit_target


def conjugated_fixture():
    molecule = Chem.AddHs(Chem.MolFromSmiles("NC(=[NH+](C))N(C)C"))
    assert AllChem.EmbedMolecule(molecule, randomSeed=20261025) == 0
    group = conjugated_nh2_groups(molecule)[0]
    group_index, center = group["nitrogen"], group["center"]
    other_n, hydrogens = group["other_nitrogens"], group["hydrogens"]
    points = molecule.GetConformer().GetPositions().copy()
    points[center] = [0., 0., 0.]
    points[group_index] = [1.35, 0., 0.]
    points[other_n[0]] = [-.675, 1.17, 0.]
    points[other_n[1]] = [-.675, -1.17, 0.]
    points[hydrogens[0]] = [1.85, np.sqrt(.75), 0.]
    points[hydrogens[1]] = [1.85, -np.sqrt(.75), 0.]
    distorted = points.copy()
    distorted[hydrogens[0]] = [1.85, 0., np.sqrt(.75)]
    return molecule, points, distorted, hydrogens


def test_strain_magnitude_never_selects_a_score_or_geometry_repair():
    molecule, original, _, _ = conjugated_fixture()
    high = audit_local_strain_geometry(molecule, original, original,
                                       local_strain_kcal_mol=1000.0)
    negative = audit_local_strain_geometry(molecule, original, original,
                                           local_strain_kcal_mol=-1000.0)
    assert high["decision"] == "NO_REPAIR_TRIGGER_TPIE_REMAINS_PRIMARY"
    assert negative["decision"] == "REVIEW_REQUIRED_NO_AUTOMATIC_SCORE_CHANGE"
    for report in (high, negative):
        assert report["repair_candidates"] == []
        assert report["primary_score"] == "pure_tPIE"
        assert report["local_strain_score_enabled"] is False
        assert report["labels_used_for_trigger"] is False


def test_k10_induced_nh2_flag_is_candidate_not_automatic_repair():
    molecule, original, distorted, _ = conjugated_fixture()
    report = audit_local_strain_geometry(molecule, original, distorted)
    assert report["decision"] == "REVIEW_REQUIRED_NO_AUTOMATIC_SCORE_CHANGE"
    assert len(report["repair_candidates"]) == 1
    assert report["repair_candidates"][0]["status"] == (
        "EXPERIMENTAL_ONLY_REQUIRES_FULL_K10_AND_QM_RERUN")
    assert report["coordinates_modified"] is False


def test_experimental_candidate_changes_only_flagged_hydrogens():
    molecule, original, distorted, hydrogens = conjugated_fixture()
    report = audit_local_strain_geometry(molecule, original, distorted)
    repaired, verification = apply_experimental_nh2_repair(molecule, distorted, report)
    heavy = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1]
    changed = {index for index in range(molecule.GetNumAtoms())
               if not np.array_equal(repaired[index], distorted[index])}
    assert changed <= set(hydrogens)
    assert np.array_equal(repaired[heavy], distorted[heavy])
    assert verification["heavy_atom_coordinates_unchanged"] is True
    assert verification["requires_full_k10_and_qm_rerun"] is True


def test_hard_geometry_gate_fails_closed():
    molecule, original, _, _ = conjugated_fixture()
    report = audit_local_strain_geometry(molecule, original, original,
                                        geometry_status="REPAIR_REQUIRED")
    assert report["decision"] == "REPAIR_REQUIRED_OR_REVIEW"
    assert report["findings"][0]["action"] == "STOP_AND_REPAIR_OR_REVIEW"
    with pytest.raises(ValueError, match="does not authorize"):
        apply_experimental_nh2_repair(molecule, original, report)


def test_nonsevere_geometry_warning_is_recorded_without_repair_trigger():
    molecule, original, _, _ = conjugated_fixture()
    report = audit_local_strain_geometry(molecule, original, original,
                                         geometry_status="REVIEW")
    assert report["decision"] == "NO_REPAIR_TRIGGER_TPIE_REMAINS_PRIMARY"
    assert report["findings"][0]["action"] == "RECORD_WARNING_NO_AUTOMATIC_REPAIR"


def _pdb_line(index, element, xyz):
    x, y, z = xyz
    atom_name = f"{element}{index:>3}"[-4:]
    return (f"HETATM{index:5d} {atom_name:4s} UNK X   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          "
            f"{element:>2s}\n")


def test_target_audit_writes_isolated_repair_candidates(tmp_path):
    molecule, original, distorted, _ = conjugated_fixture()
    molecule.SetProp("_Name", "guanidino")
    conformer = molecule.GetConformer()
    for index, point in enumerate(original):
        conformer.SetAtomPosition(index, point)
    sdf = tmp_path / "ligands.sdf"
    with Chem.SDWriter(str(sdf)) as writer:
        writer.write(molecule)
    minimized = tmp_path / "minimized" / "guanidino"
    minimized.mkdir(parents=True)
    pose = minimized / "minimized_complex.pdb"
    pose.write_text("".join(_pdb_line(index + 1, atom.GetSymbol().upper(), distorted[index])
                              for index, atom in enumerate(molecule.GetAtoms())))
    geometry = tmp_path / "post_k10_geometry.json"
    geometry.write_text(json.dumps({"cases": [{"ligand": "guanidino",
                                                "minimized_pose": {"quality_gate": "PASS"}}]}))
    output = tmp_path / "repair_audit.json"
    report = audit_target(sdf, tmp_path / "minimized", geometry, output,
                          emit_experimental_candidates=True)
    assert report["status"] == "COMPLETE"
    assert report["n_experimental_repair_candidates"] == 1
    assert report["ranking_modified"] is False
    candidate = tmp_path / "repair_audit_experimental_nh2_repair_candidates" / "guanidino" / "minimized_complex.pdb"
    assert candidate.is_file()
    assert pose.read_text() != candidate.read_text()
    assert json.loads(output.read_text())["labels_used_for_trigger"] is False
