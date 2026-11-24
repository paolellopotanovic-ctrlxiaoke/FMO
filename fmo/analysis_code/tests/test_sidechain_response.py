import copy

import pytest
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from fmotpiefinal.core.sidechain_response import decide_target, severe_overlaps
from fmotpiefinal.core.structure import Atom, Protein


def repaired_case():
    return {"baseline_new_amide_flags": [{"repaired": True}],
        "new_candidate_amide_flags": [], "severe_overlap_count_not_increased": True}


def test_known_artifact_with_verified_repair_is_accepted_without_scores():
    result = decide_target({"lig": repaired_case()}, ["lig"], geometry_verified=True, qm_verified=True)
    assert result["action"] == "ACCEPT_LOCAL_REPAIR_RESEARCH_CONTROL"
    assert not result["affinity_labels_used"]


@pytest.mark.parametrize("failure", ["geometry", "qm", "amide", "new_amide", "overlap"])
def test_incomplete_or_adverse_repair_is_not_accepted(failure):
    case = copy.deepcopy(repaired_case())
    if failure == "amide":
        case["baseline_new_amide_flags"][0]["repaired"] = False
    if failure == "new_amide":
        case["new_candidate_amide_flags"] = [{"new_deformation_flag": True}]
    if failure == "overlap":
        case["severe_overlap_count_not_increased"] = False
    result = decide_target({"lig": case}, ["lig"],
        geometry_verified=failure != "geometry", qm_verified=failure != "qm")
    assert result["action"] == "KEEP_BASELINE_REPAIR_REQUIRES_REVIEW"


def test_no_artifact_does_not_justify_unconditional_flexibility():
    case = repaired_case()
    case["baseline_new_amide_flags"] = []
    result = decide_target({"lig": case}, ["lig"], geometry_verified=True, qm_verified=True)
    assert result["action"] == "KEEP_BASELINE_MECHANISM_NOT_DETECTED"


def test_full_cohort_identity_required():
    with pytest.raises(ValueError):
        decide_target({"lig": repaired_case()}, ["lig", "missing"], geometry_verified=True, qm_verified=True)


def test_overlap_screen_includes_water_and_backbone_not_only_sidechains():
    molecule = Chem.AddHs(Chem.MolFromSmiles("C"))
    AllChem.EmbedMolecule(molecule, randomSeed=7)
    positions = np.asarray(molecule.GetConformer().GetPositions())
    positions -= positions[0]
    protein = Protein((Atom(0, "O", (2., 0., 0.), "A", "HOH", "1", "HETATM", "O"),
        Atom(1, "C", (-2., 0., 0.), "A", "ALA", "2", "ATOM", "CA")))
    result = severe_overlaps(protein, molecule, positions)
    assert {row["residue"][2] for row in result} == {"HOH", "ALA"}


def test_ordinary_polar_contact_is_not_a_gross_overlap():
    molecule = Chem.AddHs(Chem.MolFromSmiles("O"))
    AllChem.EmbedMolecule(molecule, randomSeed=7)
    positions = np.asarray(molecule.GetConformer().GetPositions())
    positions -= positions[0]
    protein = Protein((Atom(0, "N", (2.5, 0., 0.), "A", "ASN", "1", "ATOM", "ND2"),))
    assert severe_overlaps(protein, molecule, positions) == []
