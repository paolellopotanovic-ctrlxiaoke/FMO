import json
from types import SimpleNamespace

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from fmotpiefinal.core.structure import Atom, Ligand, Protein
from fmotpiefinal.pipeline import sidechain_response as response


@pytest.fixture
def setup_response(tmp_path, monkeypatch):
    molecule = Chem.AddHs(Chem.MolFromSmiles("C"))
    AllChem.EmbedMolecule(molecule, randomSeed=7)
    molecule.SetProp("_Name", "lig")
    sdf = tmp_path / "input.sdf"
    with Chem.SDWriter(str(sdf)) as writer:
        writer.write(molecule)
    protein_path = tmp_path / "protein.pdb"
    protein_path.write_text("source receptor\n")
    complex_path = tmp_path / "complex.pdb"
    complex_path.write_text("baseline complex\n")
    protein = Protein((Atom(0, "C", (10., 0., 0.), "A", "VAL", "1", "ATOM", "CB"),))
    points = np.asarray(molecule.GetConformer().GetPositions())
    ligand = Ligand("lig", tuple(Atom(index, atom.GetSymbol(), tuple(points[index]))
        for index, atom in enumerate(molecule.GetAtoms())))
    monkeypatch.setattr("fmotpiefinal.pipeline.target_stages.receptor_at_minimized_geometry", lambda *args: protein)
    rows = [{"ligand": "lig", "tpie": -10., "status": "SUCCEEDED"}]
    arguments = (protein_path, sdf, {"lig": complex_path}, protein, protein, [ligand],
                 SimpleNamespace(merge_disulfides=True, repair_peptide_numbering=True), tmp_path / "run", rows,
                 (protein, [], []))
    return arguments, molecule, ligand


def test_no_trigger_preserves_rows_and_never_minimizes(setup_response, monkeypatch):
    arguments, _, _ = setup_response
    monkeypatch.setattr("fmotpiefinal.core.minimization.minimize_complex",
                        lambda *args, **kwargs: pytest.fail("no-trigger path must not minimize"))
    assert response.apply_response(*arguments) is arguments[-2]
    decision = json.loads((arguments[-3] / "sidechain_response/decision.json").read_text())
    assert decision["action"] == "KEEP_BASELINE_MECHANISM_NOT_DETECTED"
    assert not decision["candidate_executed"]


def trigger(monkeypatch, *, repaired):
    def diagnose(*args):
        return {"baseline_new_amide_flags": [{"repaired": repaired}], "new_candidate_amide_flags": [],
                "severe_overlap_count_not_increased": True}
    monkeypatch.setattr(response, "diagnose_case", diagnose)
    monkeypatch.setattr(response, "audit_sidechains", lambda *args: {})
    monkeypatch.setattr(response, "add_displaced_clash_seeds", lambda *args: {})
    monkeypatch.setattr(response, "extend_one_contact_shell", lambda *args, **kwargs:
                        {"selected_residues": [["A", "1", "VAL"]]})


def test_failed_candidate_keeps_verified_baseline_and_records_failure(setup_response, monkeypatch):
    arguments, _, _ = setup_response
    trigger(monkeypatch, repaired=False)
    def fail(*args, **kwargs):
        raise ValueError("unsupported candidate chemistry")
    monkeypatch.setattr("fmotpiefinal.core.minimization.minimize_complex", fail)
    assert response.apply_response(*arguments) is arguments[-2]
    record = json.loads((arguments[-3] / "sidechain_response/decision.json").read_text())
    assert not record["accepted"]
    assert "unsupported candidate chemistry" in record["error"]
    assert response.apply_response(*arguments) is arguments[-2]


@pytest.mark.parametrize("qm_fails", [False, True])
def test_candidate_requires_complete_geometry_and_qm(setup_response, monkeypatch, qm_fails):
    arguments, _, ligand = setup_response
    trigger(monkeypatch, repaired=True)
    monkeypatch.setattr("fmotpiefinal.core.minimization.minimize_complex", lambda *args, **kwargs:
                        SimpleNamespace(complex_pdb=arguments[2]["lig"]))
    monkeypatch.setattr(response, "load_bound_ligand", lambda *args: ligand)
    monkeypatch.setattr("fmotpiefinal.core.fragmentation.make_facio_like_fragments", lambda *args, **kwargs:
                        (arguments[3], [SimpleNamespace()], {"hop_pairs_gamess_indices": []}))
    def score(*args):
        if qm_fails:
            raise ValueError("SCF failed")
        return {"ligand": "lig", "tpie": 100., "status": "SUCCEEDED"}
    monkeypatch.setattr("fmotpiefinal.pipeline.target_stages.score_one", score)
    actual = response.apply_response(*arguments)
    assert actual == (arguments[-2] if qm_fails else [{"ligand": "lig", "tpie": 100., "status": "SUCCEEDED"}])
    assert response.apply_response(*arguments) == actual


def test_changed_source_fails_closed_on_resume(setup_response):
    arguments, _, _ = setup_response
    response.apply_response(*arguments)
    arguments[0].write_text("changed receptor\n")
    with pytest.raises(ValueError, match="provenance changed"):
        response.apply_response(*arguments)


def test_duplicate_cohort_identity_is_rejected(setup_response):
    arguments, _, _ = setup_response
    arguments[-2].append(dict(arguments[-2][0]))
    with pytest.raises(ValueError, match="unique response inventory"):
        response.apply_response(*arguments)


def test_changed_cached_candidate_evidence_fails_closed(setup_response):
    arguments, _, _ = setup_response
    response.apply_response(*arguments)
    path = arguments[-3] / "sidechain_response/decision.json"
    record = json.loads(path.read_text())
    record["candidate_sha256"] = {str(arguments[2]["lig"]): "changed hash"}
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="cached response evidence changed"):
        response.apply_response(*arguments)


def test_candidate_cannot_pull_mm_only_component_into_qm_region(setup_response, monkeypatch):
    arguments, _, ligand = setup_response
    trigger(monkeypatch, repaired=True)
    monkeypatch.setattr("fmotpiefinal.core.minimization.minimize_complex", lambda *args, **kwargs:
                        SimpleNamespace(complex_pdb=arguments[2]["lig"]))
    monkeypatch.setattr(response, "load_bound_ligand", lambda *args: ligand)
    def invalid_region(*args):
        raise ValueError("retained MM component intersects actual ligand-defined QM region")
    monkeypatch.setattr("fmotpiefinal.pipeline.target_stages.validate_retained_component_region", invalid_region)
    monkeypatch.setattr("fmotpiefinal.pipeline.target_stages.score_one",
                        lambda *args: pytest.fail("QM must not start with an unsupported component"))
    assert response.apply_response(*arguments) is arguments[-2]
    decision = json.loads((arguments[-3] / "sidechain_response/decision.json").read_text())
    assert not decision["accepted"]
    assert "MM component intersects" in decision["error"]
