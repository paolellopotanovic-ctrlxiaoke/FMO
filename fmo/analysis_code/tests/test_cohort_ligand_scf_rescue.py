import json
import pytest


def test_failed_cohort_names_use_same_sanitized_key_as_qm_runner(tmp_path):
    from fmotpiefinal.pipeline.cohort_rescue import eligible_failure

    qm = tmp_path / "qm"
    run = qm / "full_target_fmo/qm/crystal_one/run"
    run.mkdir(parents=True)
    (run / "result.json").write_text(json.dumps({"status": "FAILED_SCF_CONVERGENCE"}))
    (qm / "run_status.json").write_text(json.dumps({"status": "REVIEW_INCOMPLETE", "names": ["crystal one"],
                                                  "failure": "CalledProcessError"}))
    assert eligible_failure(qm, ["crystal one"])["failed_cases"] == {"crystal one": "FAILED_SCF_CONVERGENCE"}


def test_numerical_rescue_keeps_original_pocket_and_geometry():
    from fmotpiefinal.pipeline.cohort_rescue import cohort_arm_controls, numerical_soscf_eligible

    assert numerical_soscf_eligible({"failed_cases": {"one": "FAILED_SCF_CONVERGENCE"}})
    assert not numerical_soscf_eligible({"failed_cases": {}})
    assert not numerical_soscf_eligible({"failed_cases": {"one": "FAILED_TIMEOUT"}})
    controls = cohort_arm_controls("whole_cohort_fmo_soscf", None)
    assert controls == ["--scf-solver", "fmo_soscf", "--scc-mixer", "gamess_soscf"]
    assert not any(flag in controls for flag in ("--pocket-heavy-only", "--exclude-pocket-residue", "--repair-clashes"))


@pytest.mark.parametrize('graph_enabled,repeated_gate', [(False, False), (True, False), (True, True)])
def test_numerical_route_without_pocket_membership_differences(tmp_path, monkeypatch, graph_enabled, repeated_gate):
    from types import SimpleNamespace
    from fmotpiefinal.pipeline import cohort_rescue as module

    qm = tmp_path / "qm"
    strict = qm / "full_target_fmo"
    strict.mkdir(parents=True)
    prepared = strict / "example_prepared.pdb"
    poses = strict / "ligands_no_properties.sdf"
    prepared.write_text("source receptor")
    poses.write_text("source ligands")
    protocol = strict / "frozen_protocol.json"
    protocol.write_text(json.dumps({"target_name": "Example", "ligand_names": ["a", "b"],
                                    "repair_peptide_numbering": True}))
    pocket = qm / "pocket_selection_audit.json"
    pocket.write_text(json.dumps({"ligand_names": ["a", "b"], "selection_production": "all_atom",
        "membership_differences": [], "source_sha256": {"receptor": module.digest(prepared),
        "ligands": module.digest(poses), "frozen_protocol": module.digest(protocol)}}))
    status = qm / "run_status.json"
    status.write_text(json.dumps({"repair_peptide_numbering": True, "pocket_heavy_only": False,
        "stage_completed": ["pdbfixer+geometry", "prepare", "minimize"],
        "pocket_selection_audit_sha256": module.digest(pocket)}))
    if graph_enabled:
        graph = qm / 'graph_preparation/graph_preparation.json'
        graph.parent.mkdir()
        graph.write_text(json.dumps({'status': 'UNCHANGED_VERIFIED', 'output_sha256': module.digest(prepared)}))
        record = json.loads(status.read_text())
        record['stage_completed'].insert(1, 'chemical_graph+repair+reaudit')
        record['graph_preparation_sha256'] = module.digest(graph)
        status.write_text(json.dumps(record))
    if repeated_gate:
        (strict / 'cohort_failure_gate.json').write_text(json.dumps({'fragment': 'AGLN61'}))
    monkeypatch.setattr(module, "eligible_failure", lambda *args: {
        "status_sha256": module.digest(status), "failed_cases": {"a": "FAILED_SCF_CONVERGENCE"}})
    monkeypatch.setattr(module, "ligand_scf_failures", lambda *args: False)
    commands = []

    def run(command, **kwargs):
        from pathlib import Path
        commands.append(command)
        output = Path(command[command.index("--output") + 1])
        if command[command.index("--stage") + 1] == "score":
            (output / "frozen_protocol.json").write_text(json.dumps({"scf_solver": "fmo_soscf",
                "scc_mixer": "gamess_soscf", "repair_peptide_numbering": True,
                "soscf_onset": 1.0,
                "nguess": 3 if '--nguess' in command else 2,
                "cnvdmp": 0 if '--cnvdmp' in command else 50,
                "prepared_protein_override_sha256": module.digest(prepared),
                "poses_protocol_sha256": module.digest(protocol)}))
            (output / "blind_predictions.json").write_text(json.dumps({"n_requested": 2, "n_succeeded": 2,
                "predictions": [{"ligand": name, "status": "SUCCEEDED", "tpie": -10.0} for name in ["a", "b"]]}))
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(module.subprocess, "run", run)
    report = module.rescue_cohort("Example", qm, ["a", "b"], tmp_path, 2, 900, 2, tmp_path, tmp_path)
    assert report["selected_arm"] == ('whole_cohort_fmo_soscf_strict_free_guess' if repeated_gate else 'whole_cohort_fmo_soscf_early')
    assert report["ranking_eligible"]
    assert len(commands) == 2
    assert all("--pocket-heavy-only" not in command and "--exclude-pocket-residue" not in command for command in commands)
    assert all(command[command.index("--poses-root") + 1] == str(strict) for command in commands)


def test_ligand_scf_arm_only_if_every_failure_is_ligand_monomer(tmp_path):
    from fmotpiefinal.pipeline.cohort_rescue import ligand_scf_failures

    qm = tmp_path / "qm"
    failures = {"crystal-one": "FAILED_SCF_CONVERGENCE",
                "crystal-two": "FAILED_SCF_CONVERGENCE"}
    for name in failures:
        run = qm / "full_target_fmo/qm" / name / "run"
        run.mkdir(parents=True)
        (run / "result.json").write_text(json.dumps({"ligand_fragment": 48}))
        (run / "gamess.log").write_text("RHF      monomer     48 corr= NONE\n"
                                         "SCF IS UNCONVERGED, TOO MANY ITERATIONS\n")
    assert ligand_scf_failures(qm, failures)
    assert not ligand_scf_failures(qm, {**failures, "other": "FAILED_TIMEOUT"})
    (qm / "full_target_fmo/qm/crystal-two/run/gamess.log").write_text(
        "RHF      monomer     47 corr= NONE\nSCF IS UNCONVERGED\n")
    assert not ligand_scf_failures(qm, failures)
