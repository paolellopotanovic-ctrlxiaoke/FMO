import json
import sys


def test_single_crystal_failure_uses_sanitized_qm_key(tmp_path):
    from fmotpiefinal.pipeline.single_case_rescue import strict_failure

    qm = tmp_path / "qm"
    run = qm / "full_target_fmo/qm/crystal_one/run"
    run.mkdir(parents=True)
    (run / "result.json").write_text(json.dumps({"status": "FAILED_SCF_CONVERGENCE"}))
    (qm / "run_status.json").write_text(json.dumps({"status": "REVIEW_INCOMPLETE"}))
    assert strict_failure(tmp_path, "crystal one") == "FAILED_SCF_CONVERGENCE"


def test_registered_scf_rescues_do_not_request_unsupported_iteration_counts():
    from fmotpiefinal.pipeline.single_case_rescue import attempts

    assert attempts("FAILED_SCF_CONVERGENCE")[0] == (
        "fmo_soscf_strict_free_guess", ["--scf-solver", "fmo_soscf", "--scc-mixer", "gamess_soscf",
            "--soscf-onset", "1.0", "--nguess", "3", "--cnvdmp", "0"])

    for name, flags in attempts("FAILED_SCF_CONVERGENCE"):
        assert name != "diis_more_iterations"
        if "--maxit" in flags:
            assert 1 <= int(flags[flags.index("--maxit") + 1]) <= 200


def test_single_crystal_rescue_uses_ligand_id_not_shared_target(tmp_path, monkeypatch):
    from fmotpiefinal.pipeline import cli, single_case_rescue, target_workflow

    geometry = tmp_path / "geometry"
    geometry.mkdir()
    (geometry / "geometry_gate.json").write_text(json.dumps({
        "target": "OpenBind", "n_ligands": 1}))
    monkeypatch.setattr(target_workflow, "prerequisites",
                        lambda *args: ({}, None, None, ["crystal-123"], "A"))
    observed = []

    def strict_failure(case, ligand, qm):
        observed.append(ligand)
        return "FAILED_SCF_CONVERGENCE" if len(observed) > 1 else None

    def execute(*args, **kwargs):
        raise RuntimeError("registered recoverable QM failure")

    def rescue(ligand, *args, **kwargs):
        observed.append(ligand)
        return {"status": "RESCUED_ALTERNATE_PROTOCOL", "ranking_eligible": True}

    monkeypatch.setattr(single_case_rescue, "strict_failure", strict_failure)
    monkeypatch.setattr(single_case_rescue, "rescue", rescue)
    monkeypatch.setattr(target_workflow, "execute", execute)
    monkeypatch.setattr(sys, "argv", ["fmotpiefinal", "adaptive-run", "--geometry-dir",
                                      str(geometry), "--output", str(tmp_path / "qm"),
                                      "--gamess-dir", str(tmp_path), "--library-dir",
                                      str(tmp_path), "--dftb-parameter-dir", str(tmp_path)])
    cli.main()
    assert observed == ["crystal-123"] * 3
    assert json.loads((tmp_path / "adaptive_status.json").read_text())["ranking_eligible"]
