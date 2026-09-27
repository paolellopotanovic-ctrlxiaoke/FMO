from pathlib import Path
import importlib.util
import json
import subprocess
import sys

from sophosqm_baseline.preparation import _filter_protein_input
from sophosqm_baseline.cohort import resolve_structure_path


def _load_runner():
    path = Path(__file__).resolve().parents[1] / "scripts" / "run_fmo87.py"
    spec = importlib.util.spec_from_file_location("run_fmo87", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_protein_input_filter_can_remove_all_crystallographic_waters(tmp_path: Path):
    source = tmp_path / "source.pdb"
    source.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N\n"
        "HETATM    2  O   HOH A 100       3.000   0.000   0.000  1.00 20.00           O\n"
        "HETATM    3  H1  HOH A 100       4.000   0.000   0.000  1.00 20.00           H\n",
        encoding="utf-8",
    )
    retained = tmp_path / "retained.pdb"
    removed = tmp_path / "removed.pdb"

    retained_count = _filter_protein_input(source, retained, keep_water=True)
    removed_count = _filter_protein_input(source, removed, keep_water=False)

    assert retained_count == 1
    assert removed_count == 1
    assert "HOH" in retained.read_text(encoding="utf-8")
    assert "HOH" not in removed.read_text(encoding="utf-8")
    assert "ATOM" in removed.read_text(encoding="utf-8")


def test_frozen_structure_paths_remap_to_local_structure_root(tmp_path: Path):
    source = tmp_path / "vendor" / "absent-checkout" / "cdk2_ligands.sdf"
    local = tmp_path / "assets" / "jacs_set" / "cdk2_ligands.sdf"
    local.parent.mkdir(parents=True)
    local.write_text("local\n", encoding="utf-8")

    assert resolve_structure_path(source, local.parent) == local


def test_existing_gamess_success_is_reused_only_for_matching_water_policy(tmp_path: Path):
    runner = _load_runner()
    job_dir = tmp_path / "job"
    minimized_dir = tmp_path / "minimized"
    job_dir.mkdir()
    minimized_dir.mkdir()
    result_path = job_dir / "run" / "result.json"
    result_path.parent.mkdir()
    result_path.write_text(json.dumps({"status": "SUCCEEDED"}), encoding="utf-8")
    (job_dir / "structure_manifest.json").write_text(
        json.dumps({"include_explicit_waters": True}), encoding="utf-8"
    )
    (minimized_dir / "minimization.json").write_text(
        json.dumps({"explicit_water_residues": 2}), encoding="utf-8"
    )

    assert runner._existing_success(
        result_path, job_dir, minimized_dir, "retained"
    ) is not None
    assert runner._existing_success(
        result_path, job_dir, minimized_dir, "removed_from_fmo_pocket"
    ) is None
    assert runner._existing_success(
        result_path, job_dir, minimized_dir, "removed_before_preparation"
    ) is None


def test_runner_rejects_preparation_water_policy_mismatch(tmp_path: Path):
    runner = _load_runner()
    metadata = tmp_path / "P38_preparation.json"
    metadata.write_text(json.dumps({"keep_water": True}), encoding="utf-8")

    try:
        runner._validate_preparation_policy(
            tmp_path, ["P38"], "removed_before_preparation"
        )
    except ValueError as error:
        assert "does not match explicit-water policy" in str(error)
    else:
        raise AssertionError("expected preparation policy mismatch")


def test_collect_fmo_shards_rebuilds_target_directories(tmp_path: Path):
    shard = tmp_path / "shards" / "0"
    job = shard / "fmo" / "CDK2_a"
    minimized = shard / "minimized" / "CDK2_a"
    job.mkdir(parents=True)
    minimized.mkdir(parents=True)
    (job / "result.json").write_text("{}", encoding="utf-8")
    (minimized / "minimization.json").write_text("{}", encoding="utf-8")
    pd = __import__("pandas")
    pd.DataFrame([{
        "target": "CDK2", "ligand_name": "a", "status": "SUCCEEDED",
        "tpie_kcal_mol": -1.0,
    }]).to_csv(shard / "fmo_features.csv", index=False)

    output = tmp_path / "targets"
    script = Path(__file__).resolve().parents[1] / "scripts" / "collect_fmo_shards.py"
    subprocess.run([sys.executable, str(script), "--shard-root", str(shard.parent),
                    "--output-root", str(output)], check=True)

    assert (output / "CDK2" / "fmo_features.csv").is_file()
    assert (output / "CDK2" / "fmo" / "CDK2_a" / "result.json").is_file()
    assert (output / "CDK2" / "minimized" / "CDK2_a" / "minimization.json").is_file()
