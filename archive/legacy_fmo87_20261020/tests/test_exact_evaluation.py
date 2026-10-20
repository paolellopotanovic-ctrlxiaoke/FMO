import numpy as np
import pandas as pd
import pytest
import hashlib
import json

from sophosqm_baseline import exact_evaluation as ee


def _frames():
    rows = []
    for target in ("A", "B"):
        for ligand in range(4):
            rows.append({
                "target": target,
                "ligand_name": str(ligand),
                "canonical_isomeric_smiles": "C",
                "experimental_dg_kcal_mol": -7.0 + ligand + (0.5 if target == "B" else 0.0),
                "rdkit_crippen_clogp": 1.0 + ligand * ligand / 10.0,
            })
    cohort = pd.DataFrame(rows)
    fmo = cohort.copy()
    fmo["status"] = "SUCCEEDED"
    fmo["tpie_kcal_mol"] = -60.0 + np.arange(len(fmo)) * 2.0
    boltz = cohort.rename(columns={"ligand_name": "ligand_id"})
    boltz["status"] = "SUCCEEDED"
    boltz["affinity_pred_value"] = 5.5 - np.arange(len(boltz)) * 0.1
    return cohort, fmo, boltz


def test_exact_evaluation_strict_and_produces_model_layers(tmp_path, monkeypatch):
    monkeypatch.setattr(ee, "EXPECTED_COUNTS", {"A": 4, "B": 4})
    monkeypatch.setattr(ee, "PAPER_BOLTZ2", {"pearson_r": 0.7, "kendall_tau": 0.5})
    monkeypatch.setattr(ee, "PAPER_FMO", {"pearson_r": 0.6, "kendall_tau": 0.4})
    cohort, fmo, boltz = _frames()
    paths = []
    for name, frame in (("cohort.csv", cohort), ("fmo.csv", fmo), ("boltz.csv", boltz)):
        path = tmp_path / name
        frame.to_csv(path, index=False)
        paths.append(path)
    result = ee.evaluate_exact_fmo87(*paths, bootstrap_repetitions=5)
    assert result.audit["strict_87_complete"]
    assert result.audit["all_required_runs_succeeded"]
    assert set(result.compounds.columns) >= {
        "boltz_dg_kcal_mol", "fmo_target_fit_dg_kcal_mol", "fmo_loo_dg_kcal_mol"
    }
    assert set(result.models["model"]) == {
        "boltz2_local", "fmo_raw_tpie", "fmo_target_fit", "fmo_loo"
    }
    assert all(np.isfinite(result.compounds["fmo_loo_dg_kcal_mol"]))
    weighted = result.models.set_index("model")["sample_weighted"].apply(pd.Series)
    assert np.isfinite(weighted[[
        "pairwise_mae_kcal_mol", "mae_noncentered_kcal_mol", "mae_centered_kcal_mol",
        "percent_within_1_noncentered_kcal_mol", "percent_within_1_centered_kcal_mol",
        "percent_within_2_noncentered_kcal_mol", "percent_within_2_centered_kcal_mol",
    ]].to_numpy()).all()

    incomplete = fmo.iloc[1:].copy()
    incomplete_path = tmp_path / "incomplete_fmo.csv"
    incomplete.to_csv(incomplete_path, index=False)
    with pytest.raises(ValueError, match="incomplete FMO"):
        ee.evaluate_exact_fmo87(paths[0], incomplete_path, paths[2])


def test_exact_evaluation_rejects_label_mismatch(tmp_path, monkeypatch):
    monkeypatch.setattr(ee, "EXPECTED_COUNTS", {"A": 4, "B": 4})
    cohort, fmo, boltz = _frames()
    fmo.loc[fmo.index[0], "experimental_dg_kcal_mol"] += 0.1
    paths = []
    for name, frame in (("cohort.csv", cohort), ("fmo.csv", fmo), ("boltz.csv", boltz)):
        path = tmp_path / name
        frame.to_csv(path, index=False)
        paths.append(path)
    with pytest.raises(ValueError, match="experimental ΔG mismatch"):
        ee.evaluate_exact_fmo87(*paths, bootstrap_repetitions=5)


def test_exact_evaluation_takes_clogp_from_fmo_when_cohort_omits_it(tmp_path, monkeypatch):
    monkeypatch.setattr(ee, "EXPECTED_COUNTS", {"A": 4, "B": 4})
    cohort, fmo, boltz = _frames()
    cohort = cohort.drop(columns=["rdkit_crippen_clogp"])
    paths = []
    for name, frame in (("cohort.csv", cohort), ("fmo.csv", fmo), ("boltz.csv", boltz)):
        path = tmp_path / name
        frame.to_csv(path, index=False)
        paths.append(path)
    result = ee.evaluate_exact_fmo87(*paths, bootstrap_repetitions=2)
    assert result.compounds["rdkit_crippen_clogp"].notna().all()


def test_exact_evaluation_verifies_boltz_input_identity(tmp_path, monkeypatch):
    monkeypatch.setattr(ee, "EXPECTED_COUNTS", {"A": 4, "B": 4})
    cohort, fmo, boltz = _frames()
    paths = []
    for name, frame in (("cohort.csv", cohort), ("fmo.csv", fmo), ("boltz.csv", boltz)):
        path = tmp_path / name
        frame.to_csv(path, index=False)
        paths.append(path)
    records = []
    for index, row in enumerate(boltz.itertuples(index=False)):
        input_path = tmp_path / f"{row.target}_{row.ligand_id}.yaml"
        input_path.write_text("input", encoding="utf-8")
        digest = hashlib.sha256(input_path.read_bytes()).hexdigest()
        boltz.at[boltz.index[index], "input"] = str(input_path)
        boltz.at[boltz.index[index], "input_sha256"] = digest
        records.append({
            "target": row.target, "ligand_id": row.ligand_id,
            "input": str(input_path), "input_sha256": digest,
            "smiles": "C", "experimental_dg_kcal_mol": row.experimental_dg_kcal_mol,
        })
    boltz_path = tmp_path / "boltz.csv"
    boltz.to_csv(boltz_path, index=False)
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({
        "defaults": {"model": "boltz2", "seed": 0, "recycling_steps": 3,
                     "sampling_steps_affinity": 200, "diffusion_samples_affinity": 5},
        "records": records,
    }), encoding="utf-8")
    result = ee.evaluate_exact_fmo87(
        paths[0], paths[1], boltz_path, bootstrap_repetitions=2,
        boltz_input_manifest=manifest,
    )
    assert result.audit["boltz2_input_identity"]["all_input_hashes_valid"]
