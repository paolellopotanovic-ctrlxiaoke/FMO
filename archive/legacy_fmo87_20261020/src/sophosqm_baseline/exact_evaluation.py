"""Strict evaluation on the frozen Boltz-2 four-target FMO cohort."""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


EXPECTED_COUNTS = {"CDK2": 16, "JNK1": 21, "P38": 34, "TYK2": 16}
PAPER_BOLTZ2 = {
    "pearson_r": 0.66,
    "kendall_tau": 0.48,
    "pairwise_mae_kcal_mol": 0.85,
    "mae_noncentered_kcal_mol": 0.75,
    "mae_centered_kcal_mol": 0.59,
    "percent_within_1_noncentered_kcal_mol": 0.69,
    "percent_within_1_centered_kcal_mol": 0.83,
    "percent_within_2_noncentered_kcal_mol": 0.97,
    "percent_within_2_centered_kcal_mol": 0.98,
}
PAPER_FMO = {"pearson_r": 0.55, "kendall_tau": 0.38}
RT_LN_10_KCAL_MOL = 0.00198720425864083 * 298.15 * math.log(10.0)


@dataclass(frozen=True)
class ExactEvaluation:
    compounds: pd.DataFrame
    models: pd.DataFrame
    metrics: dict
    audit: dict


def _key(target: object, ligand_id: object) -> str:
    return f"{str(target).upper()}|{str(ligand_id)}"


def _identity(frame: pd.DataFrame, target_column: str, ligand_column: str) -> pd.Index:
    keys = [_key(target, ligand) for target, ligand in zip(frame[target_column], frame[ligand_column])]
    return pd.Index(keys, name="cohort_key")


def _audit_cohort(frame: pd.DataFrame, label: str, *,
                  expected_keys: set[str],
                  status_column: str | None = None) -> dict:
    if frame.index.duplicated().any():
        duplicates = sorted(frame.index[frame.index.duplicated()].unique().tolist())
        raise ValueError(f"duplicate {label} cohort identities: {duplicates}")
    observed_keys = set(frame.index.tolist())
    missing = sorted(expected_keys - observed_keys)
    unexpected = sorted(observed_keys - expected_keys)
    counts = frame.groupby(level=0).size().groupby(lambda key: key.split("|", 1)[0]).sum().to_dict()
    audit = {
        "label": label,
        "n_records": len(frame),
        "counts_by_target": counts,
        "missing_keys": missing,
        "unexpected_keys": unexpected,
        "complete": not missing and not unexpected and counts == EXPECTED_COUNTS,
    }
    if status_column is not None and status_column in frame:
        audit["status_counts"] = frame[status_column].value_counts(dropna=False).astype(int).to_dict()
        failures = frame.loc[frame[status_column] != "SUCCEEDED"] if label == "FMO" else frame.loc[~frame[status_column].fillna("").str.startswith("SUCCEEDED")]
        audit["failure_keys"] = failures.index.tolist()
        audit["all_required_succeeded"] = len(failures) == 0
    if not audit["complete"]:
        raise ValueError(f"incomplete {label} cohort: {audit}")
    return audit


def _pearson(predicted: np.ndarray, observed: np.ndarray) -> float:
    predicted = np.asarray(predicted, dtype=float)
    observed = np.asarray(observed, dtype=float)
    if len(predicted) < 2 or np.std(predicted) == 0 or np.std(observed) == 0:
        return math.nan
    return float(np.corrcoef(predicted, observed)[0, 1])


def _kendall(predicted: np.ndarray, observed: np.ndarray) -> float:
    predicted = np.asarray(predicted, dtype=float)
    observed = np.asarray(observed, dtype=float)
    if len(predicted) < 2:
        return math.nan
    concordant = discordant = tied_prediction = tied_observed = 0
    for left in range(len(predicted)):
        for right in range(left + 1, len(predicted)):
            first = np.sign(predicted[left] - predicted[right])
            second = np.sign(observed[left] - observed[right])
            if first == 0 or second == 0:
                if first == 0 and second != 0:
                    tied_prediction += 1
                if second == 0 and first != 0:
                    tied_observed += 1
                if first == 0 and second == 0:
                    tied_prediction += 1
                    tied_observed += 1
            elif first == second:
                concordant += 1
            else:
                discordant += 1
    denominator = math.sqrt((concordant + discordant + tied_prediction) * (concordant + discordant + tied_observed))
    return float((concordant - discordant) / denominator) if denominator else math.nan


def _ranking_metrics(frame: pd.DataFrame, predicted: str, observed: str) -> dict:
    per_target = {}
    for target in EXPECTED_COUNTS:
        subset = frame.loc[frame["target"] == target]
        predicted_values = subset[predicted].to_numpy(float)
        observed_values = subset[observed].to_numpy(float)
        residuals = predicted_values - observed_values
        centered_residuals = residuals - residuals.mean()
        if len(predicted_values) > 1:
            pairwise_errors = np.abs(
                (predicted_values[:, None] - predicted_values[None, :])
                - (observed_values[:, None] - observed_values[None, :])
            )
            pairwise_mae = float(pairwise_errors[np.triu_indices(len(predicted_values), 1)].mean())
        else:
            pairwise_mae = math.nan
        per_target[target] = {
            "n": int(len(subset)),
            "pearson_r": _pearson(predicted_values, observed_values),
            "kendall_tau": _kendall(predicted_values, observed_values),
            "pairwise_mae_kcal_mol": pairwise_mae,
            "mae_noncentered_kcal_mol": float(np.abs(residuals).mean()) if len(residuals) else math.nan,
            "mae_centered_kcal_mol": float(np.abs(centered_residuals).mean()) if len(residuals) else math.nan,
            "percent_within_1_noncentered_kcal_mol": float((np.abs(residuals) <= 1.0).mean()) if len(residuals) else math.nan,
            "percent_within_1_centered_kcal_mol": float((np.abs(centered_residuals) <= 1.0).mean()) if len(residuals) else math.nan,
            "percent_within_2_noncentered_kcal_mol": float((np.abs(residuals) <= 2.0).mean()) if len(residuals) else math.nan,
            "percent_within_2_centered_kcal_mol": float((np.abs(centered_residuals) <= 2.0).mean()) if len(residuals) else math.nan,
        }
    weighted = {}
    valid_counts = {target: item["n"] for target, item in per_target.items()
                    if math.isfinite(item["pearson_r"]) and math.isfinite(item["kendall_tau"])}
    metric_names = (
        "pearson_r",
        "kendall_tau",
        "pairwise_mae_kcal_mol",
        "mae_noncentered_kcal_mol",
        "mae_centered_kcal_mol",
        "percent_within_1_noncentered_kcal_mol",
        "percent_within_1_centered_kcal_mol",
        "percent_within_2_noncentered_kcal_mol",
        "percent_within_2_centered_kcal_mol",
    )
    for metric in metric_names:
        weighted[metric] = (
            sum(valid_counts[target] * per_target[target][metric] for target in valid_counts) / sum(valid_counts.values())
            if valid_counts else math.nan
        )
    return {"per_target": per_target, "sample_weighted": weighted,
            "n_compounds": int(sum(valid_counts.values()))}


def _fit_sophosqm(compounds: pd.DataFrame) -> pd.DataFrame:
    predictions = []
    coefficients = []
    for target, subset in compounds.groupby("target", sort=False):
        subset = subset.sort_index()
        features = subset[["tpie_kcal_mol", "rdkit_crippen_clogp"]].to_numpy(float)
        observed = subset["experimental_dg_kcal_mol"].to_numpy(float)
        design = np.column_stack([features, np.ones(len(subset))])
        if np.linalg.matrix_rank(design) < design.shape[1]:
            raise ValueError(f"rank-deficient SophosQM fit for {target}")
        fitted_coefficients = np.linalg.lstsq(design, observed, rcond=None)[0]
        fitted = design @ fitted_coefficients
        for index, (cohort_key, row) in enumerate(subset.iterrows()):
            train = np.arange(len(subset)) != index
            if np.linalg.matrix_rank(design[train]) < design.shape[1]:
                loo = math.nan
            else:
                loo_coefficients = np.linalg.lstsq(design[train], observed[train], rcond=None)[0]
                loo = float(design[index] @ loo_coefficients)
            predictions.append({
                "cohort_key": cohort_key,
                "fmo_target_fit_dg_kcal_mol": float(fitted[index]),
                "fmo_loo_dg_kcal_mol": loo,
            })
        coefficients.append({
            "target": target,
            "alpha_tpie": float(fitted_coefficients[0]),
            "beta_clogp": float(fitted_coefficients[1]),
            "gamma": float(fitted_coefficients[2]),
            "n_training": int(len(subset)),
        })
    output = compounds.join(pd.DataFrame(predictions).set_index("cohort_key"))
    return output, pd.DataFrame(coefficients)


def _bootstrap_interval(frame: pd.DataFrame, predicted: str, seed: int,
                        repetitions: int = 10000) -> dict:
    rng = np.random.default_rng(seed)
    values = {metric: [] for metric in ("pearson_r", "kendall_tau")}
    grouped = [subset.reset_index(drop=True) for _, subset in frame.groupby("target", sort=False)]
    for _ in range(repetitions):
        sampled = pd.concat([subset.loc[rng.integers(0, len(subset), len(subset))] for subset in grouped])
        metrics = _ranking_metrics(sampled, predicted, "experimental_dg_kcal_mol")["sample_weighted"]
        for metric in values:
            if math.isfinite(metrics[metric]):
                values[metric].append(metrics[metric])
    return {metric: {
        "low_95": float(np.percentile(items, 2.5)) if items else None,
        "high_95": float(np.percentile(items, 97.5)) if items else None,
        "repetitions": len(items),
    } for metric, items in values.items()}


def _audit_boltz_input_identity(manifest_path: str | Path, cohort: pd.DataFrame,
                                 boltz: pd.DataFrame) -> dict:
    manifest = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    records = manifest.get("records", [])
    if not isinstance(records, list) or not records:
        raise ValueError("Boltz2 input manifest contains no records")
    expected_defaults = {
        "model": "boltz2", "seed": 0, "recycling_steps": 3,
        "sampling_steps_affinity": 200, "diffusion_samples_affinity": 5,
    }
    if any(manifest.get("defaults", {}).get(key) != value for key, value in expected_defaults.items()):
        raise ValueError(f"Boltz2 input manifest defaults differ from production protocol: {manifest.get('defaults')}")

    def sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        return digest.hexdigest()

    identities = {}
    for record in records:
        key = _key(record["target"], record["ligand_id"])
        if key in identities:
            raise ValueError(f"duplicate Boltz2 input identity: {key}")
        identities[key] = record
    if set(identities) != set(cohort.index):
        raise ValueError("Boltz2 input manifest identities differ from the frozen cohort")

    checked = 0
    for key, record in identities.items():
        input_path = Path(record["input"])
        if not input_path.is_file():
            raise FileNotFoundError(input_path)
        actual_sha256 = sha256(input_path)
        if actual_sha256 != record["input_sha256"]:
            raise ValueError(f"Boltz2 input file hash mismatch for {key}")
        if str(boltz.loc[key, "input"]) != str(input_path):
            raise ValueError(f"Boltz2 prediction input path mismatch for {key}")
        if boltz.loc[key, "input_sha256"] != record["input_sha256"]:
            raise ValueError(f"Boltz2 prediction input hash mismatch for {key}")
        if record["smiles"] != cohort.loc[key, "canonical_isomeric_smiles"]:
            raise ValueError(f"Boltz2 canonical SMILES mismatch for {key}")
        if not math.isclose(float(record["experimental_dg_kcal_mol"]),
                            float(cohort.loc[key, "experimental_dg_kcal_mol"])):
            raise ValueError(f"Boltz2 experimental ΔG mismatch for {key}")
        checked += 1
    return {
        "manifest": str(manifest_path), "n_records": checked,
        "all_input_hashes_valid": True,
        "all_canonical_smiles_match_cohort": True,
    }


def evaluate_exact_fmo87(cohort_csv: str | Path, fmo_csv: str | Path,
                         boltz_csv: str | Path, *,
                         bootstrap_seed: int = 20260925,
                         bootstrap_repetitions: int = 10000,
                         boltz_input_manifest: str | Path | None = None) -> ExactEvaluation:
    cohort = pd.read_csv(cohort_csv, dtype={"ligand_name": str})
    fmo = pd.read_csv(fmo_csv, dtype={"ligand_name": str})
    boltz = pd.read_csv(boltz_csv, dtype={"ligand_id": str})
    if len(cohort) != sum(EXPECTED_COUNTS.values()) or cohort["target"].value_counts().to_dict() != EXPECTED_COUNTS:
        raise ValueError("cohort CSV is not the frozen four-target 87-compound cohort")
    if cohort.duplicated(["target", "ligand_name"]).any():
        raise ValueError("duplicate cohort identity")

    cohort = cohort.copy()
    cohort["cohort_key"] = [_key(target, ligand) for target, ligand in zip(cohort["target"], cohort["ligand_name"])]
    cohort = cohort.set_index("cohort_key")
    expected_keys = set(cohort.index.tolist())
    audit_cohort = _audit_cohort(cohort, "cohort", expected_keys=expected_keys)

    fmo = fmo.copy()
    fmo["cohort_key"] = [_key(target, ligand) for target, ligand in zip(fmo["target"], fmo["ligand_name"])]
    fmo = fmo.set_index("cohort_key")
    audit_fmo = _audit_cohort(fmo, "FMO", expected_keys=expected_keys, status_column="status")
    if not audit_fmo.get("all_required_succeeded", False):
        raise ValueError(f"FMO contains failed records: {audit_fmo['failure_keys']}")
    if not np.isfinite(fmo["tpie_kcal_mol"].to_numpy(float)).all():
        raise ValueError("nonfinite FMO tPIE")

    boltz = boltz.copy()
    boltz["cohort_key"] = [_key(target, ligand) for target, ligand in zip(boltz["target"], boltz["ligand_id"])]
    boltz = boltz.set_index("cohort_key")
    audit_boltz = _audit_cohort(boltz, "Boltz2", expected_keys=expected_keys, status_column="status")
    if not audit_boltz.get("all_required_succeeded", False):
        raise ValueError(f"Boltz2 contains failed records: {audit_boltz['failure_keys']}")
    if not np.isfinite(boltz["affinity_pred_value"].to_numpy(float)).all():
        raise ValueError("nonfinite Boltz2 affinity")
    boltz_input_audit = (
        _audit_boltz_input_identity(boltz_input_manifest, cohort, boltz)
        if boltz_input_manifest is not None else None
    )

    compounds = cohort[["target", "ligand_name", "canonical_isomeric_smiles",
                        "experimental_dg_kcal_mol"]].copy()
    if "rdkit_crippen_clogp" in cohort:
        compounds["rdkit_crippen_clogp"] = cohort["rdkit_crippen_clogp"]
    else:
        compounds["rdkit_crippen_clogp"] = fmo["rdkit_crippen_clogp"].reindex(compounds.index)
    compounds["tpie_kcal_mol"] = fmo["tpie_kcal_mol"].reindex(compounds.index)
    compounds["fmo_status"] = fmo["status"].reindex(compounds.index)
    compounds["boltz_status"] = boltz["status"].reindex(compounds.index)
    compounds["boltz_affinity_pred_value"] = boltz["affinity_pred_value"].reindex(compounds.index)
    compounds["boltz_dg_kcal_mol"] = RT_LN_10_KCAL_MOL * (compounds["boltz_affinity_pred_value"] - 6.0)

    observed = compounds["experimental_dg_kcal_mol"].to_numpy(float)
    if not np.isfinite(observed).all():
        raise ValueError("nonfinite experimental ΔG")
    if not np.allclose(fmo["experimental_dg_kcal_mol"].reindex(compounds.index).to_numpy(float), observed):
        raise ValueError("experimental ΔG mismatch between cohort and FMO")
    if not np.allclose(boltz["experimental_dg_kcal_mol"].reindex(compounds.index).to_numpy(float), observed):
        raise ValueError("experimental ΔG mismatch between cohort and Boltz2")
    if "canonical_isomeric_smiles" in fmo and not (
            fmo["canonical_isomeric_smiles"].reindex(compounds.index) ==
            compounds["canonical_isomeric_smiles"]).all():
        raise ValueError("canonical SMILES mismatch between cohort and FMO")
    if "formal_charge" in fmo and not (
            fmo["formal_charge"].reindex(compounds.index).astype(int) == 0).all():
        raise ValueError("FMO contains a charged ligand in the neutral cohort")

    compounds, fit_coefficients = _fit_sophosqm(compounds)
    model_rows = []
    for label, column in (
        ("boltz2_local", "boltz_dg_kcal_mol"),
        ("fmo_raw_tpie", "tpie_kcal_mol"),
        ("fmo_target_fit", "fmo_target_fit_dg_kcal_mol"),
        ("fmo_loo", "fmo_loo_dg_kcal_mol"),
    ):
        metrics = _ranking_metrics(compounds, column, "experimental_dg_kcal_mol")
        model_rows.append({"model": label, **metrics})
    models = pd.DataFrame(model_rows)

    comparison_rows = []
    for label, row in models.set_index("model").iterrows():
        if label == "boltz2_local":
            reference = PAPER_BOLTZ2
        elif label == "fmo_raw_tpie":
            reference = PAPER_FMO
        else:
            continue
        comparison = {
            "model": label,
            "local_pearson_r": row["sample_weighted"]["pearson_r"],
            "local_kendall_tau": row["sample_weighted"]["kendall_tau"],
            "paper_pearson_r": reference["pearson_r"],
            "paper_kendall_tau": reference["kendall_tau"],
            "pearson_delta": row["sample_weighted"]["pearson_r"] - reference["pearson_r"],
            "kendall_delta": row["sample_weighted"]["kendall_tau"] - reference["kendall_tau"],
        }
        for metric in (
            "pairwise_mae_kcal_mol",
            "mae_noncentered_kcal_mol",
            "mae_centered_kcal_mol",
            "percent_within_1_noncentered_kcal_mol",
            "percent_within_1_centered_kcal_mol",
            "percent_within_2_noncentered_kcal_mol",
            "percent_within_2_centered_kcal_mol",
        ):
            comparison[f"local_{metric}"] = row["sample_weighted"][metric]
            comparison[f"paper_{metric}"] = reference.get(metric)
        comparison_rows.append(comparison)
    comparison = pd.DataFrame(comparison_rows)
    bootstrap = {
        label: _bootstrap_interval(compounds, column, bootstrap_seed + offset, bootstrap_repetitions)
        for offset, (label, column) in enumerate((
            ("boltz2_local", "boltz_dg_kcal_mol"),
            ("fmo_raw_tpie", "tpie_kcal_mol"),
            ("fmo_target_fit", "fmo_target_fit_dg_kcal_mol"),
            ("fmo_loo", "fmo_loo_dg_kcal_mol"),
        ))
    }
    metrics = {
        "paper_reference": {"Boltz-2": PAPER_BOLTZ2, "FMO": PAPER_FMO},
        "sample_weighted_target_average": models.set_index("model")["sample_weighted"].apply(pd.Series).to_dict(orient="index"),
        "per_target": models.set_index("model")["per_target"].to_dict(),
        "paper_comparison": comparison.to_dict(orient="records"),
        "bootstrap_95": bootstrap,
        "sophosqm_target_fit_coefficients": fit_coefficients.to_dict(orient="records"),
        "affinity_conversion": {
            "equation": "RT ln(10) * (affinity_pred_value - 6)",
            "RT_ln10_kcal_mol": RT_LN_10_KCAL_MOL,
            "temperature_k": 298.15,
        },
    }
    audit = {
        "cohort": audit_cohort,
        "fmo": audit_fmo,
        "boltz2": audit_boltz,
        "boltz2_input_identity": boltz_input_audit,
        "strict_87_complete": bool(audit_cohort["complete"] and audit_fmo["complete"] and audit_boltz["complete"]),
        "all_required_runs_succeeded": bool(audit_fmo.get("all_required_succeeded") and audit_boltz.get("all_required_succeeded")),
        "model_interpretation": {
            "fmo_raw_tpie": "unfitted enthalpic ranking score",
            "fmo_target_fit": "per-target in-sample SophosQM-style α*tPIE+β*clogP+γ fit",
            "fmo_loo": "per-target leave-one-variant-out prospective SophosQM-style score",
            "boltz2_local": "direct zero-shot Boltz-2 affinity prediction",
        },
    }
    return ExactEvaluation(compounds, models, metrics, audit)


def finite_json(value):
    if isinstance(value, dict):
        return {key: finite_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [finite_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    return value
