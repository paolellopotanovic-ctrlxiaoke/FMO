#!/usr/bin/env python
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from sophosqm_baseline.exact_evaluation import (
    EXPECTED_COUNTS,
    _ranking_metrics,
    finite_json,
)


SOPHOSQM_SI_SHA256 = "12e0cbf2c604466ed1599ce8a031b8dc5e1a6c534eed4d0952aec78b57ee5f13"
OFFICIAL_SOPHOSQM_COEFFICIENTS = {
    "DNA_ligase": {"alpha_tpie": 0.074, "beta_clogp": -0.456, "gamma": 1.856},
    "MUP_I": {"alpha_tpie": 0.076, "beta_clogp": -2.079, "gamma": -3.437},
    "HSP90": {"alpha_tpie": 0.092, "beta_clogp": 0.045, "gamma": 0.837},
    "P38": {"alpha_tpie": 0.030, "beta_clogp": -1.107, "gamma": 0.223},
    "JAK2": {"alpha_tpie": 0.0840, "beta_clogp": -0.643, "gamma": 0.773},
    "MCL1": {"alpha_tpie": 0.0054, "beta_clogp": -0.9418, "gamma": -1.31116},
}


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fit_prediction(frame: pd.DataFrame, feature_columns: list[str], *, leave_one_out: bool) -> np.ndarray:
    target = frame["experimental_dg_kcal_mol"].to_numpy(float)
    design = np.column_stack(
        [frame[column].to_numpy(float) for column in feature_columns] + [np.ones(len(frame))]
    )
    predictions = np.full(len(frame), np.nan, dtype=float)
    if not leave_one_out:
        coefficients = np.linalg.lstsq(design, target, rcond=None)[0]
        return design @ coefficients
    for index in range(len(frame)):
        training = np.ones(len(frame), dtype=bool)
        training[index] = False
        coefficients = np.linalg.lstsq(design[training], target[training], rcond=None)[0]
        predictions[index] = design[index] @ coefficients
    return predictions


def target_coefficients(frame: pd.DataFrame, feature_columns: list[str]) -> list[dict]:
    records = []
    names = [f"coefficient_{column}" for column in feature_columns] + ["intercept"]
    for target, subset in frame.groupby("target", sort=False):
        design = np.column_stack(
            [subset[column].to_numpy(float) for column in feature_columns] + [np.ones(len(subset))]
        )
        coefficients = np.linalg.lstsq(
            design, subset["experimental_dg_kcal_mol"].to_numpy(float), rcond=None
        )[0]
        records.append({"target": target, **dict(zip(names, map(float, coefficients)))})
    return records


def markdown_table(frame: pd.DataFrame) -> str:
    headers = [str(column) for column in frame.columns]
    rows = []
    for row in frame.itertuples(index=False, name=None):
        values = []
        for value in row:
            if isinstance(value, float):
                values.append(f"{value:.6g}")
            else:
                values.append(str(value).replace("|", "\\|"))
        rows.append(values)
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    lines.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--compounds", default="runs/fmo87/final/compounds.csv")
    parser.add_argument("--sophosqm-si", default="data/references/sophosqm_si.pdf")
    parser.add_argument("--output-dir", default="runs/fmo87/final")
    args = parser.parse_args()

    compounds_path = Path(args.compounds)
    si_path = Path(args.sophosqm_si)
    output_dir = Path(args.output_dir)
    frame = pd.read_csv(compounds_path)
    observed_counts = frame.groupby("target").size().to_dict()
    if observed_counts != EXPECTED_COUNTS or len(frame) != 87:
        raise ValueError(f"expected exact 87-compound cohort, got {observed_counts}")
    required = {
        "target", "experimental_dg_kcal_mol", "rdkit_crippen_clogp", "tpie_kcal_mol",
        "fmo_target_fit_dg_kcal_mol", "fmo_loo_dg_kcal_mol",
    }
    if missing := required - set(frame.columns):
        raise ValueError(f"missing required columns: {sorted(missing)}")
    if not np.isfinite(frame[list(required - {"target"})].to_numpy(float)).all():
        raise ValueError("calibration ablation requires finite scores and labels")
    if not si_path.exists() or file_sha256(si_path) != SOPHOSQM_SI_SHA256:
        raise ValueError(f"SophosQM SI missing or SHA-256 mismatch: {si_path}")

    frame = frame.rename(columns={"rdkit_crippen_clogp": "clogp"})
    models = {
        "raw_tpie": "tpie_kcal_mol",
        "tPIE_only_in_sample": "tpie_only_in_sample",
        "clogP_only_in_sample": "clogp_only_in_sample",
        "tPIE_clogP_in_sample": "tpie_clogp_in_sample",
        "tPIE_only_LOO": "tpie_only_loo",
        "clogP_only_LOO": "clogp_only_loo",
        "tPIE_clogP_LOO": "tpie_clogp_loo",
    }
    feature_sets = {
        "tpie_only": ["tpie_kcal_mol"],
        "clogp_only": ["clogp"],
        "tpie_clogp": ["tpie_kcal_mol", "clogp"],
    }
    for label, columns in feature_sets.items():
        frame[f"{label}_in_sample"] = np.concatenate([
            fit_prediction(subset, columns, leave_one_out=False)
            for _, subset in frame.groupby("target", sort=False)
        ])
        frame[f"{label}_loo"] = np.concatenate([
            fit_prediction(subset, columns, leave_one_out=True)
            for _, subset in frame.groupby("target", sort=False)
        ])

    official_p38 = OFFICIAL_SOPHOSQM_COEFFICIENTS["P38"]
    frame["official_p38_coefficients"] = (
        official_p38["alpha_tpie"] * frame["tpie_kcal_mol"]
        + official_p38["beta_clogp"] * frame["clogp"]
        + official_p38["gamma"]
    )
    frame["official_p38_coefficients_p38_only"] = frame["official_p38_coefficients"].where(
        frame["target"] == "P38", np.nan
    )

    metric_rows = []
    for model, column in models.items():
        metrics = _ranking_metrics(frame, column, "experimental_dg_kcal_mol")
        metric_rows.append({
            "model": model,
            "pearson_r": metrics["sample_weighted"]["pearson_r"],
            "kendall_tau": metrics["sample_weighted"]["kendall_tau"],
            "pairwise_mae_kcal_mol": metrics["sample_weighted"]["pairwise_mae_kcal_mol"],
            "mae_noncentered_kcal_mol": metrics["sample_weighted"]["mae_noncentered_kcal_mol"],
            "mae_centered_kcal_mol": metrics["sample_weighted"]["mae_centered_kcal_mol"],
        })
    metrics_frame = pd.DataFrame(metric_rows)

    p38 = frame.loc[frame["target"] == "P38"]
    p38_transfer = pd.DataFrame([
        {
            "model": model,
            "pearson_r": float(np.corrcoef(p38[column], p38["experimental_dg_kcal_mol"])[0, 1]),
        }
        for model, column in {
            "raw_tpie": "tpie_kcal_mol",
            "local_tPIE_clogP_in_sample": "tpie_clogp_in_sample",
            "local_tPIE_clogP_LOO": "tpie_clogp_loo",
            "official_p38_alpha_beta_gamma": "official_p38_coefficients",
        }.items()
    ])

    coefficients = {
        "official_sophosqm_si": OFFICIAL_SOPHOSQM_COEFFICIENTS,
        "local_tpie_only": target_coefficients(frame, feature_sets["tpie_only"]),
        "local_clogp_only": target_coefficients(frame, feature_sets["clogp_only"]),
        "local_tpie_clogp": target_coefficients(frame, feature_sets["tpie_clogp"]),
    }
    result = {
        "audit": {
            "n_compounds": len(frame),
            "counts_by_target": observed_counts,
            "sophosqm_si_sha256": SOPHOSQM_SI_SHA256,
            "clogp_implementation": "RDKit Crippen, not original VolSurf",
            "fmo_score": "GAMESS FMO2-DFTB3/3OB-3-1 water-PCM tPIE, not original FMO2-MP2/6-31G*",
        },
        "metrics": metrics_frame.to_dict(orient="records"),
        "p38_transfer_sensitivity": p38_transfer.to_dict(orient="records"),
        "coefficients": coefficients,
        "interpretation": {
            "ranking_calibration": "A positive tPIE-only affine fit does not change within-target ranking; only the clogP term can change it.",
            "prospective_test": "LOO metrics evaluate whether fitted coefficients generalize; they are not paper Table 12 comparisons.",
            "official_transfer": "Official SI coefficients are target/series-specific; P38 transfer is only sensitivity analysis because cohorts and methods differ.",
        },
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    metrics_frame.to_csv(output_dir / "calibration_ablation.csv", index=False)
    frame.to_csv(output_dir / "calibration_ablation_predictions.csv", index=False)
    (output_dir / "calibration_ablation.json").write_text(
        json.dumps(finite_json(result), indent=2) + "\n", encoding="utf-8"
    )
    report = [
        "# SophosQM 校准消融",
        "",
        "## 排名与误差",
        "",
        markdown_table(metrics_frame),
        "",
        "## p38 官方系数迁移敏感性",
        "",
        markdown_table(p38_transfer),
        "",
        "## 解释",
        "",
        "- `tPIE-only` 的正向仿射校准只改变绝对 kcal/mol 误差，不改变靶点内 Pearson/Kendall 排名。",
        "- `tPIE+clogP` in-sample 提升不能当作外推能力；LOO 才用于检查拟合系数是否稳定。",
        "- 官方 SI 系数来自 SophosQM 六个原始 series；除 p38 目标名重合外，Boltz2 cohort 的配体、pose、几何和 QM 方法均不同，因此只做迁移敏感性。",
        "- 本地 clogP 是 RDKit Crippen，原文是 VolSurf；本地 FMO 是 DFTB3，原文 SophosQM 主模型是 MP2/6-31G*。",
    ]
    (output_dir / "CALIBRATION_ABLATION_CN.md").write_text("\n".join(report) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
