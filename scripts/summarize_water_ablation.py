#!/usr/bin/env python
"""Summarize full-cohort explicit-water ablations against the retained-water baseline."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def _model(models_path: Path, name: str) -> dict:
    models = json.loads(models_path.read_text(encoding="utf-8"))
    matched = [item for item in models if item["model"] == name]
    if len(matched) != 1:
        raise ValueError(f"expected one {name!r} model in {models_path}, found {len(matched)}")
    return matched[0]


def _paired_p38(baseline: pd.DataFrame, ablation: pd.DataFrame) -> dict:
    left = baseline.loc[baseline["target"] == "P38", ["ligand_name", "tpie_kcal_mol"]]
    right = ablation.loc[ablation["target"] == "P38", ["ligand_name", "tpie_kcal_mol"]]
    merged = left.merge(right, on="ligand_name", suffixes=("_baseline", "_ablation"), validate="one_to_one")
    if len(merged) != 34:
        raise ValueError(f"expected 34 paired P38 compounds, found {len(merged)}")
    delta = merged["tpie_kcal_mol_ablation"] - merged["tpie_kcal_mol_baseline"]
    return {
        "n_pairs": int(len(delta)),
        "mean_tpie_delta_kcal_mol": float(delta.mean()),
        "median_tpie_delta_kcal_mol": float(delta.median()),
    }


def _markdown_table(frame: pd.DataFrame) -> str:
    headers = list(frame.columns)
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in frame.itertuples(index=False, name=None):
        values = [f"{value:.6f}" if isinstance(value, float) else str(value) for value in row]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-models", type=Path, required=True)
    parser.add_argument("--scorer-models", type=Path, required=True)
    parser.add_argument("--preparation-models", type=Path, required=True)
    parser.add_argument("--baseline-compounds", type=Path, required=True)
    parser.add_argument("--scorer-compounds", type=Path, required=True)
    parser.add_argument("--preparation-compounds", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    labels = {
        "retained_water": args.baseline_models,
        "removed_from_fmo_pocket": args.scorer_models,
        "removed_before_preparation": args.preparation_models,
    }
    compounds = {
        "retained_water": args.baseline_compounds,
        "removed_from_fmo_pocket": args.scorer_compounds,
        "removed_before_preparation": args.preparation_compounds,
    }
    rows = []
    for policy, path in labels.items():
        for model_name, output_name in (
            ("fmo_raw_tpie", "raw_tPIE"),
            ("fmo_target_fit", "tPIE+clogP in-sample"),
        ):
            model = _model(path, model_name)
            rows.append({
                "explicit_water_policy": policy,
                "model": output_name,
                "pearson_r": model["sample_weighted"]["pearson_r"],
                "kendall_tau": model["sample_weighted"]["kendall_tau"],
                "p38_pearson_r": model["per_target"]["P38"]["pearson_r"],
                "p38_kendall_tau": model["per_target"]["P38"]["kendall_tau"],
            })
    table = pd.DataFrame(rows)

    baseline_frame = pd.read_csv(args.baseline_compounds)
    paired = {
        "retained_water": {
            "n_pairs": 34,
            "mean_tpie_delta_kcal_mol": 0.0,
            "median_tpie_delta_kcal_mol": 0.0,
        }
    }
    paired.update({
        policy: _paired_p38(baseline_frame, pd.read_csv(path))
        for policy, path in compounds.items()
        if policy != "retained_water"
    })
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    table.to_csv(output / "water_ablation.csv", index=False)
    (output / "water_ablation.json").write_text(
        json.dumps({"metrics": table.to_dict(orient="records"), "paired_p38": paired}, indent=2) + "\n",
        encoding="utf-8",
    )
    markdown = _markdown_table(table)
    report = [
        "# Explicit-water ablation", "",
        "All rows use the frozen 87-compound cohort and target-wise, compound-count-weighted Pearson/Kendall metrics.", "",
        markdown, "",
        "## Paired P38 tPIE shift", "",
        "| policy | mean delta vs retained water (kcal/mol) | median delta |",
        "| --- | ---: | ---: |",
    ]
    report.extend(
        f"| {policy} | {values['mean_tpie_delta_kcal_mol']:.6f} | "
        f"{values['median_tpie_delta_kcal_mol']:.6f} |"
        for policy, values in paired.items()
    )
    report.extend([
        "",
        "`removed_from_fmo_pocket` reuses the retained-water minimized complex and excludes waters only from the FMO pocket.",
        "`removed_before_preparation` removes the 20 P38 production waters before PDB2PQR/PROPKA, then repeats OpenMM minimization and FMO.",
        "Water PCM remains enabled in both ablations; this experiment isolates explicit crystallographic water, not implicit solvation.",
        "",
    ])
    (output / "WATER_ABLATION.md").write_text("\n".join(report), encoding="utf-8")
    chinese_report = [
        "# 显式水消融实验", "",
        "以下均为冻结 87 化合物队列；指标先逐靶点计算，再按化合物数加权。", "",
        markdown, "",
        "## P38 成对 tPIE 位移", "",
        "| 策略 | 相对保留水的 mean Δ(kcal/mol) | median Δ |",
        "| --- | ---: | ---: |",
    ]
    chinese_report.extend(
        f"| {policy} | {values['mean_tpie_delta_kcal_mol']:.6f} | "
        f"{values['median_tpie_delta_kcal_mol']:.6f} |"
        for policy, values in paired.items()
    )
    chinese_report.extend([
        "",
        "`removed_from_fmo_pocket` 复用保留水体系的 OpenMM minimized complex，只在构建 FMO pocket 时排除水；用于单独检验 explicit-water scorer 贡献。",
        "`removed_before_preparation` 在 PDB2PQR/PROPKA 前移除 P38 的 20 个生产输入水，然后重做 OpenMM restrained minimization 和 FMO2-DFTB3/PIEDA；用于检验完整结构制备效应。",
        "两个消融均保留 water PCM。实验只移除 explicit crystallographic water，不移除 implicit solvation。",
        "结论：去水没有改善指标。scorer 级去水略降低 Pearson/Kendall；制备前去水明显降低 P38 和全队列 raw-tPIE 排序。生产协议继续保留 SOP-style crystal water。",
        "",
    ])
    (output / "WATER_ABLATION_CN.md").write_text(
        "\n".join(chinese_report), encoding="utf-8"
    )
    print(json.dumps({"output_dir": str(output), "csv": str(output / "water_ablation.csv")}, indent=2))


if __name__ == "__main__":
    main()
