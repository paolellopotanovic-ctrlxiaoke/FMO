"""Compare an unfitted positive molecular-weight score on complete cohorts."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors

from build_model_comparison import COHORTS, check_metrics, digest, format_metrics, metrics, read_csv


def build(workspace: Path, root: Path) -> dict:
    workspace = workspace.resolve()
    comparison_path = root / "docs/MODEL_COMPARISON.json"
    comparison = json.loads(comparison_path.read_text())
    evaluation_path = workspace / "FMOtpiefinal_HIV1lab/runs/all_accelerated_robust_20261020/evaluation.json"
    evaluation = json.loads(evaluation_path.read_text())
    sources = {str(comparison_path): digest(comparison_path), str(evaluation_path): digest(evaluation_path)}
    original_path = workspace / "data/fmo87_cohort.csv"
    original = read_csv(original_path)
    sources[str(original_path)] = digest(original_path)
    output = {"date": "2026-10-21", "score": "+RDKit Descriptors.MolWt, larger predicts stronger binding",
              "fitted": False, "new_QM": False, "groups": {}, "sources_sha256": sources}
    cases = []
    for target, (expected, sdf_stem, fep_file) in COHORTS.items():
        if target in ("CDK2", "JNK1", "P38", "TYK2"):
            selected = [row for row in original if row["target"] == target]
            molecules = {row["ligand_name"]: Chem.MolFromSmiles(row["canonical_isomeric_smiles"])
                         for row in selected}
            truth = {row["ligand_name"]: -float(row["experimental_dg_kcal_mol"]) for row in selected}
        else:
            source = workspace / "vendor/public_binding_free_energy_benchmark/fep_benchmark_inputs/structure_inputs" / f"{sdf_stem}_ligands.sdf"
            sources[str(source)] = digest(source)
            loaded = list(Chem.SDMolSupplier(str(source), removeHs=False))
            if any(molecule is None for molecule in loaded):
                raise ValueError(f"invalid molecule: {target}")
            molecules = {molecule.GetProp("_Name").strip(): molecule for molecule in loaded}
            if len(molecules) != len(loaded):
                raise ValueError(f"duplicate molecule: {target}")
            truth = {name: -float(molecule.GetProp("r_exp_dg")) for name, molecule in molecules.items()}
        if target == "Galectin3":
            source = workspace / "FMOtpiefinal_HIV1lab/datasets/evaluation_only/Galectin3.csv"
            sources[str(source)] = digest(source)
            truth = {row["ligand"]: float(row["experimental_score_minus_dG"]) for row in read_csv(source)}
        elif target == "FXa":
            source = workspace / "vendor/public_binding_free_energy_benchmark/21_4_results/ligand_predictions" / fep_file
            sources[str(source)] = digest(source)
            truth = {row["Ligand name"].strip(): -float(row["Exp. dG (kcal/mol)"]) for row in read_csv(source)}
        if len(truth) != expected or set(truth) != set(molecules) or any(molecule is None for molecule in molecules.values()):
            raise ValueError(f"incomplete molecular-weight cohort: {target}")
        verification_path = workspace / f"FMOtpiefinal_HIV1lab/runs/all_accelerated_robust_20261020/{target}/independent_verification.json"
        sources[str(verification_path)] = digest(verification_path)
        records = json.loads(verification_path.read_text())["records"]
        predictions = {row["ligand"]: row["raw_score"] for row in records}
        if len(predictions) != len(records) or set(predictions) != set(truth):
            raise ValueError(f"FMO cohort mismatch: {target}")
        names = sorted(truth)
        observed = [truth[name] for name in names]
        check_metrics(metrics(observed, [predictions[name] for name in names]), comparison["groups"][target]["FMO_raw"])
        weights = [float(Descriptors.MolWt(molecules[name])) for name in names]
        output["groups"][target] = {"n": expected, "MW_only": metrics(observed, weights),
                                    "FMO_raw": comparison["groups"][target]["FMO_raw"]}
        cases.extend({"cohort": target, "id": name, "label": truth[name], "mw": weight}
                     for name, weight in zip(names, weights))
    source = workspace / "runs/openbind/score_distribution_audit/paired_compounds.csv"
    sources[str(source)] = digest(source)
    paired = read_csv(source)
    latest = {str(row["compound_index"]): row for row in evaluation["groups"]["OpenBind"]["compound_predictions"]}
    if len(paired) != 490 or len({row["compound_index"] for row in paired}) != 490 or set(latest) != {row["compound_index"] for row in paired}:
        raise ValueError("OpenBind compound identities differ")
    for row in paired:
        weight = float(Descriptors.MolWt(Chem.MolFromSmiles(row["compound_smiles"])))
        if abs(weight - float(row["molecular_weight"])) > 0.05 or abs(float(row["experimental_pKD"]) - latest[row["compound_index"]]["experimental_pKD"]) > 1e-10:
            raise ValueError("OpenBind molecular weight or label differs")
        cases.append({"cohort": "OpenBind", "id": row["compound_index"], "label": float(row["experimental_pKD"]), "mw": weight})
    members = [row for row in cases if row["cohort"] == "OpenBind"]
    output["groups"]["OpenBind"] = {"n": 490, "MW_only": metrics([row["label"] for row in members], [row["mw"] for row in members]),
                                      "FMO_raw": comparison["groups"]["OpenBind"]["FMO_raw"]}
    output["original87_target_size_weighted"] = {model: {metric: sum(output["groups"][target]["n"] * output["groups"][target][model][metric]
        for target in ("CDK2", "JNK1", "P38", "TYK2")) / 87 for metric in ("pearson", "spearman", "kendall")}
        for model in ("MW_only", "FMO_raw")}
    output["cases"] = cases
    (root / "docs/MOLECULAR_WEIGHT_COMPARISON.json").write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n")
    lines = ["# 分子量无拟合基线与当前tPIE", "", "核算日期：2026-10-21。分数为+MolWt，即分子量越大预测越强；不拟合、不翻转方向、不加入FMO主流程。",
             "全部按当前完整队列身份核对。每格为Pearson / Spearman / Kendall；FMO为未软压缩tPIE的有利方向分数。", "",
             "| 队列 | n | 分子量基线 | 当前tPIE |", "|---|---:|---|---|"]
    for target, group in output["groups"].items():
        lines.append(f"| {target} | {group['n']} | {format_metrics(group['MW_only'])} | {format_metrics(group['FMO_raw'])} |")
    weighted = output["original87_target_size_weighted"]
    lines.extend([f"| 原87：靶点样本数加权 | 87 | {format_metrics(weighted['MW_only'])} | {format_metrics(weighted['FMO_raw'])} |", "",
                  "原87汇总是四个靶点的相关系数按样本数加权，不是混合87个标签计算单一相关性。OpenBind按490化合物，不重复计算639晶体。",
                  "分子量基线不是FMO模型，也不是控制分子量后的增量分析；某队列分子量相关性较高不证明因果或跨系列泛化。",
                  "同分子量异构体会产生并列。需要进一步区分尺寸效应和额外相互作用信息时，应另做同靶点分层/分箱或配对诊断。",
                  "本次未运行新QM，未调整分子量权重或按测试集选分数。精确值、逐化合物分子量和来源哈希见同名JSON。"])
    (root / "docs/MOLECULAR_WEIGHT_COMPARISON.zh-CN.md").write_text("\n".join(lines) + "\n")
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    arguments = parser.parse_args()
    report = build(arguments.workspace, Path(__file__).resolve().parents[1])
    print(json.dumps({target: values["MW_only"] for target, values in report["groups"].items()}, indent=2))
