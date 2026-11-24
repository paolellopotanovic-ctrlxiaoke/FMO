"""Publish completed controls and explicitly mark unfinished QM geometry work."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

from run_cdk8_solution_cycle import digest, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    runs = args.workspace.resolve() / "FMOtpiefinal_HIV1lab/runs"
    reference = runs / "CDK8_solution_cycle_consistent_20261021"
    bound = runs / "CDK8_bound_qm_20261021"
    input_pose = runs / "CDK8_input_pose_control_20261021/evaluation.json"
    gradient = bound / "gradient_check/evaluation.json"
    previous = json.loads((reference / "final_evaluation.json").read_text())
    sampled = json.loads((reference / "solution_sampling4/final_evaluation.json").read_text())
    raw = json.loads(input_pose.read_text())
    differences = json.loads(gradient.read_text())
    mm_path = bound / "mm_internal/evaluation.json"
    mm = json.loads(mm_path.read_text()) if mm_path.exists() else None
    amides = json.loads((reference / "amide_geometry_diagnostic.json").read_text())
    previous_rows = {row["ligand"]: row for row in previous["rows"]}
    geometry = {}
    optimization_progress = {}
    for name in previous_rows:
        path = bound / name / "summary.json"
        geometry[name] = json.loads(path.read_text()).get("status", "UNKNOWN") if path.exists() else "RUNNING_OR_QUEUED"
        log = bound / name / "optimization/gamess.log"
        if log.exists():
            last_step = None
            with log.open(errors="replace") as handle:
                for line in handle:
                    if "NSERCH:" in line:
                        last_step = line.strip()
            optimization_progress[name] = last_step
    complete = all(status == "SUCCEEDED" for status in geometry.values())
    available_metrics = {arm: previous["metrics"][arm] for arm in (
        "tpie_pcm", "three_state_pcm_bound", "three_state_pcm_free_min", "tpie_plus_ligand_strain")}
    available_metrics.update(sampled["metrics"])
    available_metrics.update(raw["metrics"])
    final_geometry = bound / "geometry_evidence.json"
    if complete and final_geometry.exists():
        available_metrics.update(json.loads(final_geometry.read_text())["metrics"])
    goal_achieved = any(value["spearman"] >= 0.6 for arm,value in available_metrics.items() if arm != "fep_archived_dg")
    truth = [row["experimental_dg"] for row in previous["rows"]]
    sources = [reference / "final_evaluation.json", reference / "solution_sampling4/final_evaluation.json",
        input_pose, gradient, reference / "amide_geometry_diagnostic.json"]
    fep_path = args.workspace.resolve() / "vendor/public_binding_free_energy_benchmark/21_4_results/ligand_predictions/misc/cdk8_koehler_213_out.csv"
    with fep_path.open() as handle:
        fep = {row["Ligand name"]: row for row in csv.DictReader(handle)}
    fep_rows = [{"ligand": name, "experimental_dg": float(row["Exp. dG (kcal/mol)"]),
        "predicted_dg": float(row["Pred. dG (kcal/mol)"]),
        "prediction_minus_experiment": float(row["Pred. dG (kcal/mol)"]) - float(row["Exp. dG (kcal/mol)"])} for name,row in fep.items()]
    sources.append(fep_path)
    if mm is not None:
        sources.append(mm_path)
    comparison_rows = [{"ligand": name, "experimental_dg": row["experimental_dg"],
        "tpie": row["tpie_pcm"], "fep": float(fep[name]["Pred. dG (kcal/mol)"])} for name,row in previous_rows.items()]
    ordered_truth = rankdata([row["experimental_dg"] for row in comparison_rows])
    ranking = {}
    for arm in ("tpie", "fep"):
        predicted_ranks = rankdata([row[arm] for row in comparison_rows])
        ranking[arm] = {"squared_rank_displacement": float(np.sum((ordered_truth-predicted_ranks)**2)),
            "correct_pairs_of_45": sum((first["experimental_dg"]-second["experimental_dg"]) * (first[arm]-second[arm]) > 0 for first,second in combinations(comparison_rows,2)),
            "rows": [{"ligand": row["ligand"], "experimental_rank": float(ordered_truth[index]),
                      "predicted_rank": float(predicted_ranks[index])} for index,row in enumerate(comparison_rows)]}
    evidence = {"generated_at_utc": datetime.now(timezone.utc).isoformat(), "n": 10,
        "production_changed": False, "target_spearman": 0.6,
        "goal_achieved_by_new_method": goal_achieved,
        "completed_controls": {"input_pose": raw, "gradient": differences, "amide_geometry": amides},
        "mm_internal_geometry_diagnostic": mm,
        "geometry_status": geometry, "geometry_complete": complete, "metrics": available_metrics,
        "optimization_last_steps": optimization_progress,
        "archived_fep_rows": fep_rows,
        "rank_diagnostics": ranking,
        "experimental_range_kcal_mol": max(truth)-min(truth),
        "original_tpie_range_kcal_mol": float(np.ptp([row["tpie_pcm"] for row in previous["rows"]])),
        "original_multistart_strain_range_kcal_mol": float(np.ptp([row["solution_min_strain"] for row in sampled["rows"]])),
        "sources_sha256": {str(path): digest(path) for path in sources},
        "limitations": ["QM bound arm changes the geometry potential, removes positional restraint, and uses a capped pocket instead of whole-protein MM; not a pure one-variable model comparison",
            "No entropy, water population free energy or converged QM/MM free energy correction computed",
            "Finite-difference checks cover two atoms/axes only", "Development set, not blind validation"]}
    write_json(args.output / "CDK8_FOLLOWUP_ANALYSIS.json", evidence)
    lines = ["# CDK8：文献、FEP和结构对照分析", "", f"生成于 {evidence['generated_at_utc']}。", "",
        "## 核心结论", "", "开发集已有达到Spearman约0.6的候选，但仍需跨集验证；生产主流程未变更。" if goal_achieved else "尚未找到Spearman约0.6的新方法，生产主流程未变更。",
        "发现并修复计算实现问题，不等于发现了足以解释全部错排的唯一物理原因。", "",
        "本报告只把已完成且全10例覆盖的对照列为总体指标；QM几何队列状态见末尾，未完成部分不冒充完成。", "",
        "## 完整10例已完成指标", "", "|方法|Pearson|Spearman|Kendall|", "|---|---:|---:|---:|"]
    lines.extend(f"|{arm}|{value['pearson']:.3f}|{value['spearman']:.3f}|{value['kendall']:.3f}|" for arm,value in available_metrics.items())
    if "fep_archived_dg" not in available_metrics:
        lines.append("|既有FEP+21-4对照，非本轮新跑|0.663|0.648|0.378|")
    lines += ["", "## 可以确定的事情", "",
        "1. 统一PCM实际项后native与独立三态闭合，但三态排序仍弱：实现一致性不能代替物理预测有效性。",
        "2. 原输入姿态绕过k10后，全10例tPIE Spearman为−0.055、三态为−0.091。不能因为17号k10酰胺扭曲，就整体取消k10。",
        "3. 17号原输入与固定受体存在约2.11 Å的C–C接触，ASP173/ TYR32片段产生大排斥；k10减少碰撞，但酰胺平面偏离约55°。两类风险并存。",
        "4. 原输入17号酰胺接近平面（约1.81°偏离），游离QM后约4.54°；支持异常应变对制备/几何敏感，不证明实验真的付出25 kcal/mol应变。",
        "5. 13号与15号在原姿态、k10及目前已完成的QM姿态中均相互作用反向，不能仅靠17号这一异常解释整个CDK8结果。",
        "6. 前5例QM结合态优化降低同模型配体应变；总体指标仅在完整10例严格审计完成后发布，不用部分队列冒充。",
        "7. 13与17的碳酰氧最大梯度方向，0.001/0.0005 Å有限差分与解析梯度差最大约1.81×10⁻⁵ Ha/bohr；局部核查支持梯度实现，非全维/驻点证明。",
        "8. 5CEI与输入受体345个CA对齐RMSD约0.0043 Å，22号21个重原子一对一匹配RMSD约0.00051 Å，排除22号原始坐标框架明显错位；不证明其它9例都是独立晶体姿态。", "",
        "## 代表性分子对与归因边界", "",
        "13比15实验强约1.477 kcal/mol；原tPIE却弱5.799，已完成QM姿态后弱约11.396。归档FEP预测13强约1.220，方向正确。",
        "QM姿态15相对13的有利差中，ASP173分片贡献约−7.527 kcal/mol，TYR32约−1.289、LYS52约−0.943；仅靠一个碳酰氢键距离不能解释整个错排。分片是电子结构划分，不等于可独立删改的因果残基。",
        "碳酰氧到LYS52 NZ距离：13约3.210→3.223 Å，15约3.376→3.129 Å；说明QM优化改变局部接触，但短距离不自动证明该氢键解释排序。",
        "FEP也并非每例正确：17号预测约−12.15，实验约−10.40，残差约−1.75 kcal/mol，并且被排到过强。",
        "独立OpenFF2.2.0配体内部单点显示17号输入→k10总内部能约升1.91 kcal/mol，其中全分子扭转项约升10.61；说明MM本身也付出了内部代价，不能直接指控它完全没有酰胺平面势垒。",
        "17号完成QM结合态优化后应变约14.38，仍未消失；酰胺N三键夹角和由输入约360°变为k10约341.5°、QM结合态约334.1°。偏离平面既有扭转也有N锥形化，不能将全部现象都说成酰胺二面角旋转。",
        "这项是同分子/同电荷的真空MM几何诊断，不能与PCM应变混为同一个物理量，更不能直接加到tPIE当溶液自由能。",
        "21号输入存在拓扑距离6的H–H约1.13 Å接触，配体内部MM非键能很大，k10后大幅下降；重原子审计不足以覆盖这种氢几何问题。",
        "实验全队列跨度仅约1.57 kcal/mol，小体系/采样/制备误差足以改变许多分子对。没有实验重复误差，不能断言全部错排都是实验噪声。", "",
        f"原tPIE正确分子对{ranking['tpie']['correct_pairs_of_45']}/45，归档FEP为{ranking['fep']['correct_pairs_of_45']}/45；排名平方偏差分别{ranking['tpie']['squared_rank_displacement']:.0f}和{ranking['fep']['squared_rank_displacement']:.0f}。",
        "在10例无并列时，Spearman=1−Σd²/165；达到≥0.6要求Σd²≤66，而当前tPIE为126。必须修正一批具体错排，而不是只改变分数尺度。",
        "严格单调软压缩不会改变Spearman/Kendall，因此不能靠软压缩实现这个目标。", "",
        "## 为什么FEP更好不能简单归因于QM不够准", "",
        "FEP对配体化学变换计算结合与溶液两条采样腿之差；当前tPIE只读一套几何的跨片段相互作用。",
        "公开CDK8运行脚本有muVT水设置与Force Field Builder；本项目没有复现这些采样成本，不能声称同等制备同等采样的公平能量层比较。",
        "自由能最佳实践全文强调结合模式、水交换、实验pH、相关时间及相空间重叠；QM/MM文献指出不能用相互作用差直接替代总势能差做自由能修正。",
        "因此DFTB、局部刚性受体、形变、局部水、未计统计贡献均有可能，现有证据不足以按原因分配百分比。", "",
        "## 下一步与默认流程的边界", "",
        "同模型QM几何10例与原始姿态对照已经完成；需要继续按预先规则区分局部受体刚性、化学态、水和关键相互作用模型误差。" if complete else "先完成同模型QM几何10例与原始姿态对照；再按预先规则检查短接触、酰胺平面性、力场扭转/电荷和关键氢键。",
        "QM几何臂同时改变势能模型、位置约束和受体范围，不能把它称为只改变QM方法的单变量消融。",
        "有结构风险可增加审计标记，但自动改构象/松动受体的默认修复须先有匹配对照与Galectin3等跨集回归，不按亲和力挑姿态。",
        "水/质子化及匹配结合态–游离态采样仍未完成；不把最低点当熵，不把快照均值当自由能。", "",
        "## QM几何队列状态", "", "|配体|状态|", "|---|---|"]
    lines.extend(f"|{name}|{status}|" for name,status in geometry.items())
    lines += ["", "详细文献与公开FEP命令：CDK8_LITERATURE_FEP_CONTROLS.zh-CN.md。",
        "原始失败、PIEDA补算、严格SCF/坐标审计及来源哈希均保留在运行目录。"]
    (args.output / "CDK8_FOLLOWUP_ANALYSIS.zh-CN.md").write_text("\n".join(lines)+"\n")


if __name__ == "__main__":
    main()
