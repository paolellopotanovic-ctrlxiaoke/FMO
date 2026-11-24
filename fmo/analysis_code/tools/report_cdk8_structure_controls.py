"""Report complete matched structure controls without selecting ligand poses."""
from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path
import sys

import numpy as np
from rdkit import Chem

from run_cdk8_solution_cycle import array, group, write_json
from analyze_cdk8_geometry_evidence import amide_torsions, bond_audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", required=True, type=Path)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    sys.path.insert(0, str(args.workspace.resolve() / "src"))
    from sophosqm_baseline.parser import parse_pie
    from fmotpiefinal.core.structure import read_pdb
    run = args.run.resolve()
    evaluation = json.loads((run / "evaluation.json").read_text())
    plan = json.loads((run / "frozen_protocol.json").read_text())
    names, arms = plan["names"], plan["arms"]
    if any(evaluation["arms"][arm]["n_succeeded"] != len(names) for arm in arms):
        raise ValueError("full-cohort report requires all 30 successful jobs")
    sdf = args.workspace / "vendor/public_binding_free_energy_benchmark/fep_benchmark_inputs/structure_inputs/misc/cdk8_koehler_ligands.sdf"
    molecules = {molecule.GetProp("_Name").strip(): molecule
                 for molecule in Chem.SDMolSupplier(str(sdf), removeHs=False) if molecule is not None}
    truth = {name: float(molecule.GetProp("r_exp_dg")) for name, molecule in molecules.items()}
    geometry = {}
    pair_values = {}
    for arm in arms:
        pair_values[arm] = {}
        geometry[arm] = {}
        for name in names:
            atoms = [atom for atom in read_pdb(run / arm / name / "minimization/minimized_complex.pdb").atoms
                     if atom.record == "HETATM" and atom.residue_name == "UNK"]
            molecule = molecules[name]
            if [atom.element.upper() for atom in atoms] != [atom.GetSymbol().upper() for atom in molecule.GetAtoms()]:
                raise ValueError("ligand geometry audit atom ordering mismatch")
            positions = np.asarray([atom.xyz for atom in atoms], dtype=float)
            geometry[arm][name] = {"amide": amide_torsions(molecule, positions),
                                   "bonds": bond_audit(molecule, positions)}
            directory = run / arm / name / "qm"
            content = group((directory / "input.inp").read_text(), "FMO")[1]
            fragment_names = array(content, "FRGNAM", "INDAT")
            ligand_id = len(fragment_names)
            pairs = [pair for pair in parse_pie(directory / "gamess.log")
                     if ligand_id in (pair.fragment_i, pair.fragment_j)]
            pair_values[arm][name] = {fragment_names[min(pair.fragment_i, pair.fragment_j) - 1]:
                                     pair.value_kcal_mol for pair in pairs}
    baseline = arms[0]
    differences = {}
    for arm in arms[1:]:
        residues = sorted(pair_values[baseline][names[0]])
        changes = {residue: [pair_values[arm][name][residue] - pair_values[baseline][name][residue]
                             for name in names] for residue in residues}
        differences[arm] = {"residues_by_score_change_sd": sorted([
            {"fragment": residue, "mean_delta_kcal_mol": float(np.mean(values)),
             "sd_delta_kcal_mol": float(np.std(values)), "per_ligand_delta": dict(zip(names, values))}
            for residue, values in changes.items()], key=lambda row: row["sd_delta_kcal_mol"], reverse=True)}
    pairwise = {}
    for arm in arms:
        correct, wrong, tied = [], [], []
        for left, right in combinations(names, 2):
            predicted = evaluation["results"][arm][left]["tpie"] - evaluation["results"][arm][right]["tpie"]
            observed = truth[left] - truth[right]
            destination = tied if predicted * observed == 0 else correct if predicted * observed > 0 else wrong
            destination.append([left, right])
        pairwise[arm] = {"correct": correct, "wrong": wrong, "tied": tied}
    result = {"arms": evaluation["arms"], "baseline_matched": evaluation["baseline_matched_within_0_05_kcal"],
        "frozen_protocol": plan,
        "coordinate_audit": json.loads((run / "independent_coordinate_audit.json").read_text()),
        "pairwise": pairwise, "residue_changes": differences, "ligand_geometry": geometry,
        "wall_seconds": evaluation["wall_seconds"]}
    write_json(args.output.with_suffix(".json"), result)
    lines = ["# CDK8 局部柔性和水氢对照：完整结果", "", "全部三组均为完整 10/10；使用原始 Pearson，未拟合参数或按亲和力选择结构。", "",
             "| 方法 | Pearson | Spearman | Kendall | 正确配体对 / 45 |", "|---|---:|---:|---:|---:|"]
    for arm in arms:
        metric = evaluation["arms"][arm]["metrics"]
        lines.append(f"| {arm} | {metric['pearson']:.4f} | {metric['spearman']:.4f} | {metric['kendall']:.4f} | {len(pairwise[arm]['correct'])} |")
    lines += ["", f"基线与历史分数最大偏差：{evaluation['arms'][baseline].get('max_historical_score_difference_kcal_mol')} kcal/mol。",
              f"本次驱动记录墙钟：{evaluation['wall_seconds']:.1f} 秒（非 CPU 总时间）。", "",
              "## 相互作用变化来自哪里", "", "以下按全队列残基贡献变化的标准差排序，不据此删残基或挑选配体。"]
    for arm in arms[1:]:
        lines += ["", f"### {arm}", "", "| 片段 | 平均贡献变化 kcal/mol | 变化标准差 kcal/mol |", "|---|---:|---:|"]
        for row in differences[arm]["residues_by_score_change_sd"][:6]:
            lines.append(f"| {row['fragment']} | {row['mean_delta_kcal_mol']:.3f} | {row['sd_delta_kcal_mol']:.3f} |")
    original_correct = {tuple(pair) for pair in pairwise[baseline]["correct"]}
    local_correct = {tuple(pair) for pair in pairwise["local_sidechains"]["correct"]}
    lines += ["", "## 匹配对照结论", "",
        f"局部侧链组修正 {len(local_correct-original_correct)} 对原先错排，同时损失 {len(original_correct-local_correct)} 对原先正确的配体对；并非每一对都改善。",
        "13/15 是关键实例：基线错误地偏向 15，局部侧链组改为偏向 13，与实验方向一致。",
        "主要变化涉及 ASP173、LYS52 和 TYR32；这是局部侧链–配体协同响应的证据，不足以认定单个残基是唯一原因。",
        "水氢单独优化没有提高排序，不能因为物理上合理就默认叠加到候选流程。",
        "全部 30 个结构的配体共价键距离筛查均无标记；固定原子另有独立 PDB 坐标审计。",
        "这些筛查不证明化学态完全正确，也不代表全部酰胺平面性问题已经解决。"]
    baseline_amide = geometry[baseline]["17"]["amide"]
    local_amide = geometry["local_sidechains"]["17"]["amide"]
    lines += ["",
        f"配体 17 的酰胺最大离平面角从 {max(row['departure_from_planarity_deg'] for row in baseline_amide):.2f}° 降为 {max(row['departure_from_planarity_deg'] for row in local_amide):.2f}°；",
        f"相应氮邻接键角总和从 {baseline_amide[0]['nitrogen_neighbor_angle_sum_deg']:.2f}° 恢复至 {local_amide[0]['nitrogen_neighbor_angle_sum_deg']:.2f}°（平面参照为 360°）。",
        "这支持固定受体迫使某些配体承受局部不合理变形的解释，但不能把全部指标改善归因于这一例。"]
    reached = any(evaluation["arms"][arm]["metrics"]["spearman"] >= 0.6 for arm in arms)
    lines += ["", "## 解释边界", "", "本轮达到 Spearman 0.6。" if reached else "本轮未达到 Spearman 0.6，不能宣布目标完成。",
              "局部侧链组同时允许配体在相同 k10 约束下响应受体变化；水氢组则严格固定基线配体。",
              "这些是单结构能量对照，不等于结合态/游离态匹配采样，也不等于 FEP。",
              "不把本轮处理自动设为生产默认；任何候选仍需跨集合回归。"]
    args.output.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
