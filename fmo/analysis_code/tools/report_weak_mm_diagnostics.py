"""Describe matched ligand geometry costs without adding them to tPIE."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from rdkit import Chem

from fmotpiefinal.core.bound_geometry import load_bound_ligand
from run_cdk8_solution_cycle import digest, write_json
from run_local_sidechain_transfer import directory_key


def angle_changes(molecule, before, after):
    changes = []
    for center in molecule.GetAtoms():
        neighbors = list(center.GetNeighbors())
        for index, first in enumerate(neighbors):
            for last in neighbors[index+1:]:
                indices = [first.GetIdx(), center.GetIdx(), last.GetIdx()]
                angles = []
                for positions in (before, after):
                    first_vector = positions[indices[0]]-positions[indices[1]]
                    second_vector = positions[indices[2]]-positions[indices[1]]
                    norm = np.linalg.norm(first_vector)*np.linalg.norm(second_vector)
                    if norm == 0:
                        raise ValueError("overlapping bonded atoms")
                    angles.append(float(np.degrees(np.arccos(np.clip(np.dot(first_vector, second_vector)/norm, -1, 1)))))
                changes.append({"atoms_zero_based": indices,
                    "elements": [atom.GetSymbol() for atom in (first, center, last)],
                    "baseline_angle_deg": angles[0], "candidate_angle_deg": angles[1],
                    "absolute_change_deg": abs(angles[1]-angles[0])})
    return sorted(changes, key=lambda row: row["absolute_change_deg"], reverse=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    status = json.loads((args.run / "final_status.json").read_text())
    if status["status"] != "COMPLETE":
        raise ValueError("complete diagnostic inventory required")
    rows = json.loads((args.run / "progress.json").read_text())
    plan = json.loads((args.source / "frozen_protocol.json").read_text())
    results = {}
    lines = ["# 弱集合配体内部几何代价：匹配 MM 诊断", "",
        "OpenFF 2.2.0 / NAGL；每个配体相同原子、参数和电荷，只换坐标。能量单位 kcal/mol。",
        "这是无受体、无溶剂的诊断，不是结合自由能；绝不把不同分子的绝对总能量直接作排序。", "",
        "| 集合 | n | B−基线内部能差 中位数 | 最小 / 最大差 | 局部最小化未达力阈值例数 |",
        "|---|---:|---:|---|---:|"]
    for target, group in rows.items():
        cohort = plan["cohorts"][target]
        molecules = {molecule.GetProp("_Name").strip(): molecule for molecule in Chem.SDMolSupplier(cohort["sdf"], removeHs=False) if molecule is not None}
        if set(group) != set(cohort["names"]):
            raise ValueError("cohort identity mismatch")
        cases = {}
        for name in cohort["names"]:
            row = group[name]
            if row["status"] != "SUCCEEDED":
                raise ValueError("incomplete MM result")
            before = row["energies"]["baseline_k10"]
            after = row["energies"]["sidechain_B"]
            delta = after["total_kcal_mol"]-before["total_kcal_mol"]
            source_path = Path(cohort["complexes"][name])
            candidate_path = args.source / target / directory_key(name) / "minimization/minimized_complex.pdb"
            if not cohort["local_residues"]:
                candidate_path = source_path
            before_positions = np.asarray([atom.xyz for atom in load_bound_ligand(cohort["sdf"], name, source_path).atoms])
            after_positions = np.asarray([atom.xyz for atom in load_bound_ligand(cohort["sdf"], name, candidate_path).atoms])
            changes = angle_changes(molecules[name], before_positions, after_positions)
            cases[name] = {"candidate_minus_baseline_internal_mm_kcal_mol": delta,
                "force_component_differences_kcal_mol": {key: after["components_kcal_mol"][key]-value for key, value in before["components_kcal_mol"].items()},
                "largest_angle_changes": changes[:5],
                "local_reference_force_tolerance_passed": row["isolated_local_minimum_force_tolerance_passed"]}
        deltas = [row["candidate_minus_baseline_internal_mm_kcal_mol"] for row in cases.values()]
        unconverged = [name for name, row in cases.items() if not row["local_reference_force_tolerance_passed"]]
        results[target] = {"n": len(cases), "cases": cases, "local_reference_not_converged": unconverged,
            "median_internal_mm_difference_kcal_mol": float(np.median(deltas)),
            "range_internal_mm_difference_kcal_mol": [min(deltas), max(deltas)]}
        lines.append(f"| {target} | {len(cases)} | {np.median(deltas):.3f} | {min(deltas):.3f} / {max(deltas):.3f} | {len(unconverged)} |")
    for target, group in results.items():
        lines += ["", f"## {target}：最明显的内部能上升例", ""]
        name, row = max(group["cases"].items(), key=lambda item: item[1]["candidate_minus_baseline_internal_mm_kcal_mol"])
        lines += [f"{name}：B−基线 = {row['candidate_minus_baseline_internal_mm_kcal_mol']:.3f} kcal/mol。",
            "分项差：" + str({key: round(value, 3) for key, value in row["force_component_differences_kcal_mol"].items()}) + "。"]
        if row["largest_angle_changes"]:
            angle = row["largest_angle_changes"][0]
            lines += [f"最大角度变化：原子索引 {angle['atoms_zero_based']}（{angle['elements']}），{angle['baseline_angle_deg']:.2f}° → {angle['candidate_angle_deg']:.2f}°。"]
        if group["local_reference_not_converged"]:
            lines += ["孤立配体参照未达力阈值：" + ", ".join(group["local_reference_not_converged"]) + "；不能称其为已收敛局部最低点。"]
    lines += ["", "## 解读", "",
        "BACE CAT-13i 出现明显内部角度代价，说明仅审计重原子碰撞和酰胺平面仍可能漏掉含氢键角异常。",
        "该异常不能独自解释整个 BACE 集合错排；需要完整姿态保持与氢约束对照，不删除异常配体。",
        "CHK1 内部代价下降也没有同步带来排序提升：局部模型几何改善不是亲和力排序的充分条件。",
        "六例孤立参照未达力阈值时仅保留近似参照状态；本表直接比较两个绑定几何，差值不依赖该参照。",
        "不按内部 MM 最低点选择 tPIE，更不利用实验标签筛分支。"]
    args.output.write_text("\n".join(lines)+"\n")
    write_json(args.output.with_suffix(".json"), {"groups": results, "status": status,
        "diagnostic_progress_sha256": digest(args.run / "progress.json"),
        "prediction_score_changed": False, "labels_used": False})
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
