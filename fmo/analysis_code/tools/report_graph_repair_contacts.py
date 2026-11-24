"""Audit residue contributions and remaining pair errors without changing tPIE."""
from __future__ import annotations

import argparse
from collections import defaultdict
from itertools import combinations
import json
from pathlib import Path

import numpy as np

from fmotpiefinal.core.parser import parse_fragment_names, parse_pie
from report_local_sidechain_transfer import load_truth
from run_cdk8_solution_cycle import digest, write_json


def contact_map(records, names, ligand_fragment):
    contacts = {}
    for record in records:
        if ligand_fragment not in (record.fragment_i, record.fragment_j):
            continue
        other = record.fragment_j if record.fragment_i == ligand_fragment else record.fragment_i
        if other == ligand_fragment or other not in names or names[other] in contacts:
            raise ValueError("nonunique or unidentified ligand contact")
        contacts[names[other]] = record.value_kcal_mol
    if len(contacts) != ligand_fragment - 1:
        raise ValueError("incomplete ligand pair inventory")
    return contacts


def ordering_errors(names, truth, before, after):
    rows = []
    for first, second in combinations(names, 2):
        gap = truth[first] - truth[second]
        if abs(gap) < 1.e-10:
            continue
        old = gap * (before[second] - before[first])
        new = gap * (after[second] - after[first])
        category = ("corrected" if old < 0 < new else "regressed" if new < 0 < old
                    else "still_wrong" if old < 0 and new < 0 else "other")
        rows.append({"first": first, "second": second, "experimental_gap_kcal_mol": gap,
                     "baseline_energy_gap": before[first] - before[second],
                     "candidate_energy_gap": after[first] - after[second], "category": category})
    return sorted(rows, key=lambda row: abs(row["experimental_gap_kcal_mol"]), reverse=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", required=True, type=Path)
    parser.add_argument("--registry", required=True, type=Path)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--verified", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    frozen = json.loads((args.run / "frozen_protocol.json").read_text())
    verified = json.loads(args.verified.read_text())
    if verified["status"] != "COMPLETE_INDEPENDENTLY_VERIFIED" or digest(args.registry) != frozen["registry_sha256"]:
        raise ValueError("unverified repair or changed registry")
    cohort = json.loads(args.registry.read_text())["cohorts"][frozen["target"]]
    names = sorted(frozen["names"])
    truth = load_truth(args.workspace.resolve(), frozen["target"])
    if set(names) != set(cohort["names"]) or set(names) != set(truth) or set(names) != set(verified["evidence"]):
        raise ValueError("matched complete cohort required")
    per_case, residue_deltas = {}, defaultdict(list)
    for name in names:
        task = frozen["tasks"][name]
        original = Path(cohort["sources"][name])
        candidate = Path(task["input"])
        old_log, new_log = original.parent / "run/gamess.log", candidate.parent / "run/gamess.log"
        evidence = verified["evidence"][name]
        if (digest(original) != task["original_sha256"] or digest(old_log) != cohort["baseline"][name]["log_sha256"]
                or digest(candidate) != evidence["input_sha256"] or digest(new_log) != evidence["log_sha256"]):
            raise ValueError("contact evidence changed")
        ligand = task["runtime"]["ligand_fragment"]
        old = contact_map(parse_pie(old_log), parse_fragment_names(old_log), ligand)
        new = contact_map(parse_pie(new_log), parse_fragment_names(new_log), ligand)
        if old.keys() != new.keys():
            raise ValueError("contact residue inventory changed")
        before, after = sum(old.values()), sum(new.values())
        if abs(before - cohort["baseline"][name]["tpie"]) > .002 or abs(after - evidence["tpie"]) > .002:
            raise ValueError("residue contribution sum does not close to verified tPIE")
        deltas = {residue: new[residue] - old[residue] for residue in old}
        for residue, value in deltas.items():
            residue_deltas[residue].append(value)
        per_case[name] = {"before": before, "after": after, "delta_kcal_mol": after - before,
                          "residue_delta_kcal_mol": deltas,
                          "original_log_sha256": digest(old_log), "candidate_log_sha256": digest(new_log)}
    residues = [{"residue": name, "n": len(values), "mean_delta": float(np.mean(values)),
                 "median_delta": float(np.median(values)), "median_abs_delta": float(np.median(np.abs(values))),
                 "min_delta": min(values), "max_delta": max(values)} for name, values in residue_deltas.items()]
    residues.sort(key=lambda row: row["median_abs_delta"], reverse=True)
    before = {name: row["before"] for name, row in per_case.items()}
    after = {name: row["after"] for name, row in per_case.items()}
    errors = ordering_errors(names, truth, before, after)
    report = {"status": "COMPLETE_VERIFIED", "n": len(names), "per_case": per_case,
              "residues_by_median_abs_change": residues, "pair_ordering": errors,
              "counts": {category: sum(row["category"] == category for row in errors)
                         for category in ("corrected", "regressed", "still_wrong", "other")},
              "registry_sha256": digest(args.registry), "verified_sha256": digest(args.verified),
              "driver_sha256": digest(__file__), "production_default_changed": False,
              "warning": "Post-hoc attribution, not a new score or proof of causality. Pair errors are dependent; experimental gaps are not independent significance tests."}
    write_json(args.output.with_suffix(".json"), report)
    lines = ["# BACE图修复：残基贡献与尚未解决的错排", "", "36例保持原生产模型，只改变被图标记的共轭NH2氢；以下仍是纯tPIE，不选择ES等分量重新打分。",
             "每例前后PIEDA残基加和均与已独立验收tPIE闭合；日志、输入及参照指纹重新核验。", "",
             "## 变化最大的残基", "", "正的Δ表示修复后该残基相互作用较不利；不是该残基被删除。按跨配体中位绝对变化排序。", "",
             "| 残基 | 例数 | Δ中位 | Δ范围，kcal/mol |", "|---|---:|---:|---:|"]
    lines.extend(f"| {row['residue']} | {row['n']} | {row['median_delta']:.3f} | {row['min_delta']:.3f}～{row['max_delta']:.3f} |" for row in residues[:12])
    lines.extend(["", "## 排序还错在哪里", "", f"修正{report['counts']['corrected']}对、变差{report['counts']['regressed']}对、仍错{report['counts']['still_wrong']}对。实验完全并列对不纳入此表。",
                  "下表仅用于诊断，不能剔除这些样本或根据实验标签改变修复规则。真值差是−ΔG之差；能量差为tPIE之差。两者同号意味着排序相反。", "",
                  "| 配体对 | 实验差，kcal/mol | 原能量差 | 修复后能量差 |", "|---|---:|---:|---:|"])
    lines.extend(f"| {row['first']} / {row['second']} | {row['experimental_gap_kcal_mol']:.3f} | {row['baseline_energy_gap']:.3f} | {row['candidate_energy_gap']:.3f} |"
                 for row in [row for row in errors if row["category"] == "still_wrong"][:12])
    lines.extend(["", "## 解释边界", "", "局部氢修复会通过整个FMO/PCM计算影响不止一个残基，不能把所有变化称为两个H的孤立键能。",
                  "残基贡献仅解释分数如何变化，不能证明缺失的真实自由能来自某个特定残基。",
                  "仍需在结合环境下检查k10的势能项和溶液/结合态差别；不能把这份事后诊断变成按实验标签选流程。", ""])
    args.output.with_suffix(".md").write_text("\n".join(lines))
    print(json.dumps({"status": report["status"], "counts": report["counts"], "top_residues": residues[:4]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
