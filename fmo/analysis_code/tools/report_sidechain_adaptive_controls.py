"""Report frozen sidechain controls and case-level structural/ranking evidence."""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np
from rdkit import Chem
from scipy.stats import rankdata

from fmotpiefinal.core.bound_geometry import load_bound_ligand
from report_local_sidechain_transfer import load_truth
from run_cdk8_solution_cycle import digest, metrics, write_json


def pair_errors(observed, predicted):
    errors = np.zeros(len(observed), dtype=int)
    comparable = np.zeros(len(observed), dtype=int)
    for first in range(len(observed)):
        for second in range(first + 1, len(observed)):
            actual = observed[first] - observed[second]
            if actual == 0:
                continue
            comparable[first] += 1
            comparable[second] += 1
            if actual * (predicted[first] - predicted[second]) <= 0:
                errors[first] += 1
                errors[second] += 1
    return errors.tolist(), comparable.tolist()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", required=True, type=Path)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    root = args.run.resolve()
    plan = json.loads((root / "frozen_protocol.json").read_text())
    result = {"schema": "adaptive-sidechain-diagnostic-v1", "groups": {},
        "protocol_sha256": digest(root / "frozen_protocol.json"),
        "interpretation": "Post-hoc associations, not causal attribution or rule-selection inputs",
        "production_default_changed": False}
    progress = {arm: json.loads((root / arm / "progress.json").read_text())
        if (root / arm / "progress.json").exists() else {} for arm in plan["arms"]}
    lines = ["# 扩展侧链图：冻结对照与逐例诊断", "",
        "A：严格种子＋两层邻接；B：再加入原始碰撞与配体位移共同证据。",
        "以下指标仅在完整集合通过时报告；未触发集合是验收基线复用，不算新 QM。", "",
        "| 集合 | n | 原始碰撞/位移关联配体数 | 基线 P / S / K | A P / S / K | B P / S / K |",
        "|---|---:|---:|---|---|---|"]
    display = lambda values: " / ".join(f"{values[key]:.3f}" for key in ("pearson", "spearman", "kendall"))
    for target in plan["targets"]:
        cohort = plan["cohorts"][plan["arms"][0]][target]
        names = cohort["names"]
        truth = load_truth(args.workspace.resolve(), target)
        if set(names) != set(truth):
            raise ValueError(f"full label identity mismatch: {target}")
        observed = [truth[name] for name in names]
        baseline = [-cohort["baseline"][name]["tpie"] for name in names]
        events = json.loads((root / target / "motion_seed_graph.json").read_text())["motion_events"]
        counts = Counter(event["ligand"] for event in events)
        molecules = {molecule.GetProp("_Name").strip(): molecule for molecule in
            Chem.SDMolSupplier(cohort["sdf"], removeHs=False) if molecule is not None}
        if set(molecules) != set(names):
            raise ValueError(f"full molecular identity mismatch: {target}")
        errors, comparable = pair_errors(observed, baseline)
        group = {"n": len(names), "motion_event_count": len(events),
            "motion_ligand_count": len(counts), "baseline": metrics(observed, baseline),
            "arms": {}, "cases": []}
        ranks = rankdata(observed)
        baseline_ranks = rankdata(baseline)
        for index, name in enumerate(names):
            molecule = molecules[name]
            original = np.asarray(molecule.GetConformer().GetPositions())
            bound = np.asarray([atom.xyz for atom in load_bound_ligand(
                cohort["sdf"], name, cohort["complexes"][name]).atoms])
            heavy = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1]
            displacement = np.linalg.norm(bound[heavy] - original[heavy], axis=1)
            group["cases"].append({"ligand": name, "experimental_minus_dG": truth[name],
                "baseline_minus_tpie": baseline[index], "baseline_rank_error": float(baseline_ranks[index]-ranks[index]),
                "baseline_wrong_or_tied_pairs": errors[index], "comparable_pairs": comparable[index],
                "n_motion_events": counts[name], "max_heavy_displacement_A": float(max(displacement)),
                "rms_heavy_displacement_A": float(np.sqrt(np.mean(displacement**2))),
                "motion_evidence": [event for event in events if event["ligand"] == name]})
        columns = []
        for arm in plan["arms"]:
            rows = progress[arm].get(target, {})
            failures = {name: rows.get(name, {"status": "MISSING"}) for name in names
                if rows.get(name, {}).get("status") != "SUCCEEDED"}
            detail = {"n_expected": len(names), "n_succeeded": len(names)-len(failures),
                "failures": failures, "selected_residues": plan["cohorts"][arm][target]["local_residues"],
                "n_new_qm": sum(bool(row.get("new_qm_executed")) for row in rows.values())}
            if failures:
                columns.append(f"未完整通过 {detail['n_succeeded']}/{len(names)}")
            else:
                candidate = [-rows[name]["tpie"] for name in names]
                detail["metrics"] = metrics(observed, candidate)
                candidate_errors, _ = pair_errors(observed, candidate)
                candidate_ranks = rankdata(candidate)
                for index, case in enumerate(group["cases"]):
                    case[arm] = {"minus_tpie": candidate[index],
                        "rank_error": float(candidate_ranks[index]-ranks[index]),
                        "wrong_or_tied_pairs": candidate_errors[index],
                        "execution_path": rows[names[index]]["execution_path"]}
                columns.append(display(detail["metrics"]))
            group["arms"][arm] = detail
        lines.append(f"| {target} | {len(names)} | {len(counts)} | {display(group['baseline'])} | {' | '.join(columns)} |")
        result["groups"][target] = group
    lines += ["", "## 如何读逐例证据", "",
        "配套 JSON 同时保存每例原始碰撞、位移、基线错排对数、各对照错排对数和实际执行路径。",
        "错排对数只比较实验值不同的配体；预测并列按未正确排序计，不等同 Kendall 的定义。",
        "位移是在同一坐标系原位比较，不做配体自身拟合；不是自由能、应变能或因果证明。",
        "B 相对 A 同时增加种子及其邻接残基，因此是种子检测策略的整体对照，不能归因到单个残基。",
        "全部集合使用相同规则；不因哪个对照在某个集合更高而逐集合挑选。",
        "尚未触发不代表结构正确，也不代表溶剂、质子化、采样等遗漏已经解决。",
        "软压缩的 Pearson 及独立固定坐标审计见各对照的独立验收报告。"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines)+"\n")
    write_json(args.output.with_suffix(".json"), result)
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
