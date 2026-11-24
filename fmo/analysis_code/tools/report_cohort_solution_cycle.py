"""Verify full cohort coverage and compare frozen PCM cycle readouts."""
from __future__ import annotations

import argparse
import io
import json
from pathlib import Path
import sys

import numpy as np

from fmotpiefinal.core.fmo_convergence import final_fmo_convergence
from report_bace_extensions import paired_rank_summary
from report_local_sidechain_transfer import load_truth
from run_cdk8_solution_cycle import digest, write_json
from run_cohort_solution_cycle import evaluate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repaired-receptor", action="store_true")
    parser.add_argument("--repaired-minima", action="store_true")
    args = parser.parse_args()
    root, workspace = args.run.resolve(), args.workspace.resolve()
    sys.path.insert(0, str(workspace / "src"))
    from sophosqm_baseline.relative_total_energy import parse_converged_ligand_optimization, parse_full_energy
    evaluation_root = root / "receptor_data_title_repair" if args.repaired_receptor else root
    if args.repaired_minima:
        evaluation_root = root / "ligand_minimum_rfo_recovery"
    protocol_path = evaluation_root / "recovered_protocol.json" if (args.repaired_receptor or args.repaired_minima) else root / "frozen_protocol.json"
    frozen = json.loads(protocol_path.read_text())
    prior = json.loads((evaluation_root / "evaluation.json").read_text())
    if set(prior["rows"]) != set(frozen["tasks"]):
        raise ValueError("not all frozen tasks accounted for")
    verified = {}
    for key, task in frozen["tasks"].items():
        row = prior["rows"][key]
        if row["status"] != "SUCCEEDED":
            raise ValueError("incomplete cycle; do not publish complete-cohort metrics")
        path = Path(task["input"])
        log = path.parent / "gamess.log"
        if digest(path) != row["input_sha256"] or digest(log) != row["log_sha256"]:
            raise ValueError("cycle input/log fingerprint mismatch")
        if digest(Path(task["runtime"]["gamess_dir"]) / "gamess.00.x") != row["executable_sha256"]:
            raise ValueError("QM executable changed")
        text = log.read_text()
        if task["state"] == "free_ligand_minimum":
            _, energy = parse_converged_ligand_optimization(text, ligand_deck=task["bound_ligand_deck"])
        else:
            energy = parse_full_energy(text, input_text=path.read_text())
        if task["state"] in {"complex", "receptor"}:
            final_fmo_convergence(path.read_text(), io.StringIO(text))
        if abs(energy-row["energy_hartree"]) > 1.e-9:
            raise ValueError("independent energy extraction mismatch")
        verified[key] = row
    report = evaluate(verified, frozen, workspace)
    truth = load_truth(workspace, frozen["target"])
    names = sorted(frozen["names"])
    observed = [truth[name] for name in names]
    baseline = [-frozen["baseline"][name] for name in names]
    report["paired_comparisons"] = {score: paired_rank_summary(observed, baseline,
        [-report["scores"][name][score] for name in names])
        for score in report["metrics"] if score != "baseline_tpie"}
    report["strain_vs_matched_tpie"] = paired_rank_summary(observed,
        [-report["scores"][name]["matched_tpie"] for name in names],
        [-report["scores"][name]["tpie_plus_ligand_strain"] for name in names])
    report["matched_tpie_max_deviation_from_selective_baseline_kcal_mol"] = max(
        abs(report["scores"][name]["matched_tpie"]-frozen["baseline"][name]) for name in names)
    report["ligand_strain_distribution"] = {
        "min": min(report["scores"][name]["ligand_strain"] for name in names),
        "median": float(np.median([report["scores"][name]["ligand_strain"] for name in names])),
        "max": max(report["scores"][name]["ligand_strain"] for name in names)}
    report["scope"] = "Full-cohort development diagnostics; single PCM local minima are not free-state ensembles or binding free energies"
    write_json(args.output.with_suffix(".json"), report)
    lines = [f"# {frozen['target']} 完整队列一致三态结果", "",
        f"{len(names)}例全部覆盖，{len(verified)}项独立QM验收。严格保留tPIE与三态分数的命名区别。", "",
        "| 预先声明的方法 | Pearson | Spearman | Kendall |",
        "|---|---:|---:|---:|"]
    for score, values in report["metrics"].items():
        lines.append(f"| {score} | {values['pearson']:.4f} | {values['spearman']:.4f} | {values['kendall']:.4f} |")
    lines += ["", "- 三态均为DFTB3/D3/PCM，统一IDISP=0，使用全FMO2总能量；绝不把选择性FMO的伪总能量用于相减。",
        "- 单片段配体FMO1不是只算静电：它仍包含完整配体片段内部的DFTB3量子计算；没有可额外枚举的跨片段二聚体。",
        "- 同几何三态不包含独立结构重组；游离配体优化只增加单个局部最小值的参照，不代表溶液人口、熵或自由能。",
        "- tPIE+strain、tPIE+cavity为预声明单位系数诊断组合，不是拟合，不据此自动改变主流程。",
        "- 原始共享受体标题空行错误的输入/失败日志完整保留；修复仅移除一空行，物理模型和几何不变。" if args.repaired_receptor or args.repaired_minima else "",
        "- 配体优化数值恢复保留原始失败记录，未降低收敛阈值；不根据实验或最终能量挑选结果。" if args.repaired_minima else "",
        "", f"配体应变范围与中位数（kcal/mol）：{report['ligand_strain_distribution']}。",
        f"匹配全AB与选择性基线tPIE最大偏差：{report['matched_tpie_max_deviation_from_selective_baseline_kcal_mol']:.6f} kcal/mol。",
        "这里既改变了全AB/选择性计算方式，又统一了PCM非极性开关，不能把差异全归为数值误差；应变作用同时与匹配全AB的tPIE比较。",
        "", "配对bootstrap和逐例完整证据见同名JSON；反复开发同一集合的选择偏差不由bootstrap消除。"]
    args.output.with_suffix(".md").write_text("\n".join(lines)+"\n")


if __name__ == "__main__":
    main()
