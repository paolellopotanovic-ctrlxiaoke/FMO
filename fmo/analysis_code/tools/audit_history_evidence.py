"""Build or verify a portable, hashed index of local historical engineering evidence."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


REPORTS = {
    "H01": ("docs/FMO_TOTAL_ENERGY_RESULTS_CN.md", ["## 排序对照", "## 解释和边界"]),
    "H02": ("docs/FMO_OFFICIAL_BINDING_AUDIT_CN.md", ["## 官方 GAMESS 原生数值交叉验证", "## 官方允许的自由 B 构象"]),
    "H03": ("docs/FMO87_ABLATION_RESULTS_2026-09-30_CN.md", ["**协议更正", "### E18", "### E23", "### E19", "### E24"]),
    "H04": ("docs/FMO_TPIE_ONLY_PROTOCOL_AUDIT_20261015_CN.md", ["| 配对比较", "## 重要更正"]),
    "H05": ("docs/FMO_TPIE_CLEAN_WATER_K10_VALIDATION_CN.md", ["## 统一结果", "## 原始87的 target"]),
    "H06": ("docs/FMO87_NEW_FULL_PIPELINE_TPIE_20261015_CN.md", ["**同一套执行步骤**", "旧 87"]),
    "H07": ("docs/FMO_PDBFIXER_FULL_PIPELINE_RESULTS_CN.md", ["## 结论", "## 结构审计"]),
    "H08": ("docs/FMO_UNIFIED_SOURCE_WATER_RUN_LEDGER_20261015_CN.md", ["## 冷启动运行记录"]),
    "H09": ("docs/FMO_FXA_REGISTERED_FAILURE_RESCUE_20261016_CN.md", ["#", "**最终冷启动已通过**"]),
    "H10": ("docs/FMO_OPENBIND_SOFTCAP_ABLATION_20261016_CN.md", ["## 先分清两种", "## 只用本轮"]),
    "H11": ("docs/FMO_PERCENTILE_TRANSFER_20261013_CN.md", ["## 全量结果", "## 重要解释"]),
    "H12": ("docs/FMO_UNIFIED_MW_AND_OPENBIND_ABLATION_20261016_CN.md", ["## 当前默认臂", "## 统一固定公式"]),
    "H13": ("docs/FMO_COMPLETE_RESULTS_REPORT_CN.md", ["## 4.", "## 5."]),
    "H14": ("docs/OPENBIND_FAILURE_ROOT_CAUSE_AND_AUTO_PIPELINE_CN.md", ["#", "## 固化的自动流程"]),
    "H15": ("FMOtpiefinal_HIV1lab/reports/DEFAULT_GRAPH_PIPELINE_20261019.md", ["## CDK8", "## BACE", "## Galectin3"]),
    "H16": ("FMOtpiefinal_HIV1lab/reports/GRAPH_REPAIR_ENGINEERING_CONCLUSION_20261019.md", ["## 为什么之前失败", "坐标匹配2×2", "## 是否提高了排序"]),
    "H17": ("FMOtpiefinal_HIV1lab/reports/ROBUST_ACCELERATED_REPAIR_20261020.md", ["## 先纠正错误归因", "## 默认能力的修复"]),
    "H18": ("FMOtpiefinal_HIV1lab/reports/NON_QM_ACCELERATION_RESULTS_20261020.md", ["## Implementation", "## Full-workflow paired"]),
    "H19": ("FMOtpiefinal_HIV1lab/reports/STRUCTURE_UNIVERSALITY_AND_WEAK_TARGETS_20261019.md", ["## 最终完成结果", "## 仍未完成的普适性"]),
    "H20": ("FMOtpiefinal_HIV1lab/reports/WEAK_GRAPH_WATER_TRANSFER_20261019.md", ["## 结果", "## 为什么这轮没有提分"]),
    "H21": ("FMOtpiefinal_HIV1lab/reports/THREE_TARGET_TRANSFER_20261018.md", ["## 全量QM结果"]),
    "H22": ("FMOtpiefinal_HIV1lab/reports/RECEPTOR_CHEMICAL_GRAPH_20261019.md", ["#", "##"]),
    "H23": ("FMOtpiefinal_HIV1lab/reports/ORIGINAL87_GRAPH_PIPELINE_REPRODUCTION_20261019.md", ["## 最终结果", "### 成本"]),
    "H24": ("FMOtpiefinal_HIV1lab/reports/SPECIAL_CHEMISTRY_EXTENSION_20261019.md", ["#", "## FXa 历史结果"]),
    "H25": ("FMOtpiefinal_HIV1lab/reports/GRAPH_CHEMICAL_ANNOTATION_20261019.md", ["## 基线结果", "## 对流程的具体建议"]),
    "H26": ("FMOtpiefinal_HIV1lab/reports/PROTON_WATER_FACTORIAL_20261018.md", ["## 实验控制", "## 全量指标"]),
    "H27": ("FMOtpiefinal_HIV1lab/NEGATIVE_TRANSFER_AUDIT_20261017.md", ["#"]),
    "H28": ("FMOtpiefinal_HIV1lab/reports/QM_ROBUSTNESS_20261019.md", ["#"]),
    "H29": ("docs/TANKYRASE_GAS_PBSA_87_RESULTS_CN.md", ["## 全量执行和原始日志", "## 独立 87/87 诊断代理"]),
    "H30": ("docs/FMO_FINAL_ADAPTIVE_METHOD_CN.md", ["## 8. 独立 HF 三态复核"]),
    "H31": ("docs/FMO_ENGINEERING_PIPELINE_STATUS_CN.md", ["### 独立 HF 三态 QM 队列结项", "### 独立 HF 三态 QM（2026-10-04）"]),
}
JSON_SOURCES = {
    "J01": "runs/fmo87/ablation/tpie_protocol_transfer_20261014/evaluation.json",
    "J02": "runs/fmo87/ablation/E23_restraint_dose_analysis/summary.json",
    "J03": "runs/fmo87/ablation/fxa_hydrogen_full_ablation_20261015.json",
    "J04": "runs/fmo87/new_full_pipeline_20261015/complete_87_tpie.json",
    "J05": "runs/fmo87/robust_tpie_four_cohorts_score_only_20261016/evaluation.json",
    "J06": "FMOtpiefinal_HIV1lab/runs/stage_reuse_benchmark_20261020/verification.json",
    "J07": "FMOtpiefinal_HIV1lab/runs/carbonyl_graph_control_20261019/evaluation.json",
    "J08": "FMOtpiefinal_HIV1lab/runs/all_accelerated_robust_20261020/evaluation.json",
    "J09": "FMOtpiefinal_HIV1lab/runs/all_accelerated_robust_20261020/final_status.json",
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def excerpts(text: str, anchors: list[str]) -> list[dict]:
    lines = text.splitlines()
    selected = []
    for anchor in anchors:
        matched = False
        for index, line in enumerate(lines):
            if anchor in line:
                selected.append({"anchor": anchor, "start_line": index + 1,
                                 "end_line": min(index + 18, len(lines)),
                                 "text": "\n".join(lines[index:index + 18])})
                matched = True
                break
        if not matched:
            raise ValueError(f"historical evidence anchor missing: {anchor}")
    return selected


def build(workspace: Path, root: Path) -> dict:
    sources = {}
    values = {}
    for identifier, (relative, anchors) in REPORTS.items():
        path = workspace / relative
        sources[identifier] = {"workspace_relative_path": relative, "sha256": digest(path),
                               "kind": "historical_report", "excerpts": excerpts(path.read_text(), anchors)}
    for identifier, relative in JSON_SOURCES.items():
        path = workspace / relative
        sources[identifier] = {"workspace_relative_path": relative, "sha256": digest(path),
                               "kind": "machine_result"}
        values[identifier] = json.loads(path.read_text())
    metrics = {
        "schema": "fmotpiefinal-history-metrics-v1",
        "scope": "Read-only extraction from historical results, not new QM or causal reanalysis",
        "paired_protocol_comparisons": values["J01"]["comparisons"],
        "restraint_dose": {key: values["J02"][key] for key in
                           ("n", "median_heavy_rmsd_vs_k100_angstrom", "metrics")},
        "fxa_fixed_pose_hydrogen_arms": values["J03"]["arms"],
        "original87_full_preparation_20261015": {key: values["J04"][key] for key in ("n", "per_target", "size_weighted")},
        "softcap_20261016": values["J05"]["per_target"],
        "acceleration_pairs": [{key: row[key] for key in
                                ("name", "variant", "end_to_end_wall_seconds", "qm_wall_seconds",
                                 "non_qm_wall_seconds", "tpie", "independently_verified")}
                               for row in values["J06"]["records"]],
        "carbonyl_fixed_pose_control": {key: values["J07"][key] for key in
                                      ("n_crystals", "n_compounds", "n_new_qm", "no_case_excluded", "results")},
        "latest_2004_groups": {name: {key: group[key] for key in
                              ("n", "n_crystals", "n_compounds", "raw", "softcap") if key in group}
                              for name, group in values["J08"]["groups"].items()},
        "latest_2004_final_status": values["J09"],
    }
    cohort_names = ("CDK2", "JNK1", "P38", "TYK2")
    metrics["latest_original87_weighted"] = {
        score: {metric: sum(values["J08"]["groups"][name]["n"] *
                           values["J08"]["groups"][name][score][metric] for name in cohort_names) / 87
                for metric in ("pearson", "spearman", "kendall")}
        for score in ("raw", "softcap")}
    wet = metrics["paired_protocol_comparisons"]["wet_k10_vs_k100"]
    if wet["n"] != 87 or wet["weighted_after"]["pearson"] <= wet["weighted_before"]["pearson"]:
        raise ValueError("historical k10 comparison does not match the recorded evidence")
    if values["J09"]["n_verified"] != 909 or values["J09"]["n_requested"] != 909:
        raise ValueError("historical final status does not cover all 909 structures")
    if digest(workspace / JSON_SOURCES["J08"]) != values["J09"]["evaluation_sha256"]:
        raise ValueError("historical final evaluation fingerprint changed")
    metrics_path = root / "docs/HISTORY_METRICS.json"
    metrics_path.write_text(json.dumps(metrics, indent=2, ensure_ascii=False) + "\n")
    report = {"schema": "fmotpiefinal-history-evidence-v1", "date_utc": "2026-10-20",
              "workspace": "sophosqm-style-fmo-baseline (not bundled)",
              "git_history_limit": "Research workspace has no continuous Git history; this traces saved reports, results and hashes, not invented commits",
              "sources": sources, "metrics_file": "docs/HISTORY_METRICS.json",
              "metrics_sha256": digest(metrics_path)}
    (root / "docs/HISTORY_EVIDENCE.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    return {"status": "PASS", "source_count": len(sources), "metrics_sha256": report["metrics_sha256"]}


def verify(root: Path, workspace: Path | None = None) -> dict:
    evidence = json.loads((root / "docs/HISTORY_EVIDENCE.json").read_text())
    if digest(root / evidence["metrics_file"]) != evidence["metrics_sha256"]:
        raise ValueError("historical metric snapshot changed")
    if workspace is not None:
        for identifier, source in evidence["sources"].items():
            if digest(workspace / source["workspace_relative_path"]) != source["sha256"]:
                raise ValueError(f"historical source fingerprint changed: {identifier}")
    return {"status": "PASS", "metric_snapshot": "verified", "source_count": len(evidence["sources"]),
            "workspace_sources_verified": workspace is not None}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path)
    parser.add_argument("--refresh", action="store_true", help="explicitly regenerate from the original research workspace")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    if args.refresh and args.workspace is None:
        parser.error("--refresh requires --workspace")
    result = build(args.workspace, root) if args.refresh else verify(root, args.workspace)
    print(json.dumps(result, indent=2))
