"""Publish full-cohort cycle evidence without claiming an unmet accuracy goal."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from run_cdk8_solution_cycle import digest, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.run.resolve()
    report = json.loads((root / "final_evaluation.json").read_text())
    native = json.loads((root / "native/execution.json").read_text())
    if len(report["rows"]) != 10 or len(native) != 10 or not all(value["status"] == "SUCCEEDED" for value in native.values()):
        raise ValueError("complete ten-ligand native and independent PCM states required")
    evidence = {"date_utc": "2026-10-21", "n": 10,
        "goal_spearman": 0.6, "goal_achieved": False, "production_changed": False,
        "metrics": report["metrics"], "posthoc_component_diagnostic_metrics": report["posthoc_component_diagnostic_metrics"],
        "rows": report["rows"], "failed_gas_jobs": report["failed_jobs"],
        "native_max_abs_reconciliation_kcal_mol": max(abs(value["native_minus_independent"]) for value in native.values()),
        "sources_sha256": {str(root / "final_evaluation.json"): digest(root / "final_evaluation.json"),
            str(root / "native/execution.json"): digest(root / "native/execution.json")},
        "limitations": report["limitations"], "sampling": None}
    sample = root / "solution_sampling4/final_evaluation.json"
    if not sample.exists():
        sample = root / "solution_sampling4/evaluation.json"
    if sample.exists():
        sampling = json.loads(sample.read_text())
        evidence["sampling"] = sampling
        evidence["sources_sha256"][str(sample)] = digest(sample)
        evidence["metrics"].update(sampling["metrics"])
        geometry = root / "solution_sampling4/geometry_audit.json"
        if geometry.exists():
            evidence["sampling_geometry_audit"] = json.loads(geometry.read_text())
            evidence["sources_sha256"][str(geometry)] = digest(geometry)
    attempt_records = []
    initial = root.parent / "CDK8_solution_cycle_20261021"
    gas_root = root.parent / "CDK8_solution_cycle_gas_corrected_20261021"
    for directory in (initial, root, gas_root):
        for path in directory.rglob("result.json"):
            value = json.loads(path.read_text())
            attempt_records.append({"path": str(path), "status": value.get("status"),
                "wall_seconds": value.get("wall_seconds"), "ncores": 2})
    evidence["runtime"] = {"n_recorded_attempt_results": len(attempt_records),
        "sum_recorded_job_wall_seconds": sum(value["wall_seconds"] or 0 for value in attempt_records),
        "requested_core_hours_with_recorded_wall_not_measured_cpu": sum((value["wall_seconds"] or 0) * value["ncores"] for value in attempt_records) / 3600,
        "warning": "Includes failed and exploratory attempts; native driver did not record per-job wall time, so this is incomplete accounting, not end-to-end preparation or measured CPU usage"}
    evidence["goal_achieved"] = any(value["spearman"] >= 0.6 for value in evidence["metrics"].values())
    truth = np.array([row["experimental_dg"] for row in evidence["rows"]])
    baseline = np.array([row["tpie_pcm"] for row in evidence["rows"]])
    random = np.random.default_rng(20261021)
    permutations = [random.permutation(len(truth)) for _ in range(10000)]
    observed = float(spearmanr(truth, baseline).statistic)
    null = [float(spearmanr(truth, baseline[indices]).statistic) for indices in permutations]
    evidence["baseline_spearman_permutation_p_two_sided"] = (1 + sum(abs(value) >= abs(observed) - 1e-12 for value in null)) / (len(null) + 1)
    lines = ["# CDK8：溶液态FMO与三态参照验证", "", "## 结论", "",
        "**已完成一致PCM三态和native交叉验证，但尚未达到Spearman约0.6的目标，不能称为提分成功。**",
        "当前生产主流程与分数没有被替换；没有删配体、拟合权重、加入分子量或按实验选择构象。",
        "本轮复用现行完整结构处理流水线的制备结果和k10姿态，属于匹配能量实验，而非重新制备原始结构。", "",
        "## 完整10例指标", "", "|方案|Pearson|Spearman|Kendall|", "|---|---:|---:|---:|"]
    lines.extend(f"|{arm}|{value['pearson']:.3f}|{value['spearman']:.3f}|{value['kendall']:.3f}|" for arm, value in evidence["metrics"].items())
    lines.extend(["|发布FEP+21-4对照|0.663|0.648|0.378|", "",
        "FEP+行为不同且并非本轮重跑；仅为既有匹配集合参照。native interaction与原tPIE几乎一致，",
        "因此native方法名本身并没有带来一个更强的新相互作用分数。", "", "## 找到并修正的问题", "",
        "### 1. PCM开关相同，实际计算项却不一致", "",
        "GAMESS在多空腔时可自动关闭IDISP。原AB和受体触发此逻辑，独立配体未触发。",
        "以12为例，独立配体带溶剂色散−35.973与排斥+5.866 kcal/mol，native配体两者均为0；",
        "约−30.107 kcal/mol的差异解释初轮三态不闭合，不能把初轮结果解释为物理上的负结论。",
        "源码证据为GAMESS source/pcmcav.src的NCAV2分支，日志明确打印Setting IDISP=0。",
        "修订版对全部三态显式IDISP=0，保留PCM静电与空腔项。31项PCM作业全部成功，",
        f"再跑10项native，最大闭合误差为{evidence['native_max_abs_reconciliation_kcal_mol']:.6f} kcal/mol（汇总值三位小数）。",
        "新增通用实际能量库存审计，拒绝输入开关相同而实际组成不同的三态。",
        "这是计算一致性修复，不是已经证明该溶剂模型准确，也不意味着tPIE主流程存在同样的三态做差问题。", "",
        "### 2. 初轮气相提取遗漏色散控制", "",
        "已修复：移除PCM时必须保留$DFT DC/IDCVER。重新运行匹配21项气相作业，",
        "但10个完整AB仍有二体SCF不收敛；不把初轮模型不匹配的成功作业混入。",
        "因此本轮没有发布完整气相三态或依赖它的溶剂增量排序。失败不代表数值为0。", "",
        "### 3. 游离优化不是自动改善", "",
        "单起点游离配体优化全部完成。17号局部应变约25.08 kcal/mol，其余约5.44–7.79；",
        "17号对齐后全原子位移约0.93 Å。已检查预期共价键，未发现超过粗筛阈值的明显断键，",
        "但这不是全面的化学态/驻点验证，也不足以断言这25 kcal/mol是真实结合应变。",
        "复合物使用MM k10几何，游离态使用QM最低点，混合几何模型可引入额外应变偏差；",
        "因此不能把大应变简单归因于实验配体本身。", "", "## 具体错排例子", "",
        "13与15：实验DG分别约−11.20与−9.72，13强约1.48 kcal/mol；但tPIE分别−61.23与−67.03，方向相反。",
        "native去溶剂化分别+24.499与+16.430，加入后进一步偏向15；对应三态约−60.89与−74.08。",
        "相互作用差主要来自ASP173片段约−4.740、LYS52片段约+3.055等贡献，不能只看一种接触。",
        "这些是数值分解证据，不足以证明ASP173或某个水应该被删除，也不能证明差异全由采样或DFTB误差造成。", "",
        "## 成分消融诊断", "", "下列为事后能量成分检查，不是新增盲验证的主要候选。", "",
        "|组合|Pearson|Spearman|Kendall|", "|---|---:|---:|---:|"])
    lines.extend(f"|{arm}|{value['pearson']:.3f}|{value['spearman']:.3f}|{value['kendall']:.3f}|" for arm, value in evidence["posthoc_component_diagnostic_metrics"].items())
    lines.extend(["", "均未达到目标。PCM极化状态上减去空腔项只是成分诊断，不等于重新计算另一种严格自由能。", "",
        "## 多起点与未完成项", ""])
    if evidence["sampling"]:
        sampling = evidence["sampling"]
        lines.append(f"固定种子4起点/配体，共40项游离态QM优化；失败起点：{sampling['failed_starts']}。")
        if not sampling["failed_starts"]:
            lines.append("40个起点全部成功，其中13/c3通过保留原失败的数值重试收敛；最低能量构象代理不是熵或自由能。")
            lines.append("tPIE加多起点应变Spearman为0.261，仅比基线0.236略高，未接近0.6；三态多起点为−0.030。")
            lines.append("40个优化结果的预期共价键距离粗筛无明显断键；这不是完整化学状态或稳定性证明。")
        else:
            lines.append("失败起点未通过时不发布该组完整排序；不以其它成功起点替代预定覆盖。最低能量构象代理不是熵或自由能。")
    lines.extend(["", "尚未得到Spearman≥0.6方案，因此未进行候选默认变更及Galectin3跨集提分宣称。",
        "下一步需要同模型的结合态/游离态几何对照，以及关键极性接触的更高层级QM验证，",
        "随后才有依据决定是否投入匹配的溶液/结合态采样或自由能修正。不能靠继续随意叠加能量项保证0.6。", "",
        "## 复核材料", "", "协议：CDK8_SOLUTION_CYCLE_PLAN.zh-CN.md；执行、数值重试、native、多起点与报告脚本位于tools。",
        "CDK8_SOLUTION_CYCLE_RESULTS.json保留逐例分数、完整指标、原始结果SHA256和覆盖状态。"])
    lines.extend(["", "## 成本口径", "", f"记录到{len(attempt_records)}个结果文件（含失败和探索性重试），有墙钟字段的作业累计{evidence['runtime']['sum_recorded_job_wall_seconds']/60:.2f}分钟。",
        "native执行器缺逐作业墙钟字段，故不能把这个不完整作业累计当作全流程耗时或实测CPU核时。"])
    args.output.mkdir(parents=True, exist_ok=True)
    write_json(args.output / "CDK8_SOLUTION_CYCLE_RESULTS.json", evidence)
    (args.output / "CDK8_SOLUTION_CYCLE_RESULTS.zh-CN.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({"goal_achieved": evidence["goal_achieved"], "n": 10, "sampling": evidence["sampling"]["metrics"] if evidence["sampling"] else None}))


if __name__ == "__main__":
    main()
