"""Validate and summarize frozen multi-frame QuanPol QM/MM vs matched-MM runs."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re
import statistics

ENERGY_TERMS = (
    "bond_stretching_kcal_mol",
    "bond_angle_bending_kcal_mol",
    "dihedral_rotation_kcal_mol",
    "lj_repulsion_kcal_mol",
    "lj_dispersion_kcal_mol",
    "charge_kcal_mol",
)


def describe(values):
    values = [float(value) for value in values]
    return {
        "n": len(values),
        "mean": statistics.mean(values),
        "sample_sd": statistics.stdev(values) if len(values) > 1 else None,
        "median": statistics.median(values),
        "min": min(values),
        "max": max(values),
    }


def single_state_weight_ess(delta_u_values, temperature_kelvin=298.15):
    beta = 1.0 / (0.00198720425864083 * temperature_kelvin)
    minimum = min(delta_u_values)
    weights = [math.exp(-beta * (float(value) - minimum)) for value in delta_u_values]
    effective_sample_size = sum(weights) ** 2 / sum(weight * weight for weight in weights)
    finite_sample_free_energy = minimum - math.log(sum(weights) / len(weights)) / beta
    return effective_sample_size, beta, finite_sample_free_energy


def summarize(input_paths):
    reports = [json.loads(path.read_text()) for path in input_paths]
    if not reports or any(report.get("status") != "COMPLETE" for report in reports):
        raise ValueError("all input runs must be complete")
    if any(report.get("labels_used") for report in reports):
        raise ValueError("label-dependent frame selection is not permitted")
    if len({report.get("target") for report in reports}) != 1:
        raise ValueError("input runs must have the same target")
    if len({tuple(report.get("times_iterations", [])) for report in reports}) != 1:
        raise ValueError("input runs must use identical fixed frame times")
    swrb_values = {float(report.get("swrb_angstrom", 1001.0)) for report in reports}
    if len(swrb_values) != 1:
        raise ValueError("input runs must use the same SWRB setting")
    shell_values = {float(report.get("shell_angstrom", 12.0)) for report in reports}
    if len(shell_values) != 1:
        raise ValueError("input runs must use the same water/ion shell setting")

    times = list(reports[0]["times_iterations"])
    results = [row for report in reports for row in report["results"]]
    by_key = {}
    for row in results:
        key = (int(row["repeat"]), int(row["iteration"]), row["leg"], int(row["endpoint"]))
        if key in by_key:
            raise ValueError(f"duplicate frame result: {key}")
        if row["qm"].get("status") != "SUCCEEDED" or row["mm"].get("status") != "SUCCEEDED":
            raise ValueError(f"failed QM/MM or MM calculation: {key}")
        if row["qm"].get("scf_unconverged") or row["mm"].get("scf_unconverged"):
            raise ValueError(f"unconverged calculation: {key}")
        parity = float(row["mm_parity_difference_kcal_mol"])
        if abs(parity) > 0.1:
            raise ValueError(f"MM parity exceeds 0.1 kcal/mol: {key} ({parity})")
        by_key[key] = row

    repeats = sorted({key[0] for key in by_key})
    expected = {(repeat, iteration, leg, endpoint)
                for repeat in repeats for iteration in times
                for leg in ("bound", "solution") for endpoint in (0, 10)}
    if set(by_key) != expected:
        missing = sorted(expected - set(by_key))
        extra = sorted(set(by_key) - expected)
        raise ValueError(f"incomplete frame matrix; missing={missing}, extra={extra}")

    qm_center_by_key = {}
    for key, row in by_key.items():
        log_path = Path(row["qm"].get("log", ""))
        center_values = re.findall(
            r"QM CENTER\s+ENERGY =\s*([-+0-9.Ee]+)\s+KCAL/MOL",
            log_path.read_text(errors="replace")) if log_path.is_file() else []
        if not center_values:
            raise ValueError(f"missing QM-center kcal/mol record: {key} ({log_path})")
        qm_center_by_key[key] = float(center_values[-1])

    cycles = []
    for repeat in repeats:
        for iteration in times:
            def delta(leg, endpoint):
                return float(by_key[(repeat, iteration, leg, endpoint)]["delta_u_kcal_mol"])

            correction = delta("bound", 10) - delta("bound", 0) - delta("solution", 10) + delta("solution", 0)
            def term_difference(leg, endpoint, term):
                row = by_key[(repeat, iteration, leg, endpoint)]
                return float(row["qm"].get(term, 0.0)) - float(row["mm"].get(term, 0.0))

            term_cycle = {
                term: term_difference("bound", 10, term) - term_difference("bound", 0, term)
                      - term_difference("solution", 10, term) + term_difference("solution", 0, term)
                for term in ENERGY_TERMS
            }

            def residual(leg, endpoint, model):
                row = by_key[(repeat, iteration, leg, endpoint)][model]
                total = float(row["potential_kcal_mol"])
                return total - sum(float(row.get(term, 0.0)) for term in ENERGY_TERMS)

            residual_cycle = ((residual("bound", 10, "qm") - residual("bound", 10, "mm"))
                              - (residual("bound", 0, "qm") - residual("bound", 0, "mm"))
                              - (residual("solution", 10, "qm") - residual("solution", 10, "mm"))
                              + (residual("solution", 0, "qm") - residual("solution", 0, "mm")))
            center_cycle = (qm_center_by_key[(repeat, iteration, "bound", 10)]
                            - qm_center_by_key[(repeat, iteration, "bound", 0)]
                            - qm_center_by_key[(repeat, iteration, "solution", 10)]
                            + qm_center_by_key[(repeat, iteration, "solution", 0)])
            cycles.append({"repeat": repeat, "iteration": iteration,
                           "correction_kcal_mol": correction,
                           "term_breakdown_kcal_mol": term_cycle,
                           "qm_center_cycle_kcal_mol": center_cycle,
                           "qm_hidden_residual_kcal_mol": residual_cycle,
                           "breakdown_closure_kcal_mol": sum(term_cycle.values()) + residual_cycle - correction})

    repeat_means = {
        str(repeat): statistics.mean(row["correction_kcal_mol"] for row in cycles if row["repeat"] == repeat)
        for repeat in repeats
    }
    by_iteration = {
        str(iteration): describe(row["correction_kcal_mol"] for row in cycles if row["iteration"] == iteration)
        for iteration in times
    }
    by_state = {}
    weight_ess = {}
    for leg in ("bound", "solution"):
        for endpoint in (0, 10):
            key = f"{leg}_{endpoint}"
            values = [float(row["delta_u_kcal_mol"]) for row in results
                      if row["leg"] == leg and int(row["endpoint"]) == endpoint]
            by_state[key] = describe(values)
            ess, beta, finite_sample_free_energy = single_state_weight_ess(values)
            weight_ess[key] = {"n_frames": len(values), "effective_sample_size": ess,
                               "temperature_kelvin": 298.15, "beta_per_kcal_mol": beta,
                               "naive_finite_sample_zwanzig_kcal_mol": finite_sample_free_energy}

    naive_zwanzig_correction = (
        weight_ess["bound_10"]["naive_finite_sample_zwanzig_kcal_mol"]
        - weight_ess["bound_0"]["naive_finite_sample_zwanzig_kcal_mol"]
        - weight_ess["solution_10"]["naive_finite_sample_zwanzig_kcal_mol"]
        + weight_ess["solution_0"]["naive_finite_sample_zwanzig_kcal_mol"]
    )
    repeat_zwanzig_corrections = {}
    for repeat in repeats:
        state_free_energies = {}
        for leg in ("bound", "solution"):
            for endpoint in (0, 10):
                state = f"{leg}_{endpoint}"
                values = [float(row["delta_u_kcal_mol"]) for row in results
                          if int(row["repeat"]) == repeat and row["leg"] == leg
                          and int(row["endpoint"]) == endpoint]
                _, _, state_free_energies[state] = single_state_weight_ess(values)
        repeat_zwanzig_corrections[str(repeat)] = (
            state_free_energies["bound_10"] - state_free_energies["bound_0"]
            - state_free_energies["solution_10"] + state_free_energies["solution_0"]
        )

    corrections = [row["correction_kcal_mol"] for row in cycles]
    parities = [float(row["mm_parity_difference_kcal_mol"]) for row in results]
    job_walls = [float(row[model].get("wall_seconds", 0.0))
                 for row in results for model in ("qm", "mm")]
    zero_step_messages = 0
    energy_at_zero_step = 0
    audited_logs = 0
    for row in results:
        for model in ("qm", "mm"):
            log_path = Path(row[model].get("log", ""))
            if not log_path.is_file():
                continue
            text = log_path.read_text(errors="replace")
            audited_logs += 1
            at_zero = bool(re.search(r"OPT STEP\s*=\s*0\b", text))
            has_energy = "POTENTIAL                     ENERGY =" in text or "POTENTIAL                     ENERGY =" in text.replace("\r", "")
            failed_zero_step = "QUANPOL OPTIMIZATION FAILED" in text
            zero_step_messages += int(failed_zero_step)
            energy_at_zero_step += int(at_zero and has_energy and failed_zero_step)
    worst_cycle = max(cycles, key=lambda row: abs(row["correction_kcal_mol"]))
    report = {
        "schema": "quanpol-multiframe-summary-v1",
        "status": "COMPLETE",
        "target": reports[0]["target"],
        "shell_angstrom": next(iter(shell_values)),
        "swrb_angstrom": next(iter(swrb_values)),
        "labels_used": False,
        "frame_times_iterations": times,
        "independent_repeat_indices": repeats,
        "paired_single_points": len(results),
        "gamess_jobs": 2 * len(results),
        "sum_gamess_wall_seconds": sum(job_walls),
        "mean_paired_qm_mm_wall_seconds": statistics.mean(
            float(row["qm"].get("wall_seconds", 0.0)) + float(row["mm"].get("wall_seconds", 0.0))
            for row in results
        ),
        "fixed_frame_cycles": len(cycles),
        "all_scf_converged_and_gamess_normal": True,
        "all_matched_mm_parity_passed_0p1_kcal_mol": True,
        "mm_parity_abs_kcal_mol": describe(abs(value) for value in parities),
        "no_bond_warning_occurrences": sum(
            int(row[model].get("warnings_no_bond", 0)) for row in results for model in ("qm", "mm")
        ),
        "gamess_log_audit": {
            "logs_found": audited_logs,
            "zero_step_optimizer_fail_messages": zero_step_messages,
            "potential_recorded_at_zero_step_before_message": energy_at_zero_step,
        },
        "cycle_correction_kcal_mol": describe(corrections),
        "cycle_correction_by_iteration": by_iteration,
        "cycle_mean_by_repeat_kcal_mol": repeat_means,
        "between_repeat_sd_of_mean_kcal_mol": statistics.stdev(repeat_means.values()) if len(repeat_means) > 1 else None,
        "term_cycle_mean_kcal_mol": {
            term: statistics.mean(row["term_breakdown_kcal_mol"][term] for row in cycles)
            for term in ENERGY_TERMS
        },
        "term_cycle_statistics_kcal_mol": {
            term: describe(row["term_breakdown_kcal_mol"][term] for row in cycles)
            for term in ENERGY_TERMS
        },
        "qm_hidden_residual_cycle_kcal_mol": describe(row["qm_hidden_residual_kcal_mol"] for row in cycles),
        "qm_center_cycle_kcal_mol": describe(row["qm_center_cycle_kcal_mol"] for row in cycles),
        "qm_center_nonzero_single_points": sum(
            abs(value) > 1e-8 for value in qm_center_by_key.values()
        ),
        "max_abs_breakdown_closure_kcal_mol": max(abs(row["breakdown_closure_kcal_mol"]) for row in cycles),
        "largest_absolute_cycle": worst_cycle,
        "delta_u_by_state_kcal_mol": by_state,
        "single_state_zwanzig_weight_ess_diagnostic": weight_ess,
        "naive_finite_sample_zwanzig_binding_correction_kcal_mol": naive_zwanzig_correction,
        "naive_finite_sample_zwanzig_correction_by_repeat_kcal_mol": repeat_zwanzig_corrections,
        "cycles": cycles,
        "source_runs": [{"path": str(path.resolve()), "driver_sha256": item.get("driver_sha256")}
                        for path, item in zip(input_paths, reports, strict=True)],
        "interpretation_limits": [
            "These are sparse fixed-frame potential-energy double differences, not a binding free energy or corrected FEP result.",
            f"Single-state exponential-weight ESS was computed from {len(repeats)} repeats and {len(times)} fixed frames per repeat, without autocorrelation correction; this is only a diagnostic, not validated free-energy reweighting.",
            "The finite-shell correction is not a correction to full periodic QM/MM. Reweighting a PME trajectory would require explicitly defining a subtractive target H_target = H_MM,PME + delta_local, then validating ensemble overlap and convergence for that approximate model.",
            "The matched MM parity validates the reconstructed finite-shell FFDATA model, not equivalence to periodic PME/OpenFE.",
            "FFDATA POL values are zero and IDOPOL=0; the calculation does not include mutual induced-dipole polarization in the MM environment.",
            "No cycle or repeat was removed; all frozen repeats and frame times are included.",
        ],
    }
    return report


def render_markdown(report):
    cycle = report["cycle_correction_kcal_mol"]
    lines = [
        f"# {report['target']} QM/MM 多帧匹配 MM 诊断",
        "",
        "此报告汇总预先固定的各条轨迹与采样时点；不使用实验亲和力选帧，也不剔除离群帧。",
        "",
        "## 完成与能量闭合",
        "",
        f"- 轨迹重复：{len(report['independent_repeat_indices'])} 条；采样时点：`{report['frame_times_iterations']}`。",
        f"- 水/离子环境半径：{report['shell_angstrom']:g} Å（蛋白残基全保留）。",
        f"- 完成成对 QM/MM 与 MM 单点：{report['paired_single_points']}；四态双差值：{report['fixed_frame_cycles']}。",
        f"- QuanPol `SWRB`：{report['swrb_angstrom']:g} Å；QM-center 罚能非零的 QM 单点数：{report['qm_center_nonzero_single_points']}。",
        f"- GAMESS 作业数：{report['gamess_jobs']}；各作业 wall time 累计 {report['sum_gamess_wall_seconds'] / 60:.1f} 分钟（串行总运行时间另含输入生成与分析）。",
        f"- QM/MM 与 MM 均正常结束、SCF 均收敛：{report['all_scf_converged_and_gamess_normal']}。",
        f"- FFDATA 匹配 MM 对照最大绝对误差：{report['mm_parity_abs_kcal_mol']['max']:.4f} kcal/mol。",
        "- 计算对象：DFTB3 QM 配体 + 固定电荷 FFDATA MM 蛋白/环境；`IDOPOL=0` 且极化参数 `POL=0`，不含 MM 诱导偶极的互相极化。",
        "- `exp(−βΔU)` 单向权重 ESS 仅作诊断（298.15 K，未校正帧自相关）。来源轨迹使用周期 PME，本修正不等于完整周期 QM/MM；若明确采用 `H_target = H_MM,PME + δ_local`，可定义局部替换目标，但仍必须验证重叠与收敛。本轮不能把直接 Zwanzig 值当作已验证的自由能。",
        f"- 日志复核：{report['gamess_log_audit']['logs_found']} 份；其中 {report['gamess_log_audit']['potential_recorded_at_zero_step_before_message']} 份先记录 `OPT STEP=0` 势能，随后打印 `QUANPOL OPTIMIZATION FAILED`。这些是零步单点，不是几何优化成功；判读仅使用已记录的单点能量。",
        f"- `FORMS NO BOND` 告警累计 {report['no_bond_warning_occurrences']:,} 次；该数是 QM 与 MM 日志中的原子级累计告警，不等于失败作业数；包括零力常数约束键与离子等拓扑表示，不能简化称为‘水被丢失’。",
        "",
        "## 固定帧双差值",
        "",
        "| 汇总 | n | 均值 | SD | 中位数 | 范围 |",
        "|---|---:|---:|---:|---:|---:|",
        f"| 全部双差值 (kcal/mol) | {cycle['n']} | {cycle['mean']:.3f} | {cycle['sample_sd']:.3f} | {cycle['median']:.3f} | {cycle['min']:.3f} 至 {cycle['max']:.3f} |",
    ]
    for iteration, values in report["cycle_correction_by_iteration"].items():
        lines.append(f"| iteration {iteration} | {values['n']} | {values['mean']:.3f} | {values['sample_sd']:.3f} | {values['median']:.3f} | {values['min']:.3f} 至 {values['max']:.3f} |")
    lines.extend([
        "",
        "### 单向权重诊断（不是自由能结果）",
        "",
        "| 状态 | ESS / 帧数 |",
        "|---|---:|",
    ])
    lines.extend(
        f"| {state} | {values['effective_sample_size']:.2f}/{values['n_frames']} |"
        for state, values in report["single_state_zwanzig_weight_ess_diagnostic"].items()
    )
    state_frame_count = next(iter(report["single_state_zwanzig_weight_ess_diagnostic"].values()))["n_frames"]
    repeats = len(report["independent_repeat_indices"])
    ess_values = [entry["effective_sample_size"]
                  for entry in report["single_state_zwanzig_weight_ess_diagnostic"].values()]
    lines.extend([
        "",
        f"用每状态 {state_frame_count} 个未作自相关校正的固定帧直接计算的有限样本 Zwanzig 四态差为 `{report['naive_finite_sample_zwanzig_binding_correction_kcal_mol']:.3f}` kcal/mol；各重复分别为 `{report['naive_finite_sample_zwanzig_correction_by_repeat_kcal_mol']}`。状态 ESS 范围 `{min(ess_values):.2f}–{max(ess_values):.2f}/{state_frame_count}`（{repeats} 条重复、每条 {len(report['frame_times_iterations'])} 帧），显示权重集中且样本很稀疏。这些数值仅用于暴露采样/重加权风险，不是已收敛的自由能修正；即使明确采用局部替换 Hamiltonian，这一采样限制仍然存在。",
    ])
    lines.extend(["", "### 每条重复的全部双差值", "",
                  "| repeat | " + " | ".join(str(iteration) for iteration in report["frame_times_iterations"]) + " |",
                  "|---|" + "---:|" * len(report["frame_times_iterations"])])
    cycle_lookup = {(row["repeat"], row["iteration"]): row["correction_kcal_mol"] for row in report["cycles"]}
    for repeat in report["independent_repeat_indices"]:
        cells = [f"{cycle_lookup[(repeat, iteration)]:.3f}" for iteration in report["frame_times_iterations"]]
        lines.append(f"| {repeat} | " + " | ".join(cells) + " |")
    lines.extend(["", "能量项分解（先逐态取 QM/MM−MM，再对四态做双差值）：", "",
                  "| 项 | 跨帧平均贡献 (kcal/mol) | SD |", "|---|---:|---:|"])
    lines.extend(
        f"| {term} | {values['mean']:.3f} | {values['sample_sd']:.3f} |"
        for term, values in report["term_cycle_statistics_kcal_mol"].items()
    )
    residual = report["qm_hidden_residual_cycle_kcal_mol"]
    lines.append(f"| QM/未显式列项残差 | {residual['mean']:.3f} | {residual['sample_sd']:.3f} |")
    worst = report["largest_absolute_cycle"]
    lines.extend([
        "",
        f"最大绝对双差值出现在 repeat {worst['repeat']} / iteration {worst['iteration']}：{worst['correction_kcal_mol']:.3f} kcal/mol；其中未显式列项残差为 {worst['qm_hidden_residual_kcal_mol']:.3f} kcal/mol。该残差合并 DFTB 电子能及未展开的 QM/MM 相互作用等项；若 QM-center 罚能非零，它也包含在残差中，因此本项必须与单独的 QM-center 审计一起解读。此帧保留在全部统计中。",
    ])
    lines.extend(["", "## 重复间均值", ""])
    lines.extend(f"- repeat {repeat}: {value:.3f} kcal/mol" for repeat, value in report["cycle_mean_by_repeat_kcal_mol"].items())
    lines.extend([
        "",
        "## 解释",
        "",
        "结果只能说明在当前有限壳层 QM/MM 与匹配 FFDATA MM 势能模型下，稀疏抽帧双差值有多大波动。数值不是 ΔG，也不能直接作为 Spearman 排序特征。即使简单均值接近零，也不能解释为 QM 修正已抵消；均值会掩盖明显的帧间与重复间差异。",
        "",
        "要得到可与 FEP 对照的 QM/MM 自由能修正，必须先定义参考与目标 Hamiltonian（包括周期边界/长程静电、溶剂与 QM/MM 阻尼），再在 bound 与 solution 两端检查能量分布重叠、有效样本量、时间相关及重复收敛；优先用双向计算进一步交叉验证。本次每态仅有稀疏固定帧，不能支撑该自由能结论。",
        "",
        "原始逐帧输入、QM/MM/MM 日志、能量闭合值和坐标清单保存在对应 `runs/CDK8_qmmm_multiframe_20261027*` 目录。",
        "",
    ])
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, nargs="+", required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path, required=True)
    args = parser.parse_args()
    report = summarize(args.input)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(report, indent=2) + "\n")
    args.output_md.parent.mkdir(parents=True, exist_ok=True)
    args.output_md.write_text(render_markdown(report))
    print(json.dumps({key: report[key] for key in (
        "status", "paired_single_points", "fixed_frame_cycles", "cycle_correction_kcal_mol",
        "mm_parity_abs_kcal_mol")}, indent=2))


if __name__ == "__main__":
    main()
