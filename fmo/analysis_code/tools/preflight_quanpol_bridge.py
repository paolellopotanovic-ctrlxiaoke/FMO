"""Estimate first-step weight concentration, never substitute it for bridge sampling."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from report_quanpol_multiframe import single_state_weight_ess, summarize
from run_matched_openfe import digest, write_json


def weights(values, temperature_kelvin, coupling=1.0):
    if not values or not all(math.isfinite(value) for value in values):
        raise ValueError("finite nonempty energies required")
    if not math.isfinite(temperature_kelvin) or temperature_kelvin <= 0:
        raise ValueError("finite positive temperature required")
    if not math.isfinite(coupling) or not 0 <= coupling <= 1:
        raise ValueError("coupling must lie in [0, 1]")
    beta = 1 / (.00198720425864083 * temperature_kelvin)
    minimum = min(values)
    raw = [math.exp(-beta * coupling * (value - minimum)) for value in values]
    total = sum(raw)
    return [value / total for value in raw]


def first_step(values, temperature_kelvin, minimum_ess_fraction=.8):
    if not 0 < minimum_ess_fraction <= 1:
        raise ValueError("ESS fraction must be in (0, 1]")

    def fraction(coupling):
        normalized = weights(values, temperature_kelvin, coupling)
        return 1 / sum(value**2 for value in normalized) / len(values)

    if fraction(1) >= minimum_ess_fraction:
        return 1.0
    lower, upper = 0.0, 1.0
    for _ in range(60):
        middle = (lower + upper) / 2
        if fraction(middle) >= minimum_ess_fraction:
            lower = middle
        else:
            upper = middle
    return lower


def analyze(rows, temperature_kelvin, minimum_ess_fraction=.8):
    states = {}
    for leg in ("bound", "solution"):
        for endpoint in (0, 10):
            selected = [row for row in rows if row["leg"] == leg and row["endpoint"] == endpoint]
            values = [float(row["delta_u_kcal_mol"]) for row in selected]
            normalized = weights(values, temperature_kelvin)
            by_repeat = {}
            for repeat in sorted({row["repeat"] for row in selected}):
                subset = [row["delta_u_kcal_mol"] for row in selected if row["repeat"] == repeat]
                ess, _, _ = single_state_weight_ess(subset, temperature_kelvin)
                by_repeat[str(repeat)] = {"frames": len(subset), "weight_ess_uncorrected": ess,
                                         "pooled_weight_fraction": sum(weight for row, weight in zip(selected, normalized, strict=True) if row["repeat"] == repeat)}
            states[f"{leg}_{endpoint}"] = {
                "frames": len(values), "weight_ess_uncorrected": 1 / sum(value**2 for value in normalized),
                "largest_normalized_weight": max(normalized), "repeat_diagnostics": by_repeat,
                "initial_coupling_for_requested_ess_fraction": first_step(values, temperature_kelvin, minimum_ess_fraction),
                "all_frame_weights": [{**{field: row[field] for field in ("repeat", "iteration")}, "weight": weight}
                                      for row, weight in zip(selected, normalized, strict=True)],
            }
    return states


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--temperature", type=float, required=True)
    parser.add_argument("--ess-fraction", type=float, default=.8)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    validation = summarize([args.source])
    if validation["qm_center_nonzero_single_points"]:
        raise ValueError("center-energy contamination in source")
    source = json.loads(args.source.read_text())
    states = analyze(source["results"], args.temperature, args.ess_fraction)
    report = {"status": "DIAGNOSTIC_ONLY_NOT_CONVERGED_FREE_ENERGY", "labels_used": False,
              "source_sha256": digest(args.source), "driver_sha256": digest(Path(__file__)),
              "temperature_kelvin": args.temperature, "requested_initial_ess_fraction": args.ess_fraction,
              "states": states,
              "conservative_common_first_coupling": min(state["initial_coupling_for_requested_ess_fraction"] for state in states.values()),
              "cautions": [
                  "Coupling values are diagnostic hints for the first step from the reference ensemble, not a complete lambda schedule.",
                  "No intermediate/target trajectories were sampled; this does not validate overlap for any subsequent bridge step.",
                  "ESS is not autocorrelation-corrected. The ESS-fraction threshold is an engineering heuristic, not a universal accuracy criterion.",
                  "Required before bridge dynamics: validate target forces, fixed-constraint manifold, periodic reference energy, and boundary smoothness.",
              ]}
    args.report.parent.mkdir(parents=True, exist_ok=True)
    write_json(args.report.with_suffix(".json"), report)
    lines = ["# CDK8 局部 QM/MM 桥接采样前置诊断", "",
             f"温度 {args.temperature:g} K；只使用现有全部冻结帧，不使用实验亲和力。", "",
             f"| 状态 | 帧数 | 总修正 ESS | 最大单帧权重 | 首步 λ 提示 (ESS/N≥{args.ess_fraction:g}) |", "|---|---:|---:|---:|---:|"]
    for state, details in states.items():
        lines.append(f"| {state} | {details['frames']} | {details['weight_ess_uncorrected']:.3f} | {details['largest_normalized_weight']:.1%} | {details['initial_coupling_for_requested_ess_fraction']:.4f} |")
    lines.extend(["", "这些 λ 仅表示从现有参考样本出发的第一步权重集中程度；不是已验证的完整桥接方案。不能把同一批 MM 帧按多个 λ 重新加权，当成在中间态进行了采样。",
                  f"阈值 {args.ess_fraction:g} 是事先设定的工程启发式，不是普适收敛判据；本轮未新增中间态轨迹。", "",
                  "下一步动力学必须先验证局部目标力、约束空间、周期 MM 基底与壳层边界连续性。当前硬壳层单点脚本不能直接当作有效的 QM/MM MD 引擎。", ""])
    args.report.write_text("\n".join(lines))


if __name__ == "__main__":
    main()
