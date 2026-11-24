"""Preserve full numerical-control failures and same-geometry disagreement."""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import re

from run_cdk8_solution_cycle import digest, write_json


def residual_cycles(path):
    pattern = re.compile(r"Iter=\s*(\d+) iFrag=\s*(\d+) EFMOc=\s*([-\d.]+) EFMOu=\s*([-\d.]+)")
    records = defaultdict(list)
    last = 0
    with path.open(errors="replace") as handle:
        for line in handle:
            match = pattern.search(line)
            if match:
                iteration, fragment, energy = int(match[1]), int(match[2]), float(match[3])
                last = max(last, iteration)
                records[fragment].append((iteration, energy))
    rows = []
    for fragment, values in records.items():
        recent = [(iteration, energy) for iteration, energy in values if iteration > last-6]
        if len(recent) < 2:
            continue
        energies = [energy for _, energy in recent]
        rows.append({"fragment": fragment, "last_six_cycles": recent,
            "energy_range_hartree": max(energies)-min(energies)})
    return sorted(rows, key=lambda row: row["energy_range_hartree"], reverse=True)[:5]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    arms = ("fxa_strict_numeric", "fxa_strict_diis", "fxa_tight_soscf")
    groups = {}
    lines = ["# FXa 宽侧链候选：数值收敛对照", "",
        "这是新侧链候选的数值诊断，不是原生产基线失败。所有对照保持几何、原子、GLN61、口袋、分片与电荷。", "",
        "| 数值对照 | 完成 | 原成功例最大同几何 tPIE 差 kcal/mol | 0.1 一致性门槛 | 失败例 |",
        "|---|---:|---:|---|---|"]
    for arm in arms:
        root = args.run_root / arm
        status = json.loads((root / "final_status.json").read_text())
        differences = status["same_geometry_tpie_differences_kcal_mol"]
        maximum = max(map(abs, differences.values())) if differences else None
        maxima = "—" if maximum is None else f"{maximum:.3f}"
        lines.append(f"| {arm} | {status['n_succeeded']}/{status['n_requested']} | {maxima} | {'通过' if status['same_geometry_agreement_passed'] else '未通过'} | {', '.join(status['failures'])} |")
        groups[arm] = {"status": status, "protocol_sha256": digest(root / "frozen_protocol.json"),
            "comparison_uses_all_source_successful_cases": True,
            "eligible_as_identical_numerical_substitute": status["status"] == "COMPLETE" and status["same_geometry_agreement_passed"]}
    root = args.run_root / "fxa_tight_soscf"
    cycles = residual_cycles(root / "5b/run/gamess.log")
    groups["fxa_tight_soscf"]["5b_last_cycle_fragment_fluctuations"] = cycles
    lines += ["", "## 已确认与尚未确认", "",
        "- 原 B 侧链候选 5e/5g 的失败可被 CNVDMP=0 消除，但同一设置在 5b 又产生周期残差，不是完整修复。",
        "- DIIS+DAMPH 在 27 例均被严格验收拒绝；日志记录 GLN61 单体持续不收敛，不删除该残基绕过。",
        "- 收紧内层 SCF 至 1e−10 仍有外层周期，说明不能简单归因为原内层阈值偏松。",
        "- 紧 SCF 的 5b 最近周期中 GLN61（片段4）能量反复切换约 1.1e−4 Hartree，明显大于外层 1e−7 阈值。",
        "- 数值设置还使部分原成功例 tPIE 改变超过预设一致性门槛，不能把不同数值结果无条件拼接。",
        "- 片段变化是关联证据；尚不足以判定是电子态切换、PCM响应、收敛器实现还是其它机制的唯一根因。",
        "- 不报告成功子集亲和力指标，不放宽阈值，不取相关系数最高的一组；生产默认未改变。"]
    args.output.write_text("\n".join(lines)+"\n")
    write_json(args.output.with_suffix(".json"), {"groups": groups, "production_default_changed": False})
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
