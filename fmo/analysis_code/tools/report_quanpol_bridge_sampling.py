"""Report actual bridge trajectories without mistaking short pilots for convergence."""
from __future__ import annotations

import argparse
import fcntl
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import time

import numpy as np

from run_matched_openfe import write_json
from sample_quanpol_bridge import GAS_KJ, KJ_PER_KCAL


def normalized_log_weights(values):
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or len(values) == 0 or not np.isfinite(values).all():
        raise ValueError("finite nonempty vector required")
    shifted = values - np.max(values)
    weights = np.exp(shifted)
    weights /= weights.sum()
    return weights, float(np.max(values) + math.log(np.exp(shifted).mean()))


def adjacent_diagnostic(first, second, first_lambda, second_lambda, temperature, offset):
    if not first or not second:
        return None
    thermal = GAS_KJ * temperature / KJ_PER_KCAL
    spacing = second_lambda - first_lambda
    forward, forward_log_mean = normalized_log_weights(-spacing * (np.asarray(first) - offset) / thermal)
    reverse, reverse_log_mean = normalized_log_weights(spacing * (np.asarray(second) - offset) / thermal)
    forward_centered = -thermal * forward_log_mean
    reverse_centered = thermal * reverse_log_mean
    return {"lambda_from": first_lambda, "lambda_to": second_lambda,
            "frames_forward": len(first), "frames_reverse": len(second),
            "forward_weight_ess_uncorrected": float(1 / np.sum(forward**2)),
            "reverse_weight_ess_uncorrected": float(1 / np.sum(reverse**2)),
            "forward_largest_weight": float(np.max(forward)), "reverse_largest_weight": float(np.max(reverse)),
            "forward_step_centered_kcal_mol": forward_centered, "reverse_step_centered_kcal_mol": reverse_centered,
            "difference_forward_reverse_kcal_mol": forward_centered - reverse_centered,
            "constant_offset_step_kcal_mol": spacing * offset,
            "diagnostic_only_not_equilibrated_or_autocorrelation_corrected": True}


def summarize(root):
    summary = {"generated_utc": datetime.now(timezone.utc).isoformat(), "status": "RUNNING",
               "labels_used": False, "free_energy_converged": False,
               "full_four_state_free_energy_correction_kcal_mol": None, "new_target_ranking_metrics": None,
               "scope": "CDK8 edge_000 four-state local-QMMM bridge pilot, not full-target affinity ranking",
               "states": {}, "warnings": ["Weight ESS is not time-series effective sample size.",
                    "Short-window high acceptance/overlap does not establish equilibration or convergence.",
                    "The target is a local-QM correction to periodic MM at fixed volume; not periodic QM/MM or final NPT binding free energy.",
                    "Rejected Metropolis proposals retain the current state and remain in the sample sequence."]}
    for state in ("b0", "b10", "s0", "s10"):
        directory = root / state
        protocol_path = directory / "protocol.json"
        if not protocol_path.exists():
            summary["states"][state] = {"status": "PREPARING_OR_QUEUED"}
            continue
        protocol = json.loads(protocol_path.read_text())
        if "lambdas" not in protocol:
            summary["states"][state] = {"status": protocol.get("status", "ERROR"), "error": protocol.get("error", "incomplete startup protocol")}
            continue
        samples_path = directory / "samples.jsonl"
        lines = samples_path.read_text().splitlines(keepends=True) if samples_path.exists() else []
        samples = [json.loads(line) for line in lines if line.endswith("\n")]
        windows = []
        production_by_lambda = {}
        for coupling in protocol["lambdas"]:
            rows = [row for row in samples if row["lambda"] == coupling]
            production = [row for row in rows if not row["burnin"]]
            values = [row["current_delta_kcal_mol"] for row in production]
            production_by_lambda[coupling] = values
            thermal_kj = GAS_KJ * protocol["temperature_kelvin"]
            mm_errors = np.array([row["mm_integration_energy_error_kj_mol"] for row in rows])
            quantum_changes = np.array([coupling * KJ_PER_KCAL * (row["proposal_delta_kcal_mol"] - row["previous_delta_kcal_mol"]) for row in rows])
            windows.append({"lambda": coupling, "attempted": len(rows), "accepted": sum(row["accepted"] for row in rows),
                            "acceptance_fraction": sum(row["accepted"] for row in rows) / len(rows) if rows else None,
                            "production_samples_including_rejections": len(values),
                            "production_accepted_transitions": sum(row["accepted"] for row in production),
                            "production_distinct_correction_values_not_independent_samples": len(set(values)),
                            "delta_mean_kcal_mol": float(np.mean(values)) if values else None,
                            "delta_sample_sd_kcal_mol": float(np.std(values, ddof=1)) if len(values) > 1 else None,
                            "mm_proposal_energy_error_rms_kj_mol": float(np.sqrt(np.mean([row["mm_integration_energy_error_kj_mol"]**2 for row in rows]))) if rows else None,
                            "qm_coupling_change_rms_kj_mol": float(np.sqrt(np.mean([(coupling * KJ_PER_KCAL * (row["proposal_delta_kcal_mol"] - row["previous_delta_kcal_mol"]))**2 for row in rows]))) if rows else None,
                            "mean_mm_only_acceptance_probability_for_same_proposals": float(np.mean(np.exp(np.minimum(0., -mm_errors / thermal_kj)))) if rows else None,
                            "mean_joint_acceptance_probability_for_same_proposals": float(np.mean(np.exp(np.minimum(0., -(mm_errors + quantum_changes) / thermal_kj)))) if rows else None,
                            "mean_proposal_rms_move_angstrom": float(np.mean([row["proposal_ligand_rms_move_angstrom"] for row in rows])) if rows else None})
        adjacent = [adjacent_diagnostic(production_by_lambda[first], production_by_lambda[second], first, second,
                                       protocol["temperature_kelvin"], protocol.get("delta_offset_kcal_mol", 0))
                    for first, second in zip(protocol["lambdas"], protocol["lambdas"][1:])]
        timing = [json.loads(path.read_text()) for path in directory.glob("bridge_*.energy.json")]
        record = {"status": protocol["status"], "windows": windows, "adjacent_diagnostics": adjacent,
                  "periodic_target_force_gate_passed": protocol.get("composed_periodic_mm_local_qmmm_force_validation_passed"),
                  "proposal_reversibility_max_angstrom": protocol.get("proposal_reversibility_max_angstrom"),
                  "qm_evaluations_completed": len(timing), "qm_wall_seconds_sum": sum(row["wall_seconds"] for row in timing),
                  "qm_cpu_seconds_sum": sum(row.get("qm_child_cpu_seconds", 0) for row in timing),
                  "mm_platform": protocol["platform"],
                  "mm_proposal_wall_seconds_sum_not_gpu_kernel_profiler_time": sum(row["md_wall_seconds"] for row in samples),
                  "proposal_integration_length_ps_not_physical_trajectory_time": len(samples) * protocol["md_steps_per_proposal"] * protocol["timestep_fs"] / 1000,
                  "lambda_one_reached": bool(samples) and max(row["lambda"] for row in samples) == 1,
                  "largest_ligand_displacement_from_start_angstrom": max((row["ligand_rms_displacement_from_start_angstrom"] for row in samples), default=None)}
        record["incomplete_jsonl_tail_waiting_for_writer"] = bool(lines and not lines[-1].endswith("\n"))
        if "error" in protocol:
            record["error"] = protocol["error"]
        summary["states"][state] = record
    statuses = [row["status"] for row in summary["states"].values()]
    if any(status == "ERROR" for status in statuses):
        summary["status"] = "ERROR"
    elif all(status == "PILOT_COMPLETE_NOT_CONVERGED_FREE_ENERGY" for status in statuses):
        summary["status"] = "PILOT_COMPLETE_NOT_CONVERGED_FREE_ENERGY"
    return summary


def markdown(summary):
    lines = ["# CDK8 QM/MM 桥接采样试运行", "", f"更新时间：{summary['generated_utc']}",
             f"状态：`{summary['status']}`。不是收敛自由能或全靶点排序结果。", "",
             "## 实际进度", "", "| 状态 | λ | 接受/尝试（含预热） | 生产样本/生产接受数 | Δδ 样本 SD (kcal/mol) |", "|---|---:|---:|---:|---:|"]
    for state, record in summary["states"].items():
        if not record.get("windows"):
            lines.append(f"| {state} | — | {record['status']} | — | — |")
        for window in record.get("windows", []):
            deviation = window["delta_sample_sd_kcal_mol"]
            deviation = f"{deviation:.4f}" if deviation is not None else "—"
            lines.append(f"| {state} | {window['lambda']:g} | {window['accepted']}/{window['attempted']} | {window['production_samples_including_rejections']}/{window['production_accepted_transitions']} | {deviation} |")
    coverage = "- 四态均已产生 λ=1 样本，但短链本身仍不足以证明平衡与收敛。" if all(row.get("lambda_one_reached") for row in summary["states"].values()) else "- 尚未把四态都连接至 λ=1，不能输出完整四态自由能修正。"
    lines += ["", "## 解释限制", "", coverage,
              "- 接受率、前后向重加权差和权重 ESS 仅作短试验诊断；单条短链没有证明平衡或独立采样。",
              "- MM 提案使用 GPU；QM/MM 接受率计算使用 CPU GAMESS。保存全部提案、拒绝状态、输入、日志、坐标与检查点。",
              "- 使用原 MM 参考、固定体积和固定壳层库存；没有修改现有 tPIE 评分或按实验标签调参。",
              "- 详细逐态计时、相邻窗口诊断和错误见同名 JSON；QM 子进程 CPU 秒与墙钟秒分别报告。"]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--wait", action="store_true")
    args = parser.parse_args()
    base = args.report.with_suffix("") if args.report.suffix in (".json", ".md") else args.report
    json_path, markdown_path = Path(str(base) + ".json"), Path(str(base) + ".md")
    deadline = time.monotonic() + 7200
    while True:
        summary = summarize(args.root)
        json_path.parent.mkdir(parents=True, exist_ok=True)
        with (args.root / ".bridge_report.lock").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            write_json(json_path, summary)
            temporary = Path(str(markdown_path) + ".tmp")
            temporary.write_text(markdown(summary))
            temporary.replace(markdown_path)
            fcntl.flock(lock, fcntl.LOCK_UN)
        if not args.wait or summary["status"] != "RUNNING" or time.monotonic() > deadline:
            break
        time.sleep(60)
    print(json.dumps({"status": summary["status"], "report": str(args.report)}, indent=2))


if __name__ == "__main__":
    main()
