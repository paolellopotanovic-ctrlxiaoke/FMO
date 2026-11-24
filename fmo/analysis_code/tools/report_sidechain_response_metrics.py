"""Evaluate frozen structural decisions, never choose branches by affinity."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from run_cdk8_solution_cycle import digest, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--decisions", required=True, type=Path)
    parser.add_argument("--verification", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    decisions = json.loads(args.decisions.read_text())
    verification = json.loads(args.verification.read_text())
    expected = decisions["source_sha256"].get(str(args.verification.resolve()))
    if expected != digest(args.verification):
        raise ValueError("decision verification fingerprint mismatch")
    rows = {}
    lines = ["# 冻结结构响应规则：开发集合回归", "",
        "先保存结构决定再评估；不按每个集合指标择优。这些是既有开发集合，不是独立盲测。", "",
        "| 集合 | 决定 | 完成 | 基线 软P / S / K | 规则输出 软P / S / K |",
        "|---|---|---:|---|---|"]
    for target, group in decisions["groups"].items():
        decision = group["decision"]
        accepted = decision["action"] == "ACCEPT_LOCAL_REPAIR_RESEARCH_CONTROL"
        source = verification["groups"][target]
        if accepted:
            if source["failures"] or "graph_guided" not in source:
                raise ValueError("cannot accept incomplete candidate")
            baseline = source["baseline"]
            baseline_soft = source["baseline_soft"]
            candidate, soft = source["graph_guided"], source["graph_soft"]
        else:
            baseline = source.get("baseline")
            baseline_soft = source.get("baseline_soft")
            if baseline is None:
                protocol_path = next(Path(path) for path in decisions["source_sha256"] if Path(path).name == "frozen_protocol.json")
                cohort = json.loads(protocol_path.read_text())["cohorts"][target]
                from report_local_sidechain_transfer import load_truth
                from run_cdk8_solution_cycle import metrics
                from fmotpiefinal.scoring import choose_reference, FACTOR
                import math
                original_workspace = Path(cohort["protein"])
                workspace = next(parent for parent in original_workspace.parents if (parent / "data/fmo87_cohort.csv").exists())
                truth = load_truth(workspace, target)
                names = cohort["names"]
                observed = [truth[name] for name in names]
                scores = [-cohort["baseline"][name]["tpie"] for name in names]
                settings, _ = choose_reference(target, scores, None)
                width = FACTOR * settings["scaled_mad"]
                baseline = metrics(observed, scores)
                baseline_soft = metrics(observed, [settings["median"]+width*math.tanh((score-settings["median"])/width) for score in scores])
            candidate, soft = baseline, baseline_soft
        rows[target] = {"action": decision["action"], "baseline": baseline, "baseline_soft": baseline_soft,
            "output": candidate, "output_soft": soft,
            "execution_path": "REUSE_COMPLETED_GRAPH_CANDIDATE" if accepted else "REUSE_VERIFIED_BASELINE",
            "n": len(group["cases"])}
        display = lambda raw, transformed: f"{transformed['pearson']:.3f} / {raw['spearman']:.3f} / {raw['kendall']:.3f}"
        lines.append(f"| {target} | {'有限机制修复' if accepted else '保留基线'} | {len(group['cases'])}/{len(group['cases'])} | {display(baseline, baseline_soft)} | {display(candidate, soft)} |")
    lines += ["", "本报告只重用已有计算，并未再次执行 270 次 QM。FXa 宽侧链候选失败，规则输出保留独立验收基线，不将其冒充候选修复成功。",
        "仅 CDK8 满足本轮机制，不代表该规则能修复 CHK1/BACE 的低指标；它们仍保留原分数。",
        "规则是已知 CDK8 现象启发的回顾性开发；不能用保留其它集合不变声称普适提分已解决。",
        "生产默认未改变，仍需新集合的前瞻验证。"]
    args.output.write_text("\n".join(lines)+"\n")
    write_json(args.output.with_suffix(".json"), {"groups": rows, "decision_sha256": digest(args.decisions),
        "verification_sha256": digest(args.verification), "new_qm_in_this_report": 0,
        "production_default_changed": False})
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
