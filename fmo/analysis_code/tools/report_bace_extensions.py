"""Independently audit full-cohort BACE controls and quantify paired rank changes."""
from __future__ import annotations

import argparse
from itertools import combinations
import io
import json
from pathlib import Path

import numpy as np
from scipy.stats import kendalltau, pearsonr, spearmanr

from fmotpiefinal.core.fmo_convergence import final_fmo_convergence
from fmotpiefinal.core.parser import parse_pie
from fmotpiefinal.scoring import FACTOR, choose_reference
from report_local_sidechain_transfer import load_truth
from run_cdk8_solution_cycle import digest, metrics, write_json


def paired_rank_summary(observed, baseline, candidate, n_bootstrap=4000):
    observed, baseline, candidate = map(np.asarray, (observed, baseline, candidate))
    if not (len(observed) == len(baseline) == len(candidate)) or len(observed) < 4:
        raise ValueError("matched complete arrays required")
    counts = {"corrected": 0, "regressed": 0, "unchanged": 0, "tied": 0}
    for first, second in combinations(range(len(observed)), 2):
        truth = np.sign(observed[first]-observed[second])
        before = truth*np.sign(baseline[first]-baseline[second])
        after = truth*np.sign(candidate[first]-candidate[second])
        if truth == 0 or before == 0 or after == 0:
            counts["tied"] += 1
        elif before == after:
            counts["unchanged"] += 1
        elif after > before:
            counts["corrected"] += 1
        else:
            counts["regressed"] += 1
    rng = np.random.default_rng(20261024)
    deltas = []
    for _ in range(n_bootstrap):
        indices = rng.integers(len(observed), size=len(observed))
        sampled = [values[indices] for values in (observed, baseline, candidate)]
        if any(np.ptp(values) == 0 for values in sampled):
            continue
        deltas.append([function(sampled[0], sampled[2]).statistic-
                       function(sampled[0], sampled[1]).statistic
                       for function in (pearsonr, spearmanr, kendalltau)])
    if not deltas:
        raise ValueError("no valid paired bootstrap samples")
    leave_one_out = []
    for excluded in range(len(observed)):
        keep = np.arange(len(observed)) != excluded
        leave_one_out.append(float(spearmanr(observed[keep], candidate[keep]).statistic-
                                   spearmanr(observed[keep], baseline[keep]).statistic))
    return {"pair_changes": counts, "paired_bootstrap_95_percent_delta": {
        name: list(map(float, np.quantile(np.asarray(deltas)[:, index], [.025, .975])))
        for index, name in enumerate(("pearson", "spearman", "kendall"))},
        "bootstrap_valid_samples": len(deltas), "seed": 20261024,
        "leave_one_out_spearman_delta_range": [min(leave_one_out), max(leave_one_out)],
        "warning": "Development-set paired bootstrap; does not correct repeated method exploration or prove transfer"}


def audit(root, names):
    evaluation = json.loads((root / "evaluation.json").read_text())
    if evaluation["status"] != "COMPLETE_VERIFIED" or set(evaluation["rows"]) != set(names):
        raise ValueError("incomplete control; partial-cohort metrics prohibited")
    frozen = json.loads((root / "frozen_protocol.json").read_text())
    target = frozen.get("target", frozen.get("cohort", {}).get("target"))
    proof = {}
    for name in names:
        parent = root / target / name
        row = evaluation["rows"][name]
        deck, log = parent / "input.inp", parent / "run/gamess.log"
        runtime = json.loads((parent / "run/result.json").read_text())
        if row["status"] != "SUCCEEDED" or digest(deck) != row["input_sha256"] or digest(log) != row["log_sha256"]:
            raise ValueError("saved evaluation fingerprint mismatch")
        executable = Path(runtime["gamess_dir"]) / "gamess.00.x"
        if runtime["executable_sha256"] != digest(executable):
            raise ValueError("QM executable changed since execution")
        convergence = final_fmo_convergence(deck.read_text(), io.StringIO(log.read_text()))
        count = runtime["ligand_fragment"]
        pairs = [pair for pair in parse_pie(log) if count in (pair.fragment_i, pair.fragment_j)]
        if len(pairs) != count-1 or abs(sum(pair.value_kcal_mol for pair in pairs)-row["tpie"]) > 1.e-8:
            raise ValueError("independent pair inventory/sum check failed")
        proof[name] = {"input_sha256": digest(deck), "log_sha256": digest(log),
                       "executable_sha256": runtime["executable_sha256"], "convergence": convergence}
    return evaluation, proof


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--controls", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cohort = json.loads(args.source.read_text())["cohorts"]["BACE"]
    names = sorted(cohort["names"])
    truth = load_truth(args.workspace.resolve(), "BACE")
    if set(truth) != set(names):
        raise ValueError("full cohort identity mismatch")
    observed = [truth[name] for name in names]
    baseline = [-cohort["baseline"][name]["tpie"] for name in names]
    settings, _ = choose_reference("BACE", baseline, None)
    width = FACTOR*settings["scaled_mad"]
    soft = lambda values: settings["median"]+width*np.tanh((np.asarray(values)-settings["median"])/width)
    groups = {}
    for root in args.controls:
        evaluation, proof = audit(root.resolve(), names)
        candidate = [-evaluation["rows"][name]["tpie"] for name in names]
        groups[root.name] = {"metrics": metrics(observed, candidate),
            "soft_metrics": metrics(observed, soft(candidate)),
            "paired_comparison": paired_rank_summary(observed, baseline, candidate),
            "n": len(names), "independent_evidence": proof,
            "wall_seconds": evaluation["wall_seconds"]}
    result = {"status": "COMPLETE_INDEPENDENTLY_VERIFIED", "baseline": metrics(observed, baseline),
        "baseline_soft": metrics(observed, soft(baseline)), "controls": groups,
        "production_default_changed": False, "cohort": "BACE36", "no_partial_cohort_metrics": True}
    write_json(args.output.with_suffix(".json"), result)
    lines = ["# BACE 完整36例：口袋与输入几何对照", "",
        "每个对照重新核对全部QM输入/日志/可执行文件哈希、严格最终收敛和完整PIEDA之和。",
        "软Pearson沿用同一冻结基线尺度；没有按实验拟合参数。", "",
        "| 方法 | 覆盖 | Pearson | 软Pearson | Spearman | Kendall |",
        "|---|---:|---:|---:|---:|---:|"]
    for label, raw, capped in [("当前主流程", result["baseline"], result["baseline_soft"])]+[
            (name, values["metrics"], values["soft_metrics"]) for name, values in groups.items()]:
        lines.append(f"| {label} | 36/36 | {raw['pearson']:.4f} | {capped['pearson']:.4f} | {raw['spearman']:.4f} | {raw['kendall']:.4f} |")
    lines += ["", "## 变化稳健性", ""]
    for name, values in groups.items():
        comparison = values["paired_comparison"]
        lines += [f"### {name}", "", f"- 配体对变化：{comparison['pair_changes']}",
            f"- 配对bootstrap的Spearman变化95%区间：{comparison['paired_bootstrap_95_percent_delta']['spearman']}",
            f"- 去任意一例后的Spearman变化范围：{comparison['leave_one_out_spearman_delta_range']}", ""]
    lines += ["这里的区间不能校正反复研究同一开发集的选择偏差，也不代表跨集合提升。",
              "无论结果高低，未经结构归因与Galectin3等回归，不自动更改生产流程。"]
    args.output.with_suffix(".md").write_text("\n".join(lines)+"\n")


if __name__ == "__main__":
    main()
