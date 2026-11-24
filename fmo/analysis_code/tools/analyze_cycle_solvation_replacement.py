"""Replace, rather than double-add, pair SOLV with matched cycle solvation."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from fmotpiefinal.core.parser import parse_pie
from report_bace_extensions import paired_rank_summary
from report_local_sidechain_transfer import load_truth
from run_cdk8_solution_cycle import digest, metrics, write_json


def corrected_scores(tpie, pair_solv, cycle_electrostatic, cycle_total, strain=None):
    electronic = tpie-pair_solv
    scores = {"electronic_pair_plus_cycle_electrostatic": electronic+cycle_electrostatic,
              "electronic_pair_plus_cycle_total_solvation": electronic+cycle_total}
    if strain is not None:
        scores["electronic_pair_plus_cycle_total_solvation_and_strain"] = electronic+cycle_total+strain
    return scores


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evaluation", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    plan = {"schema": "pair-solv-replaced-by-cycle-solvation-v1", "driver_sha256": digest(__file__),
        "evaluation_sha256": digest(args.evaluation), "protocol_sha256": digest(args.protocol),
        "arms": ["electronic_pair_plus_cycle_electrostatic", "electronic_pair_plus_cycle_total_solvation",
                 "electronic_pair_plus_cycle_total_solvation_and_strain"],
        "equation": "electronic_pair = full tPIE - sum(pair SOLV); add matched A/B/AB PCM increment once, optionally add bound-minus-free-ligand-local-minimum energy at coefficient 1",
        "labels_used_for_fit": False, "production_default_changed": False,
        "scope": "Exploratory matched-geometry decomposition; pair term without SOLV is not itself full tPIE, no entropy or ensemble. After observed poor earlier controls, not a new blind test."}
    path = args.output.parent / (args.output.stem+".frozen.json")
    if path.exists() and json.loads(path.read_text()) != plan:
        raise ValueError("changed frozen diagnostic")
    write_json(path, plan)
    evaluation = json.loads(args.evaluation.read_text())
    frozen = json.loads(args.protocol.read_text())
    names = sorted(frozen["names"])
    if set(evaluation["scores"]) != set(names):
        raise ValueError("incomplete static three-state cohort")
    rows = {}
    fields = ("electrostatic", "exchange_repulsion", "charge_transfer", "dispersion", "solvation")
    for name in names:
        task = frozen["tasks"][frozen["case_tasks"][name]["complex"]]
        path = Path(task["input"])
        verified = evaluation["rows"][task["key"]]
        log = path.parent / "gamess.log"
        if digest(path) != verified["input_sha256"] or digest(log) != verified["log_sha256"]:
            raise ValueError("complex input/log changed")
        count = task["runtime"]["ligand_fragment"]
        pairs = [pair for pair in parse_pie(log) if count in (pair.fragment_i, pair.fragment_j)]
        if len(pairs) != count-1:
            raise ValueError("incomplete ligand pairs")
        if any(any(getattr(pair, field+"_kcal_mol") is None for field in fields) for pair in pairs):
            raise ValueError("missing PIEDA components")
        components = {field: sum(getattr(pair, field+"_kcal_mol") for pair in pairs) for field in fields}
        old = evaluation["scores"][name]
        tpie = sum(pair.value_kcal_mol for pair in pairs)
        if abs(tpie-old["matched_tpie"]) > 1.e-8:
            raise ValueError("pair score mismatch")
        if any(abs(sum(getattr(pair, field+"_kcal_mol") for field in fields)-pair.value_kcal_mol) > .006 for pair in pairs):
            raise ValueError("PIEDA components do not close within printed precision")
        inventory = old["pcm_audit"]["states"]
        delta = {field: inventory["complex"][field]-inventory["receptor"][field]-inventory["ligand"][field]
                 for field in ("electrostatic", "cavitation", "total")}
        rows[name] = {"components": components, "cycle_delta": delta, "matched_tpie": tpie,
            **corrected_scores(tpie, components["solvation"], delta["electrostatic"], delta["total"], old.get("ligand_strain"))}
    truth = load_truth(args.workspace.resolve(), frozen["target"])
    if set(truth) != set(names):
        raise ValueError("full evaluation identity mismatch")
    observed = [truth[name] for name in names]
    baseline = [-frozen["baseline"][name] for name in names]
    compared = {score: metrics(observed, [-rows[name][score] for name in names])
                for score in plan["arms"] if all(score in rows[name] for name in names)}
    result = {"status": "COMPLETE_MATCHED_STATIC_CYCLE", "n": len(names), "rows": rows,
        "metrics": compared, "baseline_metrics": metrics(observed, baseline),
        "paired_comparisons": {score: paired_rank_summary(observed, baseline, [-rows[name][score] for name in names])
                               for score in compared}, "plan": plan}
    write_json(args.output.with_suffix(".json"), result)
    lines = ["# 匹配三态溶剂修正：替换配对SOLV而不是重复相加", "",
        "完整队列，保留全部配体。以下是预声明单位系数诊断组合，不称为完整tPIE或结合自由能。", "",
        "| 方法 | Pearson | Spearman | Kendall |", "|---|---:|---:|---:|"]
    for score, value in [("baseline_full_tpie", result["baseline_metrics"]), *compared.items()]:
        lines.append(f"| {score} | {value['pearson']:.4f} | {value['spearman']:.4f} | {value['kendall']:.4f} |")
    lines += ["", "配对SOLV与整个结合过程的去溶剂化不是同一个量。该消融先从完整tPIE扣除配对SOLV，再加入A/B/AB一致PCM计算得到的溶剂增量一次。",
        "这仍省略受体/配体内部重组、溶液构象人口与熵；诊断组合不是严格自由能循环，不因指标较高就自动纳入默认。",
        "参数未拟合，但数据集已经反复用于开发；需跨集合回归与独立验证。完整bootstrap及逐例分解见JSON。"]
    args.output.with_suffix(".md").write_text("\n".join(lines)+"\n")


if __name__ == "__main__":
    main()
