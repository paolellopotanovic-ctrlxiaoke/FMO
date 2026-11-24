"""Decompose the fixed PCM local-strain correction without searching weights."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from fmotpiefinal.core.parser import parse_pie
from fmotpiefinal.core.solution_cycle import HARTREE_TO_KCAL_MOL, pcm_energy_inventory
from report_local_sidechain_transfer import load_truth
from run_cdk8_solution_cycle import digest, metrics, write_json


PAIR_FIELDS = ("electrostatic", "exchange_repulsion", "charge_transfer", "dispersion", "solvation")


def strain_components(bound, free, expected):
    internal = (bound["internal_hartree"]-free["internal_hartree"])*HARTREE_TO_KCAL_MOL
    solvent = bound["total"]-free["total"]
    if abs(internal+solvent-expected) > .01:
        raise ValueError("PCM strain components do not close within printed precision")
    return {"internal_response": internal, "solvent_response": solvent,
            "pcm_local_strain": expected, "closure_error": internal+solvent-expected}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--evaluation", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    evaluation = json.loads(args.evaluation.read_text())
    frozen = json.loads(args.protocol.read_text())
    cohort = json.loads(args.registry.read_text())["cohorts"][frozen["target"]]
    names = sorted(frozen["names"])
    if (evaluation["status"] != "COMPLETE_VERIFIED" or set(evaluation["scores"]) != set(names)
            or set(names) != set(cohort["names"])):
        raise ValueError("complete matched cohort required")
    plan = {"schema": "pcm-local-strain-components-v1", "target": frozen["target"],
        "weights": "All coefficients fixed to one; no label-dependent arm selection",
        "arms": ["baseline_tpie", "matched_tpie", "baseline_tpie_plus_pcm_strain",
                 "matched_tpie_plus_pcm_strain", "matched_tpie_plus_internal_response",
                 "matched_tpie_plus_solvent_response"],
        "evaluation_sha256": digest(args.evaluation), "protocol_sha256": digest(args.protocol),
        "registry_sha256": digest(args.registry), "driver_sha256": digest(__file__),
        "development_warning": "Retrospective diagnostic after complete strain scores were observed; not blind validation"}
    plan_path = args.output.with_suffix(".frozen.json")
    if plan_path.exists() and json.loads(plan_path.read_text()) != plan:
        raise ValueError("refusing to replace frozen diagnostic")
    write_json(plan_path, plan)
    rows = {}
    fields = PAIR_FIELDS
    for name in names:
        score = evaluation["scores"][name]
        inventories = {}
        for state in ("ligand", "free_ligand_minimum"):
            key = frozen["case_tasks"][name][state]
            task, verified = frozen["tasks"][key], evaluation["rows"][key]
            path = Path(task["input"])
            if digest(path) != verified["input_sha256"] or digest(path.parent/"gamess.log") != verified["log_sha256"]:
                raise ValueError("input/log fingerprint mismatch")
            inventories[state] = pcm_energy_inventory((path.parent/"gamess.log").read_text())
        decomposition = strain_components(inventories["ligand"], inventories["free_ligand_minimum"], score["ligand_strain"])
        complex_key = frozen["case_tasks"][name]["complex"]
        complex_path = Path(frozen["tasks"][complex_key]["input"])
        old_input = Path(cohort["sources"][name])
        old_log = old_input.parent/"run/gamess.log"
        pair_sets = {}
        for label, log in (("baseline", old_log), ("matched", complex_path.parent/"gamess.log")):
            runtime = cohort["runtime"][name]
            pairs = [pair for pair in parse_pie(log) if runtime["ligand_fragment"] in (pair.fragment_i, pair.fragment_j)]
            if len(pairs) != runtime["ligand_fragment"]-1:
                raise ValueError("incomplete PIEDA inventory")
            if any(getattr(pair, field+"_kcal_mol") is None for pair in pairs for field in fields):
                raise ValueError("missing PIEDA component")
            pair_sets[label] = {field: sum(getattr(pair, field+"_kcal_mol") for pair in pairs) for field in fields}
        rows[name] = {"decomposition": decomposition,
            "pair_component_delta": {field: pair_sets["matched"][field]-pair_sets["baseline"][field] for field in fields},
            "scores": {"baseline_tpie": score["baseline_tpie"], "matched_tpie": score["matched_tpie"],
                "baseline_tpie_plus_pcm_strain": score["baseline_tpie"]+score["ligand_strain"],
                "matched_tpie_plus_pcm_strain": score["matched_tpie"]+score["ligand_strain"],
                "matched_tpie_plus_internal_response": score["matched_tpie"]+decomposition["internal_response"],
                "matched_tpie_plus_solvent_response": score["matched_tpie"]+decomposition["solvent_response"]}}
    truth = load_truth(args.workspace.resolve(), frozen["target"])
    if set(truth) != set(names):
        raise ValueError("complete evaluation truth inventory required")
    results = {arm: metrics([truth[name] for name in names], [-rows[name]["scores"][arm] for name in names]) for arm in plan["arms"]}
    report = {"plan": plan, "n": len(names), "rows": rows, "metrics": results,
        "warning": "Internal response is evaluated in each PCM-polarized state, not a gas-phase strain; these are not free energies"}
    write_json(args.output.with_suffix(".json"), report)
    lines = [f"# {frozen['target']}：配体局部应变拆解", "", "| 固定系数组合 | Pearson | Spearman | Kendall |", "|---|---:|---:|---:|"]
    for arm, values in results.items():
        lines.append(f"| {arm} | {values['pearson']:.4f} | {values['spearman']:.4f} | {values['kendall']:.4f} |")
    lines += ["", "内部能响应含各自PCM极化，不等于真空形变能。溶剂响应为同一配体两种几何的PCM总项之差，不是整个结合过程的去溶剂化。",
        "本分析在看到完整应变组合后设计，属于开发集诊断；不依据每个集合的指标挑选不同方法，不自动更改生产默认。"]
    args.output.with_suffix(".md").write_text("\n".join(lines)+"\n")


if __name__ == "__main__":
    main()
