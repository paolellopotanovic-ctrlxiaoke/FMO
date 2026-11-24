"""Complete-cohort local-strain metrics, including fixed OpenBind aggregation."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from fmotpiefinal.scoring import FACTOR, parameters
from report_local_sidechain_transfer import load_truth
from run_cdk8_solution_cycle import digest, metrics, write_json
from run_local_sidechain_transfer import directory_key
from run_local_strain_transfer import verify_state


ARMS = ("baseline_tpie", "matched_tpie", "matched_tpie_plus_ligand_strain")


def aggregate_crystals(compounds, cases):
    used, output = set(), {}
    for compound in compounds:
        identifiers = [name.strip() for name in compound["fragalysis_codes"].split(";")]
        if not identifiers or len(set(identifiers)) != len(identifiers) or any(name in used for name in identifiers):
            raise ValueError("duplicate or empty OpenBind crystal mapping")
        if any(name not in cases for name in identifiers):
            raise ValueError("missing OpenBind crystal; no best-pose or partial aggregation")
        used.update(identifiers)
        key = str(compound["compound_index"])
        if key in output:
            raise ValueError("duplicate compound identity")
        output[key] = {arm: float(np.mean([cases[name][arm] for name in identifiers])) for arm in ARMS}
        output[key].update(n_crystals=len(identifiers), experimental_score=float(compound["experimental_pKD"]))
    if used != set(cases):
        raise ValueError("OpenBind crystal inventory mismatch")
    return output


def score_metrics(truth, records):
    if set(truth) != set(records):
        raise ValueError("matched full cohort required")
    names = sorted(records)
    observed = [truth[name] for name in names]
    output = {}
    for arm in ARMS:
        scores = [-records[name][arm] for name in names]
        raw = metrics(observed, scores)
        settings = parameters(scores)
        soft = None
        if settings is not None:
            width = FACTOR * settings["scaled_mad"]
            transformed = [settings["median"] + width * float(np.tanh((value - settings["median"]) / width)) for value in scores]
            if len(set(transformed)) == len(set(scores)):
                soft = metrics(observed, transformed)
                if any(abs(soft[metric] - raw[metric]) > 1.e-10 for metric in ("spearman", "kendall")):
                    raise ValueError("soft compression changed rank ordering")
        output[arm] = {"raw": raw, "soft": soft, "soft_reference": settings,
                       "soft_scope": "Same fixed tanh/3MAD formula, label-free cohort-local scale for each readout; not the bundled pure-tPIE reference or a fitted calibration"}
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", required=True, type=Path)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--verify-logs", action="store_true")
    args = parser.parse_args()
    frozen_path = args.run / "frozen_protocol.json"
    frozen = json.loads(frozen_path.read_text())
    groups = {}
    for target, expected in frozen["targets"].items():
        rows, failures = {}, {}
        for key, case in frozen["cases"].items():
            if case["target"] != target:
                continue
            path = args.run / target / directory_key(case["name"]) / "case_result.json"
            if not path.exists():
                continue
            row = json.loads(path.read_text())
            if row["status"] != "SUCCEEDED_INDEPENDENTLY_VERIFIED":
                failures[case["name"]] = row
                continue
            if args.verify_logs:
                for state in row["states"].values():
                    verify_state(state["task"], state["row"], case, str(args.workspace.resolve()))
            bound = row["states"]["ligand"]["row"]["energy_hartree"]
            free = row["states"]["free_ligand_minimum"]["row"]["energy_hartree"]
            strain = (bound - free) * 627.5094740631
            tpie = row["states"]["complex"]["row"]["tpie"]
            if abs(strain - row["local_ligand_strain"]) > 1.e-8 or abs(tpie + strain - row["matched_tpie_plus_ligand_strain"]) > 1.e-8:
                raise ValueError("reported score does not close to accepted state energies")
            rows[case["name"]] = row
        group = {"n_expected_structures": expected, "n_succeeded_structures": len(rows), "failures": failures,
                 "status": "COMPLETE" if len(rows) == expected else "INCOMPLETE", "metrics": {}}
        groups[target] = group
        if len(rows) != expected:
            continue
        if target == "OpenBind":
            mapping_path = args.workspace / "runs/openbind/score_distribution_audit/paired_compounds.csv"
            preflight = json.loads((args.run / "openbind_aggregation_preflight.json").read_text())
            if preflight["status"] != "PASSED" or digest(mapping_path) != preflight["mapping_sha256"]:
                raise ValueError("OpenBind frozen crystal-to-compound mapping changed")
            compounds = list(csv.DictReader(mapping_path.open()))
            if len(compounds) != 490:
                raise ValueError("all 490 compounds required")
            records = aggregate_crystals(compounds, rows)
            truth = {name: row["experimental_score"] for name, row in records.items()}
            group.update(mapping_sha256=digest(mapping_path), n_compounds=len(records), aggregation=frozen["openbind_aggregation"])
        else:
            truth = load_truth(args.workspace.resolve(), target)
            records = rows
        group["metrics"] = score_metrics(truth, records)
        group["scores"] = {name: {field: row[field] for field in ARMS} for name, row in records.items()}
    report = {"status": "COMPLETE" if all(group["status"] == "COMPLETE" for group in groups.values()) else "RUNNING_OR_INCOMPLETE",
              "groups": groups, "protocol_sha256": digest(frozen_path), "driver_sha256": digest(__file__),
              "logs_reverified_in_this_report": args.verify_logs, "production_default_changed": False,
              "warning": "Same frozen unit coefficients on all cohorts; no per-target best-arm selection, no new hydrogen normalization, no case deletion; retrospective transfer, not binding free energy"}
    original = ("CDK2", "JNK1", "P38", "TYK2")
    if all(groups[target]["status"] == "COMPLETE" for target in original):
        report["original87_case_weighted_target_metrics"] = {
            arm: {metric: sum(groups[target]["n_expected_structures"] * groups[target]["metrics"][arm]["raw"][metric] for target in original) / 87.
                  for metric in ("pearson", "spearman", "kendall")} for arm in ARMS}
    write_json(args.output.with_suffix(".json"), report)
    lines = ["# 局部配体应变：全部测试集固定协议对照", "", "主候选：匹配IDISP0全AB tPIE + [PCM游离配体结合几何能量 − PCM游离局部优化能量]，系数均为1。",
             "不加入额外H归一化、NH2重建、分子量或按target选择的权重；不是完整结合自由能。",
             "所有行使用完整集合；未算齐不公布子集指标。复用已验收制备/k10快照，不是原始PDB全链重跑。", "",
             "| 集合 | 结构覆盖 | 原生产纯tPIE P/S/K | 匹配纯tPIE P/S/K | 匹配tPIE＋应变 P/S/K | 匹配纯tPIE软P | 加应变后软P |", "|---|---:|---|---|---|---:|---:|"]
    for target, group in groups.items():
        values = group["metrics"]
        formatted = [" / ".join(f"{values[arm]['raw'][metric]:.3f}" for metric in ("pearson", "spearman", "kendall")) if values else "待完整验收" for arm in ARMS]
        soft = values[ARMS[-1]]["soft"] if values else None
        pure_soft = values["matched_tpie"]["soft"] if values else None
        pure_soft_text = f"{pure_soft['pearson']:.3f}" if pure_soft else "—"
        soft_text = f"{soft['pearson']:.3f}" if soft else "—"
        lines.append(f"| {target} | {group['n_succeeded_structures']}/{group['n_expected_structures']} | {' | '.join(formatted)} | {pure_soft_text} | {soft_text} |")
    lines.extend(["", "Pearson主列是未软压缩值。软压缩使用相同固定tanh/3MAD形式，但为各读数计算自身不使用标签的集合尺度，不能冒充原纯tPIE内置参照。",
                  "单调软压缩不改变Spearman/Kendall；若浮点饱和导致新并列则不输出该软压缩指标。",
                  "OpenBind固定将全部639晶体的能量分数先按化合物取算术平均，再在490化合物上计算指标；不择优晶体。",
                  "原87例如汇总，采用四个target指标按例数加权，不将不同target绝对能量混为一个池相关。", ""])
    args.output.with_suffix(".md").write_text("\n".join(lines))
    print(report["status"], {target: group["n_succeeded_structures"] for target, group in groups.items()}, flush=True)


if __name__ == "__main__":
    main()
