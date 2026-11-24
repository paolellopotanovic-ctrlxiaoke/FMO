"""Report graph-guided conditional execution, including unchanged baseline reuse."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from fmotpiefinal.core.structure import read_pdb
from fmotpiefinal.pipeline.target_stages import receptor_at_minimized_geometry
from fmotpiefinal.scoring import choose_reference, FACTOR
from report_local_sidechain_transfer import load_truth
from run_cdk8_solution_cycle import digest, metrics, write_json
from run_local_sidechain_transfer import directory_key


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", required=True, type=Path)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    workspace, root = args.workspace.resolve(), args.run.resolve()
    plan = json.loads((root / "frozen_protocol.json").read_text())
    final = json.loads((root / "final_status.json").read_text())
    if final["evaluation_sha256"] != digest(root / "evaluation.json"):
        raise ValueError("final evaluation fingerprint mismatch")
    progress = json.loads((root / "progress.json").read_text())
    groups = {}
    lines = ["# 图驱动选择性侧链处理：完整验证结果", "",
        "固定规则，不依据实验标签选残基；OpenBind 不在本轮范围。", "",
        "| 集合 | 完成 | 选择侧链数 | 新QM例数 | 基线 P / 软P / S / K | 图驱动 P / 软P / S / K |",
        "|---|---:|---:|---:|---|---|"]
    for target, cohort in plan["cohorts"].items():
        names = cohort["names"]
        rows = progress[target]
        failures = {name: rows.get(name, {"status": "MISSING"}) for name in names if rows.get(name, {}).get("status") != "SUCCEEDED"}
        group = {"n_expected": len(names), "n_succeeded": len(names)-len(failures), "failures": failures,
            "selected_residues": cohort["local_residues"],
            "n_new_qm": sum(bool(row.get("new_qm_executed")) for row in rows.values()),
            "n_unchanged_reused": sum(row.get("execution_path") == "UNCHANGED_VERIFIED_BASELINE_NO_STRUCTURAL_TRIGGER" for row in rows.values())}
        if not failures:
            truth = load_truth(workspace, target)
            if set(truth) != set(names):
                raise ValueError("label identity mismatch")
            observed = [truth[name] for name in names]
            baseline = [-cohort["baseline"][name]["tpie"] for name in names]
            candidate = [-rows[name]["tpie"] for name in names]
            settings, reference = choose_reference(target, baseline, None)
            width = FACTOR * settings["scaled_mad"]
            transform = lambda values: [settings["median"]+width*math.tanh((value-settings["median"])/width) for value in values]
            group.update(baseline=metrics(observed, baseline), graph_guided=metrics(observed, candidate),
                baseline_soft=metrics(observed, transform(baseline)), graph_soft=metrics(observed, transform(candidate)),
                softcap_reference=reference,
                rows=[{"ligand": name, "experimental_minus_dG": truth[name],
                    "baseline_tpie": cohort["baseline"][name]["tpie"], "candidate_tpie": rows[name]["tpie"],
                    "execution_path": rows[name]["execution_path"]} for name in names])
            if cohort["local_residues"]:
                protein = read_pdb(cohort["protein"])
                selected = {tuple(residue) for residue in cohort["local_residues"]}
                backbone = {"N", "CA", "C", "O", "OXT", "H", "H1", "H2", "H3", "HA", "HA2", "HA3"}
                fixed = [index for index, atom in enumerate(protein.atoms)
                    if (atom.chain, atom.residue_number, atom.residue_name) not in selected or atom.atom_name in backbone]
                audit = {}
                for name in names:
                    directory = root / target / directory_key(name)
                    receptor = receptor_at_minimized_geometry(protein, directory / "minimization/minimized_complex.pdb")
                    displacement = np.linalg.norm(np.asarray([atom.xyz for atom in receptor.atoms])-
                                                  np.asarray([atom.xyz for atom in protein.atoms]), axis=1)
                    maximum = float(max(displacement[fixed]))
                    if maximum > 0.002:
                        raise ValueError("fixed backbone/water/remote atom moved")
                    audit[name] = {"max_fixed_displacement_A": maximum, "max_receptor_displacement_A": float(max(displacement))}
                group["independent_geometry_audit"] = audit
            elif candidate != baseline:
                raise ValueError("untriggered target changed scores")
            display = lambda raw, soft: f"{group[raw]['pearson']:.3f} / {group[soft]['pearson']:.3f} / {group[raw]['spearman']:.3f} / {group[raw]['kendall']:.3f}"
            lines.append(f"| {target} | {len(names)}/{len(names)} | {len(cohort['local_residues'])} | {group['n_new_qm']} | {display('baseline', 'baseline_soft')} | {display('graph_guided', 'graph_soft')} |")
        else:
            lines.append(f"| {target} | {group['n_succeeded']}/{len(names)} | {len(cohort['local_residues'])} | {group['n_new_qm']} | — | 不报告成功子集指标 |")
        groups[target] = group
    lines += ["", "## 结论边界", "",
        "未触发集合复用已验收基线，并未重新计算 QM；相同指标是保留原路径的结果，不是新提分证据。",
        "触发集合只执行冻结的局部选择，不按结果扩大邻域或换阈值。",
        "结构规则有针对性不代表排序一定提高；只有本轮完整结果能说明该候选的实际效果。",
        "本轮为回顾性方法开发，不是外部前瞻盲测。生产默认未改变。"]
    if "CDK8" in groups and not groups["CDK8"]["failures"]:
        cdk8 = groups["CDK8"]
        lines += ["", "## 当前发现", "",
            f"CDK8 本轮 Spearman {cdk8['graph_guided']['spearman']:.3f}，基线 {cdk8['baseline']['spearman']:.3f}。",
            "本轮选择：" + ", ".join(f"{residue[0]}:{residue[2]}:{residue[1]}" for residue in cdk8["selected_residues"]) + "。",
            "能检测局部几何异常不等于已找全需要协同变化的残基；需要与固定受体和宽范围柔性对照共同比较。",
            "未触发集合不代表没有问题；只能说明本轮冻结的结构证据没有覆盖其问题。",
            "后续侧链构象和局部力学响应对照须单独预注册，不能看指标逐残基加到满意为止。"]
    result = {"schema": "graph-sidechain-transfer-result-v1", "final_status": final, "groups": groups,
        "frozen_protocol_sha256": digest(root / "frozen_protocol.json"), "production_default_changed": False,
        "total_new_qm": sum(group["n_new_qm"] for group in groups.values()),
        "total_unchanged_reused": sum(group["n_unchanged_reused"] for group in groups.values())}
    write_json(args.output.with_suffix(".json"), result)
    args.output.write_text("\n".join(lines)+"\n")
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
