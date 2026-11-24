"""Independent matched pose-control coverage, geometry and scoring report."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from rdkit import Chem

from fmotpiefinal.core.bound_geometry import load_bound_ligand
from fmotpiefinal.core.structure import read_pdb
from fmotpiefinal.pipeline.target_stages import receptor_at_minimized_geometry
from fmotpiefinal.scoring import FACTOR, choose_reference
from report_local_sidechain_transfer import load_truth
from report_weak_mm_diagnostics import angle_changes
from run_cdk8_solution_cycle import digest, metrics, write_json
from run_local_sidechain_transfer import directory_key


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.run.resolve()
    plan = json.loads((root / "frozen_protocol.json").read_text())
    final = json.loads((root / "final_status.json").read_text())
    if final["evaluation_sha256"] != digest(root / "evaluation.json"):
        raise ValueError("evaluation fingerprint mismatch")
    progress = json.loads((root / "progress.json").read_text())
    control = plan["control"] if plan["schema"] == "matched-ligand-geometry-control-v2" else "input_heavy_fixed"
    lines = [f"# 配体几何对照：{control}", "",
        "受体全部固定；保持原 QM、口袋、分片和电荷。软压缩使用同一基线参考尺度，不重新拟合。", "",
        "| 集合 | 完成 | 基线 P / 软P / S / K | 对照 P / 软P / S / K |",
        "|---|---:|---|---|"]
    groups = {}
    for target, cohort in plan["cohorts"].items():
        names = cohort["names"]
        rows = progress[target]
        if set(rows) != set(names):
            raise ValueError("full result inventory required")
        molecules = {molecule.GetProp("_Name").strip(): molecule for molecule in Chem.SDMolSupplier(cohort["sdf"], removeHs=False) if molecule is not None}
        if set(molecules) != set(names):
            raise ValueError("full molecule inventory required")
        failures = {name: row for name, row in rows.items() if row["status"] != "SUCCEEDED"}
        protein = read_pdb(cohort["protein"])
        receptor_positions = np.asarray([atom.xyz for atom in protein.atoms])
        geometry = {}
        for name in names:
            path = root / target / directory_key(name) / "minimization/minimized_complex.pdb"
            if not path.exists():
                geometry[name] = {"status": "MISSING_GEOMETRY"}
                continue
            receptor = receptor_at_minimized_geometry(protein, path)
            difference = np.asarray([atom.xyz for atom in receptor.atoms])-receptor_positions
            fixed_max = float(np.max(np.linalg.norm(difference, axis=1)))
            if fixed_max > 0.002:
                raise ValueError("fixed receptor coordinate changed")
            ligand = load_bound_ligand(cohort["sdf"], name, path)
            molecule = molecules[name]
            positions = np.asarray([atom.xyz for atom in ligand.atoms])
            original = molecule.GetConformer().GetPositions()
            heavy = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1]
            heavy_max = float(np.max(np.linalg.norm(positions[heavy]-original[heavy], axis=1)))
            if control == "input_heavy_fixed" and heavy_max > 0.001:
                raise ValueError("input-heavy-fixed ligand moved")
            baseline = np.asarray([atom.xyz for atom in load_bound_ligand(cohort["sdf"], name, cohort["complexes"][name]).atoms])
            geometry[name] = {"status": "VERIFIED", "fixed_receptor_max_displacement_A": fixed_max,
                "ligand_heavy_max_displacement_A": heavy_max,
                "largest_bond_angle_changes_from_baseline": angle_changes(molecule, baseline, positions)[:5],
                "complex_sha256": digest(path)}
        group = {"n": len(names), "n_succeeded": len(names)-len(failures), "failures": failures,
            "independent_geometry": geometry}
        truth = load_truth(args.workspace.resolve(), target)
        if set(truth) != set(names):
            raise ValueError("full label inventory required")
        observed = [truth[name] for name in names]
        baseline_scores = [-cohort["baseline"][name]["tpie"] for name in names]
        settings, reference = choose_reference(target, baseline_scores, None)
        width = FACTOR * settings["scaled_mad"]
        transform = lambda scores: [settings["median"]+width*math.tanh((score-settings["median"])/width) for score in scores]
        group.update(baseline=metrics(observed, baseline_scores), baseline_soft=metrics(observed, transform(baseline_scores)), soft_reference=reference)
        display = lambda raw, soft: f"{raw['pearson']:.3f} / {soft['pearson']:.3f} / {raw['spearman']:.3f} / {raw['kendall']:.3f}"
        if not failures:
            scores = [-rows[name]["tpie"] for name in names]
            group.update(candidate=metrics(observed, scores), candidate_soft=metrics(observed, transform(scores)))
            candidate_text = display(group["candidate"], group["candidate_soft"])
        else:
            candidate_text = "有失败，不报告成功子集指标"
        lines.append(f"| {target} | {group['n_succeeded']}/{len(names)} | {display(group['baseline'], group['baseline_soft'])} | {candidate_text} |")
        groups[target] = group
    lines += ["", "此处全部为新制备、新 QM 尝试；失败配体保留，不替换成历史成功分数。",
        "独立检查受体全部固定；input_heavy_fixed 另核对每个配体重原子与输入相同。",
        "只有全部配体通过的集合报告指标；新实验失败不能解释成生产基线也失败。",
        "对照以冻结的结构假说定义，不按亲和力标签择优部署；生产默认未变。"]
    args.output.write_text("\n".join(lines)+"\n")
    write_json(args.output.with_suffix(".json"), {"groups": groups, "final_status": final,
        "protocol_sha256": digest(root / "frozen_protocol.json"), "control": control,
        "production_default_changed": False})
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
