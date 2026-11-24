"""Full-cohort 2x2 sidechain/H-restraint comparison, no best-arm selection."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from rdkit import Chem

from fmotpiefinal.core.bound_geometry import load_bound_ligand
from fmotpiefinal.core.ligand_hydrogen_geometry import audit_hch_angles
from fmotpiefinal.core.sidechain_policy import BACKBONE
from fmotpiefinal.core.structure import read_pdb
from fmotpiefinal.pipeline.target_stages import receptor_at_minimized_geometry
from fmotpiefinal.scoring import FACTOR, choose_reference
from report_local_sidechain_transfer import load_truth
from run_cdk8_solution_cycle import digest, metrics, write_json
from run_local_sidechain_transfer import directory_key


TARGETS = ("CHK1", "BACE", "Galectin3", "CDK8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", required=True, type=Path)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    source = args.source.resolve()
    root = args.run_root.resolve()
    plan = json.loads((source / "frozen_protocol.json").read_text())
    source_progress = json.loads((source / "progress.json").read_text())
    roots = {"fixed_Hfree": root / "hydrogens_unrestrained", "sidechain_Hfree": root / "sidechain_hydrogens_unrestrained"}
    progress = {}
    statuses = {}
    for arm, directory in roots.items():
        statuses[arm] = json.loads((directory / "final_status.json").read_text())
        if statuses[arm]["evaluation_sha256"] != digest(directory / "evaluation.json"):
            raise ValueError("factorial evaluation fingerprint changed")
        progress[arm] = json.loads((directory / "progress.json").read_text())
    groups = {}
    lines = ["# 侧链柔性 × 配体氢位置约束：完整 2×2 对照", "",
        "固定同一配体重原子 k10、受体侧链 k2.5、QM/PCM、分片与电荷；不把内部 MM 项叠加到 tPIE。",
        "H约束表示配体氢有 k10 位置约束；H自由表示仅移除该位置约束，不删氢、不改质子化。", "",
        "| 集合 | 固定受体＋H约束 | 侧链＋H约束 | 固定受体＋H自由 | 侧链＋H自由 |",
        "|---|---|---|---|---|"]
    for target in TARGETS:
        cohort = plan["cohorts"][target]
        names = cohort["names"]
        truth = load_truth(args.workspace.resolve(), target)
        if set(names) != set(truth):
            raise ValueError("full label inventory required")
        observed = [truth[name] for name in names]
        baseline = {name: {"status": "SUCCEEDED", "tpie": cohort["baseline"][name]["tpie"]} for name in names}
        all_rows = {"fixed_Hrestrained": baseline, "sidechain_Hrestrained": source_progress[target],
            **{arm: values[target] for arm, values in progress.items()}}
        molecules = {molecule.GetProp("_Name").strip(): molecule for molecule in Chem.SDMolSupplier(cohort["sdf"], removeHs=False) if molecule is not None}
        if set(molecules) != set(names):
            raise ValueError("full molecule inventory required")
        protein = read_pdb(cohort["protein"])
        original_receptor = np.asarray([atom.xyz for atom in protein.atoms])
        selected = set(map(tuple, cohort["local_residues"]))
        fixed_sidechain_indices = [index for index, atom in enumerate(protein.atoms)
            if (atom.chain, atom.residue_number, atom.residue_name) not in selected or atom.atom_name in BACKBONE]
        scores = [-baseline[name]["tpie"] for name in names]
        settings, reference = choose_reference(target, scores, None)
        width = FACTOR * settings["scaled_mad"]
        transform = lambda values: [settings["median"]+width*math.tanh((value-settings["median"])/width) for value in values]
        details = {}
        columns = []
        for arm, rows in all_rows.items():
            if set(rows) != set(names):
                raise ValueError("full factorial cell inventory required")
            failures = {name: row for name, row in rows.items() if row["status"] != "SUCCEEDED"}
            detail = {"n_requested": len(names), "n_succeeded": len(names)-len(failures), "failures": failures,
                "geometry": {}, "new_qm_in_this_report": 0}
            for name in names:
                if arm == "fixed_Hrestrained":
                    path = Path(cohort["complexes"][name])
                elif arm == "sidechain_Hrestrained":
                    path = source / target / directory_key(name) / "minimization/minimized_complex.pdb"
                else:
                    path = roots[arm] / target / directory_key(name) / "minimization/minimized_complex.pdb"
                if not path.exists():
                    detail["geometry"][name] = {"status": "MISSING_GEOMETRY"}
                    continue
                positions = np.asarray([atom.xyz for atom in load_bound_ligand(cohort["sdf"], name, path).atoms])
                receptor = receptor_at_minimized_geometry(protein, path)
                displacement = np.linalg.norm(np.asarray([atom.xyz for atom in receptor.atoms])-original_receptor, axis=1)
                indices = fixed_sidechain_indices if arm.startswith("sidechain_") else list(range(len(protein.atoms)))
                maximum = float(np.max(displacement[indices]))
                if maximum > 0.002:
                    raise ValueError(f"fixed receptor coordinate moved: {target}/{arm}/{name}")
                detail["geometry"][name] = {"status": "VERIFIED", "fixed_receptor_max_displacement_A": maximum,
                    "hydrogen_angle_audit": audit_hch_angles(molecules[name], positions), "complex_sha256": digest(path)}
            if failures:
                columns.append(f"{detail['n_succeeded']}/{len(names)}；不报子集指标")
            else:
                predicted = [-rows[name]["tpie"] for name in names]
                detail.update(metrics=metrics(observed, predicted), soft_metrics=metrics(observed, transform(predicted)),
                    ligand_minus_tpie=dict(zip(names, predicted, strict=True)))
                raw, soft = detail["metrics"], detail["soft_metrics"]
                columns.append(f"{soft['pearson']:.3f} / {raw['spearman']:.3f} / {raw['kendall']:.3f}")
            details[arm] = detail
        lines.append(f"| {target} | {' | '.join(columns)} |")
        differences = {}
        if all(not cell["failures"] for cell in details.values()):
            for metric in ("pearson", "spearman", "kendall"):
                get = lambda arm: details[arm]["metrics"][metric]
                differences[metric] = {"Hfree_effect_fixed": get("fixed_Hfree")-get("fixed_Hrestrained"),
                    "Hfree_effect_sidechain": get("sidechain_Hfree")-get("sidechain_Hrestrained"),
                    "factorial_difference_of_effects": get("sidechain_Hfree")-get("sidechain_Hrestrained")-get("fixed_Hfree")+get("fixed_Hrestrained")}
        groups[target] = {"cells": details, "descriptive_metric_differences": differences, "soft_reference": reference}
    lines += ["", "表格为软 Pearson / Spearman / Kendall；所有指标都来自完整集合。",
        "指标差的差仅为描述性统计，不是物理能量分解或因果显著性检验。",
        "基线及侧链＋H约束复用已验收计算；两个 H自由对照各尝试 85 次新 QM。",
        "本轮 H–C–H 审计是在发现 CAT-13i 异常后新增的诊断，不能冒充原先预注册选择规则。",
        "不得按 target 选择表中最高数值作为统一流程；不改变生产默认。"]
    args.output.write_text("\n".join(lines)+"\n")
    write_json(args.output.with_suffix(".json"), {"groups": groups, "statuses": statuses,
        "metric_order": ["soft_pearson", "spearman", "kendall"],
        "source_protocol_sha256": digest(source / "frozen_protocol.json"), "production_default_changed": False})
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
