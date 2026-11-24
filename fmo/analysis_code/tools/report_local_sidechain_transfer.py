"""Audit completed local-sidechain transfer and report matched raw/soft metrics."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from rdkit import Chem

from build_model_comparison import COHORTS
from run_cdk8_solution_cycle import digest, metrics, write_json
from run_local_sidechain_transfer import directory_key


def load_truth(workspace, target):
    if target in {"CDK2", "JNK1", "P38", "TYK2"}:
        rows = csv.DictReader((workspace / "data/fmo87_cohort.csv").open())
        return {row["ligand_name"]: -float(row["experimental_dg_kcal_mol"])
                for row in rows if row["target"] == target}
    if target == "Galectin3":
        path = workspace / "FMOtpiefinal_HIV1lab/datasets/evaluation_only/Galectin3.csv"
        return {row["ligand"]: float(row["experimental_score_minus_dG"])
                for row in csv.DictReader(path.open())}
    sdf = workspace / "vendor/public_binding_free_energy_benchmark/fep_benchmark_inputs/structure_inputs" / (COHORTS[target][1] + "_ligands.sdf")
    return {molecule.GetProp("_Name").strip(): -float(molecule.GetProp("r_exp_dg"))
            for molecule in Chem.SDMolSupplier(str(sdf), removeHs=False) if molecule is not None}


def audit_geometry(cohort, run):
    from fmotpiefinal.core.structure import read_pdb
    from fmotpiefinal.pipeline.target_stages import receptor_at_minimized_geometry
    backbone = {"N", "CA", "C", "O", "OXT", "H", "H1", "H2", "H3", "HA", "HA2", "HA3"}
    protein = read_pdb(cohort["protein"])
    selected = {tuple(value) for value in cohort["local_residues"]}
    fixed = [index for index, atom in enumerate(protein.atoms)
             if (atom.chain, atom.residue_number, atom.residue_name) not in selected
             or atom.atom_name in backbone or atom.residue_name in {"HOH", "WAT"}]
    positions = np.asarray([atom.xyz for atom in protein.atoms])
    audit = {}
    for name in cohort["names"]:
        directory = run / cohort["target"] / directory_key(name)
        receptor = receptor_at_minimized_geometry(protein, directory / "minimization/minimized_complex.pdb")
        displacement = np.linalg.norm(np.asarray([atom.xyz for atom in receptor.atoms]) - positions, axis=1)
        maximum = float(max(displacement[fixed]))
        if maximum > 0.002:
            raise ValueError(f"fixed PDB atom changed: {cohort['target']}/{name}")
        provenance = json.loads((directory / "minimization/minimization.json").read_text())
        if provenance["max_fixed_receptor_displacement_angstrom"] > 0.0001:
            raise ValueError("full precision fixed atom changed")
        audit[name] = {"max_fixed_PDB_displacement_A": maximum,
            "max_receptor_PDB_displacement_A": float(max(displacement)),
            "max_fixed_full_precision_A": provenance["max_fixed_receptor_displacement_angstrom"],
            "n_mobile_receptor_atoms": provenance["n_mobile_receptor_atoms"]}
    return audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", required=True, type=Path)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    from fmotpiefinal.scoring import choose_reference, FACTOR
    workspace, run = args.workspace.resolve(), args.run.resolve()
    plan = json.loads((run / "frozen_protocol.json").read_text())
    final = json.loads((run / "final_status.json").read_text())
    if final["evaluation_sha256"] != digest(run / "evaluation.json"):
        raise ValueError("final evaluation fingerprint differs")
    progress = json.loads((run / "progress.json").read_text())
    report = {"final_status": final, "groups": {}, "production_default_changed": False,
        "parameter_selection_used_labels": False, "baseline_reused_not_recomputed": True}
    lines = ["# 局部侧链柔性：跨集合完整对照", "",
        "P 为原始 Pearson；软P沿用工具的冻结软压缩规则；S/K 为原始 Spearman/Kendall。", "",
        "| 集合 | 完成 | 基线 P / 软P / S / K | 侧链候选 P / 软P / S / K | ΔS |", "|---|---:|---|---|---:|"]
    for target, cohort in plan["cohorts"].items():
        names = cohort["names"]
        rows = progress.get(target, {})
        failures = {name: rows.get(name, {"status": "MISSING"}) for name in names
                    if rows.get(name, {}).get("status") != "SUCCEEDED"}
        group = {"n_expected": len(names), "n_succeeded": len(names)-len(failures), "failures": failures}
        if not failures:
            truth = load_truth(workspace, target)
            if set(truth) != set(names):
                raise ValueError("complete label identity mismatch")
            observed = [truth[name] for name in names]
            baseline = [-cohort["baseline"][name]["tpie"] for name in names]
            candidate = [-rows[name]["tpie"] for name in names]
            settings, provenance = choose_reference(target, baseline, None)
            if settings is None:
                raise ValueError("undefined softcap reference")
            center, width = settings["median"], FACTOR * settings["scaled_mad"]
            soft = lambda values: [center + width * math.tanh((value-center)/width) for value in values]
            group.update(baseline=metrics(observed, baseline), local_sidechains=metrics(observed, candidate),
                baseline_softcap=metrics(observed, soft(baseline)), local_softcap=metrics(observed, soft(candidate)),
                softcap_reference=provenance, geometry_audit=audit_geometry(cohort, run),
                rows=[{"ligand": name, "experimental_minus_dG": truth[name],
                       "baseline_minus_tpie": baseline[index], "local_minus_tpie": candidate[index]}
                      for index, name in enumerate(names)])
            display = lambda raw, transformed: f"{group[raw]['pearson']:.3f} / {group[transformed]['pearson']:.3f} / {group[raw]['spearman']:.3f} / {group[raw]['kendall']:.3f}"
            delta = group["local_sidechains"]["spearman"]-group["baseline"]["spearman"]
            lines.append(f"| {target} | {len(names)}/{len(names)} | {display('baseline', 'baseline_softcap')} | {display('local_sidechains', 'local_softcap')} | {delta:+.3f} |")
        else:
            lines.append(f"| {target} | {group['n_succeeded']}/{len(names)} | — | 有失败，不报告子集指标 | — |")
        report["groups"][target] = group
    original = ("CDK2", "JNK1", "P38", "TYK2")
    if all(target in report["groups"] and not report["groups"][target]["failures"] for target in original):
        report["original87_weighted"] = {method: {field: sum(report["groups"][target]["n_expected"] *
            report["groups"][target][method][field] for target in original) / 87
            for field in ("pearson", "spearman", "kendall")}
            for method in ("baseline", "local_sidechains", "baseline_softcap", "local_softcap")}
    lines += ["", "同一局部规则原样迁移；不按集合标签在两组中择优。",
        "基线为已独立验收的冻结结果；候选全部为新制备、新QM。",
        "软压缩使用同一基线参考尺度，不重新寻找有利参数。",
        "此表不包含 OpenBind；其多晶体聚合需要单独接续。"]
    write_json(args.output.with_suffix(".json"), report)
    args.output.write_text("\n".join(lines) + "\n")
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
