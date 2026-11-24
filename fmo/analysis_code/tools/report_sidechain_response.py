"""Evaluate frozen structural repair evidence before reading experimental labels."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from rdkit import Chem

from fmotpiefinal.core.bound_geometry import load_bound_ligand
from fmotpiefinal.core.sidechain_response import decide_target, diagnose_case
from fmotpiefinal.core.structure import read_pdb
from fmotpiefinal.pipeline.target_stages import receptor_at_minimized_geometry
from run_cdk8_solution_cycle import digest, write_json
from run_local_sidechain_transfer import directory_key


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--verification", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    source = args.source.resolve()
    plan = json.loads((source / "frozen_protocol.json").read_text())
    verification = json.loads(args.verification.read_text())
    if verification["frozen_protocol_sha256"] != digest(source / "frozen_protocol.json"):
        raise ValueError("candidate verification protocol mismatch")
    groups = {}
    source_hashes = {str(source / "frozen_protocol.json"): digest(source / "frozen_protocol.json"),
        str(args.verification.resolve()): digest(args.verification)}
    lines = ["# 无标签的侧链修复响应验收", "",
        "规则先冻结，结构决定不读取实验亲和力、tPIE 或相关系数。此处是有限机制研究对照，不改变默认流程。", "",
        "| 集合 | 基线新增酰胺扭曲数 | 修复通过数 | 新增扭曲数 | 碰撞计数增加例数 | 决定 |",
        "|---|---:|---:|---:|---:|---|"]
    for target, cohort in plan["cohorts"].items():
        protein = read_pdb(cohort["protein"])
        molecules = {molecule.GetProp("_Name").strip(): molecule for molecule in
            Chem.SDMolSupplier(cohort["sdf"], removeHs=False) if molecule is not None}
        if set(molecules) != set(cohort["names"]):
            raise ValueError("molecule inventory mismatch")
        verified = verification["groups"][target]
        geometry_verified = not cohort["local_residues"] or (
            len(verified.get("independent_geometry_audit", {})) == len(cohort["names"]))
        cases = {}
        for name in cohort["names"]:
            baseline_path = Path(cohort["complexes"][name])
            candidate_path = source / target / directory_key(name) / "minimization/minimized_complex.pdb"
            if not cohort["local_residues"]:
                candidate_path = baseline_path
            if not candidate_path.exists():
                raise ValueError(f"complete candidate geometry required: {target}/{name}")
            for path in (baseline_path, candidate_path, Path(cohort["protein"]), Path(cohort["sdf"])):
                source_hashes[str(path)] = digest(path)
            if digest(baseline_path) != cohort["input_sha256"][str(baseline_path)]:
                raise ValueError("baseline complex fingerprint changed")
            baseline_ligand = np.asarray([atom.xyz for atom in load_bound_ligand(cohort["sdf"], name, baseline_path).atoms])
            candidate_ligand = np.asarray([atom.xyz for atom in load_bound_ligand(cohort["sdf"], name, candidate_path).atoms])
            candidate_receptor = receptor_at_minimized_geometry(protein, candidate_path)
            cases[name] = diagnose_case(molecules[name], molecules[name].GetConformer().GetPositions(),
                baseline_ligand, candidate_ligand, protein, candidate_receptor)
        decision = decide_target(cases, cohort["names"], geometry_verified=geometry_verified,
            qm_verified=not verified["failures"] and verified["n_succeeded"] == len(cohort["names"]))
        count = sum(len(case["baseline_new_amide_flags"]) for case in cases.values())
        repaired = sum(flag["repaired"] for case in cases.values() for flag in case["baseline_new_amide_flags"])
        new = sum(len(case["new_candidate_amide_flags"]) for case in cases.values())
        increased = sum(not case["severe_overlap_count_not_increased"] for case in cases.values())
        lines.append(f"| {target} | {count} | {repaired} | {new} | {increased} | {decision['action']} |")
        groups[target] = {"decision": decision, "cases": cases}
    result = {"schema": "structure-response-acceptance-v1", "groups": groups,
        "source_sha256": source_hashes, "code_sha256": {
            str(Path(__file__).resolve()): digest(__file__),
            str(Path(__file__).parents[1] / "src/fmotpiefinal/core/sidechain_response.py"):
                digest(Path(__file__).parents[1] / "src/fmotpiefinal/core/sidechain_response.py")},
        "labels_read_by_decision": False, "scores_read_by_decision": False,
        "production_default_changed": False}
    lines += ["", "接纳仅说明通过预先固定的制备变形修复机制及安全检查，不保证亲和力排序改善。",
        "无触发不说明没有问题；KEEP_BASELINE 不是对基线物理正确性的认证。",
        "碰撞检查覆盖完整受体，包括骨架与保留水；不把所有较短极性接触当成碰撞。",
        "已知开发集合上的规则不能被称为新的外部盲测。"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines)+"\n")
    write_json(args.output.with_suffix(".json"), result)
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
