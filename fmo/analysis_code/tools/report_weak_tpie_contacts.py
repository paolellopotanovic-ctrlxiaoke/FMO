"""Audit existing tPIE contributions and chemical inventory without changing scores."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import re

import numpy as np
from rdkit import Chem

from fmotpiefinal.core.parser import parse_pie
from fmotpiefinal.core.structure import read_pdb
from run_cdk8_solution_cycle import digest, write_json


def fragment_charges(deck_path):
    deck = deck_path.read_text()
    names_match = re.search(r"FRGNAM\(1\)\s*=(.*?)\bINDAT\(1\)", deck, re.S)
    charges_match = re.search(r"ICHARG\(1\)\s*=(.*?)\bFRGNAM\(1\)", deck, re.S)
    if names_match is None or charges_match is None:
        raise ValueError("fragment charge/name inventory missing")
    names = [value.strip() for value in names_match[1].strip().strip(",").split(",")]
    charges = [int(value) for value in re.split(r"[\s,]+", charges_match[1].strip().strip(","))]
    if len(names) != len(charges) or len(names) != len(set(names)):
        raise ValueError("fragment charge/name counts mismatch")
    return dict(zip(names, charges, strict=True))


def named_pairs(deck_path, expected):
    deck = deck_path.read_text()
    count = int(re.search(r"\bNFRAG\s*=\s*(\d+)", deck)[1])
    match = re.search(r"FRGNAM\(1\)\s*=(.*?)\bINDAT\(1\)", deck, re.S)
    if match is None:
        raise ValueError("fragment names missing")
    names = [value.strip() for value in match[1].strip().strip(",").split(",")]
    if len(names) != count or len(set(names)) != count:
        raise ValueError("fragment name inventory mismatch")
    log = deck_path.parent / "run/gamess.log"
    pairs = [pair for pair in parse_pie(log) if count in (pair.fragment_i, pair.fragment_j)]
    values = {names[(pair.fragment_j if pair.fragment_i == count else pair.fragment_i)-1]:
        pair.value_kcal_mol for pair in pairs}
    if len(values) != count-1 or abs(sum(values.values())-expected) > 1e-6:
        raise ValueError("baseline tPIE did not reproduce from all named pairs")
    return values, digest(log)


def chemical_inventory(protein_path, molecules):
    protein = read_pdb(protein_path)
    grouped = defaultdict(list)
    for atom in protein.atoms:
        grouped[(atom.chain, atom.residue_number, atom.residue_name)].append(atom)
    ligand_positions = np.concatenate([molecule.GetConformer().GetPositions()[
        [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1]]
        for molecule in molecules.values()])
    sites = []
    for residue, atoms in grouped.items():
        if residue[2] not in {"ASP", "ASH", "GLU", "GLH", "HIS", "HID", "HIE", "HIP"}:
            continue
        heavy = np.asarray([atom.xyz for atom in atoms if atom.element.upper() != "H"])
        distance = float(np.min(np.linalg.norm(heavy[:, None, :]-ligand_positions[None, :, :], axis=2)))
        if distance > 4.0:
            continue
        oxygens = [atom for atom in atoms if atom.atom_name in {"OD1", "OD2", "OE1", "OE2"}]
        hydrogens = [atom for atom in atoms if atom.element.upper() == "H"]
        oh = [{"oxygen": oxygen.atom_name, "hydrogen": hydrogen.atom_name,
            "distance_A": float(np.linalg.norm(np.asarray(oxygen.xyz)-hydrogen.xyz))}
            for oxygen in oxygens for hydrogen in hydrogens
            if np.linalg.norm(np.asarray(oxygen.xyz)-hydrogen.xyz) < 1.3]
        sites.append({"residue": list(residue), "minimum_ligand_heavy_distance_A": distance,
            "atom_names": [atom.atom_name for atom in atoms], "short_carboxyl_oxygen_H_contacts": oh,
            "interpretation": "Atom inventory and distance evidence only; no pKa or experimental protonation assignment"})
    return {"ligand_formal_charges": {name: Chem.GetFormalCharge(molecule) for name, molecule in molecules.items()},
        "nearby_acidic_histidine_sites": sites}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    plan = json.loads((args.run / "frozen_protocol.json").read_text())
    groups = {}
    lines = ["# 弱集合 tPIE 接触与化学状态清单", "",
        "只分析已验收基线；不改分数、不挑残基作为新预测项。范围预设为 CHK1、BACE、CDK8、HIV1PROT 及 Galectin3 回归。", "",
        "以下为逐片段总 PIE，不是 ES 或 SOLV 拟合。每例加总必须复现原 tPIE。", ""]
    cohorts = plan["cohorts"][plan["arms"][0]]
    for target in ("CHK1", "BACE", "CDK8", "HIV1PROT", "Galectin3"):
        cohort = cohorts[target]
        molecules = {molecule.GetProp("_Name").strip(): molecule for molecule in
            Chem.SDMolSupplier(cohort["sdf"], removeHs=False) if molecule is not None}
        if set(molecules) != set(cohort["names"]):
            raise ValueError("molecular identity mismatch")
        rows = {}
        charge_inventory = {}
        for name in cohort["names"]:
            source = Path(cohort["sources"][name])
            if digest(source) != cohort["input_sha256"][str(source)]:
                raise ValueError("baseline deck fingerprint mismatch")
            pairs, log_hash = named_pairs(source, cohort["baseline"][name]["tpie"])
            if log_hash != cohort["baseline"][name]["log_sha256"]:
                raise ValueError("baseline log fingerprint mismatch")
            rows[name] = pairs
            charge_inventory[name] = fragment_charges(source)
        fragments = sorted(next(iter(rows.values())))
        if any(set(row) != set(fragments) for row in rows.values()):
            raise ValueError("whole cohort fragment inventory mismatch")
        summary = [{"fragment": fragment, "mean_pie": float(np.mean([row[fragment] for row in rows.values()])),
            "std_pie": float(np.std([row[fragment] for row in rows.values()])),
            "range_pie": float(np.ptp([row[fragment] for row in rows.values()]))} for fragment in fragments]
        inventory = chemical_inventory(cohort["protein"], molecules)
        receptor_charges = {fragment: charge_inventory[cohort["names"][0]][fragment] for fragment in fragments}
        for name in cohort["names"]:
            if any(charge_inventory[name][fragment] != charge for fragment, charge in receptor_charges.items()):
                raise ValueError("receptor fragment charge changed within baseline cohort")
        for site in inventory["nearby_acidic_histidine_sites"]:
            chain, number, residue = site["residue"]
            site["qm_fragment_charge"] = receptor_charges.get(chain+residue+number)
        totals = [sum(row.values()) for row in rows.values()]
        total_statistics = {"mean": float(np.mean(totals)), "std": float(np.std(totals)),
            "range": float(np.ptp(totals))}
        groups[target] = {"n": len(rows), "chemical_inventory": inventory,
            "fragment_summary": summary, "ligand_fragment_pie": rows,
            "total_tpie_statistics": total_statistics, "qm_fragment_charges": charge_inventory}
        lines += [f"## {target}", "", f"完整复现 {len(rows)}/{len(cohort['names'])} 例 tPIE。",
            "配体形式电荷分布：" + str(dict(sorted(Counter(inventory['ligand_formal_charges'].values()).items()))) + "。",
            f"总 tPIE：均值 {total_statistics['mean']:.3f}，标准差 {total_statistics['std']:.3f}，极差 {total_statistics['range']:.3f} kcal/mol。", "",
            "| 变化最大的片段 | 平均 PIE (kcal/mol) | 队列标准差 | 极差 |", "|---|---:|---:|---:|"]
        for row in sorted(summary, key=lambda row: row["std_pie"], reverse=True)[:6]:
            lines.append(f"| {row['fragment']} | {row['mean_pie']:.3f} | {row['std_pie']:.3f} | {row['range_pie']:.3f} |")
        lines += ["", "4 Å 范围酸性/组氨酸位点（编号为当前受体实际编号）："]
        for site in inventory["nearby_acidic_histidine_sites"]:
            residue = ":".join(site["residue"])
            hydrogens = site["short_carboxyl_oxygen_H_contacts"]
            lines.append(f"- {residue}；最短重原子距离 {site['minimum_ligand_heavy_distance_A']:.3f} Å；QM 片段电荷 {site['qm_fragment_charge']}；羧基短 O–H 接触：{hydrogens}。")
        lines += [""]
    lines += ["## 结论边界", "",
        "贡献变化大只说明该接触主导分值差异，不证明其错误，也不证明修改它能提分。",
        "平均吸引能很大但队列变化很小的接触主要提供共同偏移；共同常数不会改变排序。片段方差不能直接相加当总方差，因为片段间可能协变。",
        "羧基没有短 O–H 是当前输入的证据，不能独自证明真实实验态去质子化；须核对制备日志、QM 电荷与微观状态模型。",
        "同样，发现一个带氢羧基不能证明整套质子化状态正确。此诊断不会按标签翻转质子化或删除位点。",
        "需要结构/化学状态对照后再检验完整队列及回归集合，不将诊断假说当作已确认根因。"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines)+"\n")
    write_json(args.output.with_suffix(".json"), {"groups": groups, "labels_used": False,
        "scores_changed": False, "frozen_protocol_sha256": digest(args.run / "frozen_protocol.json")})
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
