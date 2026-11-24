"""Trace BACE misranking to chemistry/contact evidence, without modifying scores."""
from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import rankdata

from fmotpiefinal.core.bound_geometry import load_bound_ligand
from fmotpiefinal.core.structure import read_pdb
from report_local_sidechain_transfer import load_truth
from run_cdk8_solution_cycle import array, digest, group, metrics, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = args.source.resolve()
    cohort = json.loads((source / "frozen_protocol.json").read_text())["cohorts"]["BACE"]
    for path, expected in cohort["input_sha256"].items():
        if digest(path) != expected:
            raise ValueError("source evidence changed")
    molecules = {molecule.GetProp("_Name").strip(): molecule for molecule in
                 Chem.SDMolSupplier(cohort["sdf"], removeHs=False) if molecule is not None}
    truth = load_truth(args.workspace.resolve(), "BACE")
    names = cohort["names"]
    if set(names) != set(molecules) or set(names) != set(truth):
        raise ValueError("full BACE chemical and experimental inventory required")
    protein = read_pdb(cohort["protein"])
    predicted = {name: -cohort["baseline"][name]["tpie"] for name in names}
    ranks = {name: float(value) for name, value in zip(names,
             rankdata([predicted[name] for name in names])-rankdata([truth[name] for name in names]), strict=True)}
    fragment_values = {}
    contacts = []
    fingerprint = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fingerprints = {name: fingerprint.GetFingerprint(Chem.RemoveHs(molecule)) for name, molecule in molecules.items()}
    for name in names:
        deck = Path(cohort["sources"][name]).read_text()
        fragment_names = array(group(deck, "FMO")[1], "FRGNAM", "INDAT")
        for pair in cohort["runtime"][name]["pairs"]:
            receptor_index = pair["fragment_j"] if pair["fragment_i"] == len(fragment_names) else pair["fragment_i"]
            fragment_values.setdefault(fragment_names[receptor_index-1], {})[name] = pair["value_kcal_mol"]
        ligand = load_bound_ligand(cohort["sdf"], name, cohort["complexes"][name])
        heavy = np.asarray([atom.xyz for atom in ligand.atoms if atom.element != "H"])
        case = {"ligand": name, "experimental_minus_dG": truth[name], "predicted_minus_tpie": predicted[name],
                "predicted_minus_experimental_rank": ranks[name], "formal_charge": Chem.GetFormalCharge(molecules[name])}
        for number in ("93", "289", "502"):
            atoms = [atom for atom in protein.atoms if atom.chain == "A" and atom.residue_number == number]
            if not atoms:
                raise ValueError("expected source residue missing")
            distance = float(np.linalg.norm(np.asarray([atom.xyz for atom in atoms if atom.element != "H"])[:, None, :]
                                           - heavy[None, :, :], axis=-1).min())
            case[f"residue_{number}_nearest_heavy_A"] = distance
        contacts.append(case)
    full_fragments = {key: values for key, values in fragment_values.items() if set(values) == set(names)}
    if len(full_fragments) != len(fragment_values):
        raise ValueError("fragment identities differ across baseline cohort")
    fragments = [{"fragment": key, "mean_kcal_mol": float(np.mean(list(values.values()))),
                  "std_kcal_mol": float(np.std(list(values.values()))),
                  "range_kcal_mol": float(np.ptp(list(values.values())))} for key, values in full_fragments.items()]
    fragments.sort(key=lambda row: row["std_kcal_mol"], reverse=True)
    reversed_pairs = []
    for first, second in combinations(names, 2):
        experimental_difference = truth[first]-truth[second]
        predicted_difference = predicted[first]-predicted[second]
        similarity = DataStructs.TanimotoSimilarity(fingerprints[first], fingerprints[second])
        if similarity >= .7 and abs(experimental_difference) >= .5 and experimental_difference*predicted_difference < 0:
            contributions = [{"fragment": key, "delta_tpie_kcal_mol": values[first]-values[second]}
                             for key, values in full_fragments.items()]
            contributions.sort(key=lambda row: abs(row["delta_tpie_kcal_mol"]), reverse=True)
            reversed_pairs.append({"first": first, "second": second, "morgan_tanimoto": similarity,
                "delta_experimental_minus_dG": experimental_difference,
                "delta_predicted_minus_tpie": predicted_difference,
                "delta_tpie_components": contributions})
    reversed_pairs.sort(key=lambda row: abs(row["delta_experimental_minus_dG"]), reverse=True)
    protonation = []
    for number in ("93", "289"):
        atoms = [atom for atom in protein.atoms if atom.chain == "A" and atom.residue_number == number]
        bonds = [{"hydrogen": hydrogen.atom_name, "oxygen": oxygen.atom_name,
                  "distance_A": float(np.linalg.norm(np.asarray(hydrogen.xyz)-np.asarray(oxygen.xyz)))}
                 for oxygen in atoms if oxygen.atom_name in {"OD1", "OD2"}
                 for hydrogen in atoms if hydrogen.element == "H"
                 and np.linalg.norm(np.asarray(hydrogen.xyz)-np.asarray(oxygen.xyz)) < 1.25]
        protonation.append({"residue": f"A:ASP:{number}", "sidechain_OH_bonds": bonds})
    result = {"schema": "bace-structural-misranking-evidence-v1", "n": len(names),
        "source_protocol_sha256": digest(source / "frozen_protocol.json"),
        "metrics": metrics([truth[name] for name in names], [predicted[name] for name in names]),
        "cases": contacts, "dyad_protonation_geometry": protonation,
        "fragment_variation": fragments, "similar_reversed_pairs": reversed_pairs,
        "pair_thresholds": {"morgan_radius": 2, "tanimoto_min": .7, "abs_experimental_gap_kcal_min": .5},
        "interpretation": "Exploratory attribution, not causal validation or a new ranking score. No fragment deleted, no MW correction, no label-directed preparation."}
    write_json(args.output.with_suffix(".json"), result)
    lines = ["# BACE 结构、化学状态与错排归因", "",
        "完整36例；以下是已有分数的事后诊断，不据此给每个配体挑计算路径。", "",
        "## 已核实", "",
        "- 配体均为+1；同一受体ASP93不带羧基氢，ASP289有羧基氢。不是简单的两个Asp都去质子化。",
        "- 保留的TLA属于远离配体的已注册MM组分；下面列出距离，不把它默认删掉或强塞入QM。",
        "- 强平均相互作用不等于能解释配体间亲和力变化。下表列出片段间相互作用在36例中的变化。", "",
        "| 片段 | 平均PIE | 标准差 | 全距（kcal/mol） |", "|---|---:|---:|---:|"]
    for row in fragments[:12]:
        lines.append(f"| {row['fragment']} | {row['mean_kcal_mol']:.3f} | {row['std_kcal_mol']:.3f} | {row['range_kcal_mol']:.3f} |")
    lines += ["", "## 接触距离", "", "| 受体残基 | 36例最短接触范围/A |", "|---|---:|"]
    for number in ("93", "289", "502"):
        values = [row[f"residue_{number}_nearest_heavy_A"] for row in contacts]
        lines.append(f"| {number} | {min(values):.3f}–{max(values):.3f} |")
    lines += ["", "## 高相似配体中的错排", "",
        "Morgan半径2、Tanimoto≥0.7、实验间隔≥0.5kcal/mol，仅作诊断；不把配体名称前缀当作相同骨架。",
        "实验和预测均使用越大越强方向；末列为PIE(前者)−PIE(后者)，与−tPIE差反号。", "",
        "| 配体对 | 相似度 | 实验差 | −tPIE差 | 变化最大的三片段 |", "|---|---:|---:|---:|---|"]
    for row in reversed_pairs[:12]:
        fragments_text = "; ".join(f"{item['fragment']}:{item['delta_tpie_kcal_mol']:+.2f}" for item in row['delta_tpie_components'][:3])
        lines.append(f"| {row['first']} / {row['second']} | {row['morgan_tanimoto']:.3f} | {row['delta_experimental_minus_dG']:+.3f} | {row['delta_predicted_minus_tpie']:+.3f} | {fragments_text} |")
    lines += ["", "## 下一步的可证伪对照", "",
        "1. 只放开全队列接触羧酸侧链（k2.5），其余受体与水固定；配体仍为k10。检验催化接触的局部适应。",
        "2. 固定既有k10配体、蛋白和水氧，仅优化现有水氢。检验水氢朝向，不混入配体构象变化。",
        "两组均重做同一45个受体片段的QM、保留36例、报告全部失败；不自动采用排序较高的一组。",
        "若两组仍差，不能据此宣称质子化正确或自由能贡献已抵消；后续应在明确实验pH、微观状态与权重后，检验配对质子化状态和结合/游离态采样。",
        "当前证据仅排除了部分简单解释，尚未建立像CDK8酰胺制备变形那样被修复实验证实的BACE根因。"]
    args.output.write_text("\n".join(lines)+"\n")
    print(json.dumps({"n": len(names), "similar_reversed_pairs": len(reversed_pairs), "metrics": result['metrics']}))


if __name__ == "__main__":
    main()
