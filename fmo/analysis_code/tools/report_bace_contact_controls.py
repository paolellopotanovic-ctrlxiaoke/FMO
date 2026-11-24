"""Independently check both complete BACE controls before reporting metrics."""
from __future__ import annotations

import argparse
import io
import json
from pathlib import Path

import numpy as np

from fmotpiefinal.core.bound_geometry import load_bound_ligand
from fmotpiefinal.core.fmo_convergence import final_fmo_convergence
from fmotpiefinal.core.parser import parse_pie
from fmotpiefinal.core.sidechain_policy import BACKBONE
from fmotpiefinal.core.structure import read_pdb
from fmotpiefinal.pipeline.target_stages import receptor_at_minimized_geometry
from fmotpiefinal.scoring import FACTOR, choose_reference
from report_local_sidechain_transfer import load_truth
from run_cdk8_solution_cycle import digest, metrics, write_json
from run_cdk8_structure_controls import coordinate_array, topology_signature
from run_local_sidechain_transfer import directory_key, render_geometry


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.run.resolve()
    plan = json.loads((root / "frozen_protocol.json").read_text())
    progress = json.loads((root / "progress.json").read_text())
    truth = load_truth(args.workspace.resolve(), "BACE")
    groups = {}
    for arm, cohort in plan["cohorts"].items():
        names = cohort["names"]
        rows = progress[arm]["BACE"]
        if set(rows) != set(names) or set(truth) != set(names):
            raise ValueError("both complete 36-case controls required")
        failures = {name: row for name, row in rows.items() if row["status"] != "SUCCEEDED"}
        if failures:
            groups[arm] = {"status": "INCOMPLETE", "failures": failures}
            continue
        protein = read_pdb(cohort["protein"])
        initial = np.asarray([atom.xyz for atom in protein.atoms])
        selected = set(map(tuple, cohort["local_residues"]))
        fixed = [index for index, atom in enumerate(protein.atoms)
            if not ((arm == "contact_carboxylates" and
                (atom.chain, atom.residue_number, atom.residue_name) in selected and atom.atom_name not in BACKBONE)
                or (arm == "water_hydrogens_only" and atom.residue_name in {"HOH", "WAT"} and atom.element == "H"))]
        evidence = {}
        for name in names:
            folder = root / arm / "BACE" / directory_key(name)
            deck_path, log_path = folder / "input.inp", folder / "run/gamess.log"
            row = rows[name]
            if digest(deck_path) != row["input_sha256"] or digest(log_path) != row["log_sha256"]:
                raise ValueError("QM provenance changed")
            deck = deck_path.read_text()
            convergence = final_fmo_convergence(deck, io.StringIO(log_path.read_text()))
            receptor = receptor_at_minimized_geometry(protein, folder / "minimization/minimized_complex.pdb")
            ligand = load_bound_ligand(cohort["sdf"], name, folder / "minimization/minimized_complex.pdb")
            regenerated = render_geometry(receptor, {name: ligand}, name, set(cohort["pocket_keys"]))
            if topology_signature(deck) != topology_signature(Path(cohort["sources"][name]).read_text()):
                raise ValueError("fragment or charge topology changed")
            symbols, coordinates = coordinate_array(deck)
            expected_symbols, expected_coordinates = coordinate_array(regenerated)
            if symbols != expected_symbols or not np.array_equal(coordinates, expected_coordinates):
                raise ValueError("QM deck differs from actual optimized geometry")
            maximum = float(np.linalg.norm(np.asarray([atom.xyz for atom in receptor.atoms])-initial, axis=1)[fixed].max())
            if maximum > .002:
                raise ValueError("nonmobile receptor atom moved")
            ligand_shift = None
            if arm == "water_hydrogens_only":
                baseline = load_bound_ligand(cohort["sdf"], name, cohort["complexes"][name])
                ligand_shift = float(np.linalg.norm(np.asarray([atom.xyz for atom in ligand.atoms])-
                    np.asarray([atom.xyz for atom in baseline.atoms]), axis=1).max())
                if ligand_shift > .002:
                    raise ValueError("water-only control moved fixed ligand")
            count = row["pair_count"]+1
            pairs = [pair for pair in parse_pie(log_path) if count in (pair.fragment_i, pair.fragment_j)]
            if len(pairs) != count-1 or abs(sum(pair.value_kcal_mol for pair in pairs)-row["tpie"]) > 1e-9:
                raise ValueError("PIEDA inventory or sum mismatch")
            evidence[name] = {"max_fixed_receptor_displacement_A": maximum,
                "max_fixed_ligand_displacement_A": ligand_shift, "convergence": convergence,
                "input_sha256": digest(deck_path), "log_sha256": digest(log_path)}
        observed = [truth[name] for name in names]
        predicted = [-rows[name]["tpie"] for name in names]
        baseline = [-cohort["baseline"][name]["tpie"] for name in names]
        settings, _ = choose_reference("BACE", baseline, None)
        width = FACTOR*settings["scaled_mad"]
        soft = lambda values: settings["median"]+width*np.tanh((np.asarray(values)-settings["median"])/width)
        groups[arm] = {"status": "COMPLETE_INDEPENDENTLY_VERIFIED", "n": len(names),
            "metrics": metrics(observed, predicted), "soft_metrics": metrics(observed, soft(predicted)),
            "baseline_metrics": metrics(observed, baseline), "baseline_soft_metrics": metrics(observed, soft(baseline)),
            "geometry_and_qm": evidence}
    result = {"schema": "bace-contact-mechanism-verification-v1", "groups": groups,
        "frozen_protocol_sha256": digest(root / "frozen_protocol.json"),
        "production_rule_adopted": False, "no_partial_cohort_metrics": True}
    write_json(args.output.with_suffix(".json"), result)
    lines = ["# BACE 新增匹配对照：完整队列验收", "",
        "以下只比较完整tPIE。每组核对分片/电荷不变、真实优化坐标对应QM输入、固定原子位移、严格收敛与45对PIEDA之和。",
        "软Pearson使用同一冻结基线尺度，不拟合；Spearman/Kendall不受严格单调软压缩改变。", "",
        "| 方法 | 完整覆盖 | 原始Pearson | 软Pearson | Spearman | Kendall |", "|---|---:|---:|---:|---:|---:|"]
    first = next((value for value in groups.values() if "metrics" in value), None)
    if first:
        baseline = first["baseline_metrics"]
        lines.append(f"| 当前基线 | 36/36 | {baseline['pearson']:.3f} | {first['baseline_soft_metrics']['pearson']:.3f} | {baseline['spearman']:.3f} | {baseline['kendall']:.3f} |")
    for arm, value in groups.items():
        if "metrics" not in value:
            lines.append(f"| {arm} | 未完成，不报子集指标 | — | — | — | — |")
            continue
        values = value["metrics"]
        lines.append(f"| {arm} | 36/36 | {values['pearson']:.3f} | {value['soft_metrics']['pearson']:.3f} | {values['spearman']:.3f} | {values['kendall']:.3f} |")
    lines += ["", "接触羧酸侧链由全队列3.5Å接触选出ASP93/ASP289，不改变质子化；仅水氢组固定既有k10配体及水氧。",
        "本轮不是重新从原始PDB制备，也没有加入游离态/熵项。几何或数值检查通过不等于该方法能预测实验自由能。",
        "不会按这两个实验的较高相关系数自动换默认分支；若要推广，还需明确结构触发机制及跨集合回归。"]
    args.output.write_text("\n".join(lines)+"\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
