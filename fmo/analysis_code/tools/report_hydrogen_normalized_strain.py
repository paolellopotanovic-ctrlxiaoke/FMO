"""Independently verify constrained-QM logs before complete-cohort diagnostics."""
from __future__ import annotations

import argparse
import io
import json
from pathlib import Path
import sys

import numpy as np
from rdkit import Chem

from audit_cycle_ligand_chemistry import geometry_screen, stereochemistry
from fmotpiefinal.core.fmo_convergence import final_fmo_convergence
from fmotpiefinal.core.parser import parse_pie
from fmotpiefinal.core.runner import charge_aware_pair_outliers
from fmotpiefinal.core.solution_cycle import HARTREE_TO_KCAL_MOL, pcm_energy_inventory
from report_bace_extensions import paired_rank_summary
from report_local_sidechain_transfer import load_truth
from run_cdk8_solution_cycle import digest, write_json
from run_hydrogen_normalized_strain import evaluate, fixed_heavy_displacement, hydrogen_only_deck, transplant_ligand, xyz_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--prior", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.workspace.resolve()/"src"))
    from sophosqm_baseline.relative_total_energy import parse_converged_ligand_optimization
    frozen = json.loads((args.run/"frozen_protocol.json").read_text())
    evaluation = json.loads((args.run/"evaluation.json").read_text())
    prior = json.loads(args.prior.read_text())
    if digest(args.prior) != frozen["evaluation_sha256"] or evaluation["status"] != "COMPLETE_VERIFIED":
        raise ValueError("complete verified frozen calculation required")
    if digest(Path(__file__).with_name("run_hydrogen_normalized_strain.py")) != frozen["driver_sha256"]:
        raise ValueError("experiment driver changed")
    for name, fingerprint in frozen["dependencies"].items():
        if digest(Path(__file__).with_name(name)) != fingerprint:
            raise ValueError("experiment dependency changed")
    if "recovery" in frozen:
        if digest(Path(__file__).with_name("recover_hydrogen_normalized_strain.py")) != frozen["recovery"]["driver_sha256"]:
            raise ValueError("recovery driver changed")
    rows, geometry = {}, {}
    for name in frozen["names"]:
        task, saved = frozen["tasks"][name], evaluation["rows"][name]
        path = Path(task["input"])
        expected = hydrogen_only_deck(task["bound_ligand_deck"])
        if name in frozen.get("recovery", {}).get("names", []):
            from recover_hydrogen_normalized_strain import qa_recovery_deck
            expected = qa_recovery_deck(expected)
        if path.read_text() != expected:
            raise ValueError("frozen Cartesian constraints were not applied")
        for state, input_path in (("hydrogen", path), ("complex", path.parent.parent/"complex/input.inp")):
            proof = saved[state]
            if (digest(input_path) != proof["input_sha256"]
                    or digest(input_path.parent/"gamess.log") != proof["log_sha256"]
                    or digest(Path(task["runtime"]["gamess_dir"])/"gamess.00.x") != proof["executable_sha256"]):
                raise ValueError("QM input/log/executable fingerprint mismatch")
        optimized, energy = parse_converged_ligand_optimization((path.parent/"gamess.log").read_text(), ligand_deck=task["bound_ligand_deck"])
        if abs(energy-saved["hydrogen"]["energy_hartree"]) > 1.e-9 or optimized != saved["hydrogen"]["optimized_coordinates"]:
            raise ValueError("independent optimizer extraction mismatch")
        fixed_heavy_displacement(xyz_rows(task["bound_ligand_deck"]), optimized)
        molecule = Chem.MolFromMolBlock(task["molblock"], removeHs=False)
        initial = np.asarray([row.split()[2:] for row in xyz_rows(task["bound_ligand_deck"])], float)
        final = np.asarray([row.split()[2:] for row in optimized], float)
        chemistry = geometry_screen(molecule, final, stereochemistry(molecule, initial))
        if not chemistry["passed_screen"]:
            raise ValueError("independent chemical geometry screen failed")
        template = Path(task["complex_template"])
        if digest(template) != task["complex_template_sha256"]:
            raise ValueError("original complex changed")
        complex_path = path.parent.parent/"complex/input.inp"
        if complex_path.read_text() != transplant_ligand(template.read_text(), task["bound_ligand_deck"], optimized):
            raise ValueError("complex was modified beyond the frozen H-coordinate replacement")
        text = (complex_path.parent/"gamess.log").read_text()
        final_fmo_convergence(complex_path.read_text(), io.StringIO(text))
        fragment = task["runtime"]["ligand_fragment"]
        pairs = [pair for pair in parse_pie(complex_path.parent/"gamess.log") if fragment in (pair.fragment_i, pair.fragment_j)]
        outliers, _ = charge_aware_pair_outliers(pairs, complex_path.read_text(), 100.)
        if len(pairs) != fragment-1 or outliers or abs(sum(pair.value_kcal_mol for pair in pairs)-saved["complex"]["tpie"]) > 1.e-8:
            raise ValueError("independent ligand PIEDA quality/sum check failed")
        inventory = pcm_energy_inventory((path.parent/"gamess.log").read_text())
        source_inventory = prior["scores"][name]["pcm_audit"]["states"]["ligand"]
        internal_drop = (source_inventory["internal_hartree"]-inventory["internal_hartree"])*HARTREE_TO_KCAL_MOL
        solvent_drop = source_inventory["total"]-inventory["total"]
        source_energy = prior["rows"][frozen["source_case_tasks"][name]["ligand"]]["energy_hartree"]
        if abs(internal_drop+solvent_drop-(source_energy-energy)*HARTREE_TO_KCAL_MOL) > .01:
            raise ValueError("H-normalization PCM energy decomposition does not close")
        bond_changes = []
        for atom in molecule.GetAtoms():
            if atom.GetAtomicNum() == 1 and atom.GetDegree() == 1:
                index, parent = atom.GetIdx(), atom.GetNeighbors()[0].GetIdx()
                bond_changes.append({"hydrogen_index": index, "parent_element": molecule.GetAtomWithIdx(parent).GetSymbol(),
                    "initial_A": float(np.linalg.norm(initial[index]-initial[parent])),
                    "final_A": float(np.linalg.norm(final[index]-final[parent]))})
        geometry[name] = {"chemistry": chemistry, "hydrogen_bonds": bond_changes,
            "internal_energy_drop_kcal_mol": internal_drop, "pcm_response_drop_kcal_mol": solvent_drop}
        rows[name] = saved
    truth = load_truth(args.workspace.resolve(), frozen["target"])
    report = evaluate(rows, frozen, prior, truth)
    report.update(independent_geometry_evidence=geometry, driver_sha256=digest(__file__), frozen_protocol_sha256=digest(args.run/"frozen_protocol.json"))
    names = sorted(frozen["names"])
    observed = [truth[name] for name in names]
    comparisons = {"remove_H_component_from_original_correction": ("original_total_strain_score", "original_tpie_plus_residual_strain"),
        "recompute_H_geometry_and_residual_vs_original_correction": ("original_total_strain_score", "hydrogen_normalized_tpie_plus_residual_strain")}
    report["paired_comparisons"] = {label: paired_rank_summary(observed,
        [-report["scores"][name][before] for name in names], [-report["scores"][name][after] for name in names])
        for label, (before, after) in comparisons.items()}
    report["energy_distributions"] = {field: {"min": min(report["scores"][name][field] for name in names),
        "median": float(np.median([report["scores"][name][field] for name in names])),
        "max": max(report["scores"][name][field] for name in names)}
        for field in ("hydrogen_relaxation_energy", "residual_local_strain")}
    report["wall_seconds"] = evaluation["wall_seconds"]
    write_json(args.output.with_suffix(".json"), report)
    lines = [f"# {frozen['target']}：固定重原子、仅H的QM对照", "",
        f"完整{len(names)}例，{2*len(names)}项新QM独立验收。原受体、水、口袋、电荷、重原子坐标与数值阈值保持不变。", "",
        "| 预先冻结读数 | 原始Pearson | Spearman | Kendall |", "|---|---:|---:|---:|"]
    for arm, values in report["metrics"].items():
        lines.append(f"| {arm} | {values['pearson']:.4f} | {values['spearman']:.4f} | {values['kendall']:.4f} |")
    lines += ["", f"氢松弛与剩余局部应变分布（kcal/mol）：{report['energy_distributions']}。",
        "", "H_component + residual等于原总应变；不是新增两个独立的物理自由能项。",
        "孤立配体的H优化不代表结合环境H优化；tPIE重算用于检验移植这些H的接触代价，不直接修改默认制备。",
        "剩余项不等于纯重原子构象自由能，包含氢耦合、MM/QM模型差异和局部盆地依赖。负值也保留并解释，不截断或择优。",
        "这些都是开发集诊断，组合系数固定1，不按实验挑几何，不把新分数称为纯tPIE；完整配对不确定性见JSON。"]
    args.output.with_suffix(".md").write_text("\n".join(lines)+"\n")


if __name__ == "__main__":
    main()
