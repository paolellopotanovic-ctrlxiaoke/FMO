"""Complete-cohort evaluation with recorded numerical recovery and pair audit."""
from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from rdkit import Chem

from run_cdk8_solution_cycle import ARMS, HARTREE_TO_KCAL, digest, metrics, write_json, group
from fmotpiefinal.core.solution_cycle import require_consistent_pcm_cycle, pcm_energy_inventory


def scores(tpie, gas, pcm, strain):
    return {"tpie_pcm": tpie, "three_state_gas_bound": gas,
        "three_state_pcm_bound": pcm, "three_state_pcm_free_min": pcm + strain,
        "tpie_plus_cycle_solvation": tpie + pcm - gas,
        "tpie_plus_ligand_strain": tpie + strain,
        "tpie_plus_solvation_and_strain": tpie + pcm - gas + strain}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--gas-run", type=Path)
    args = parser.parse_args()
    workspace, root = args.workspace.resolve(), args.run.resolve()
    frozen = json.loads((root / "frozen_protocol.json").read_text())
    baseline = workspace / "FMOtpiefinal_HIV1lab/runs/all_accelerated_robust_20261020/CDK8/independent_verification.json"
    tpie = {row["ligand"]: row["tpie"] for row in json.loads(baseline.read_text())["records"]}
    sdf = workspace / "vendor/public_binding_free_energy_benchmark/fep_benchmark_inputs/structure_inputs/misc/cdk8_koehler_ligands.sdf"
    molecules = {molecule.GetProp("_Name").strip(): molecule for molecule in Chem.SDMolSupplier(str(sdf), removeHs=False)}
    if set(molecules) != set(frozen["names"]) or set(tpie) != set(molecules):
        raise ValueError("complete exact cohort required")
    results, provenance, missing = {}, {}, []
    for key in frozen["deck_sha256"]:
        original = root / key / "result.json"
        selected = original
        for recovery_name in ("recovery", "recovery_v2", "recovery_v3", "recovery_v4", "recovery_v5", "recovery_v6"):
            candidate = root / recovery_name / key / "result.json"
            if candidate.exists():
                recovered = json.loads(candidate.read_text())
                if recovered["status"] == "SUCCEEDED":
                    selected = candidate
        result = json.loads(selected.read_text())
        if result["status"] != "SUCCEEDED":
            missing.append(key)
        results[key] = result
        provenance[key] = {"result": str(selected), "sha256": digest(selected),
                           "log_sha256": digest(selected.parent / "gamess.log"), "recovery": selected != original}
    if args.gas_run:
        gas_root = args.gas_run.resolve()
        for source in gas_root.glob("gas/*/*/result.json"):
            key = str(source.parent.relative_to(gas_root))
            results[key] = json.loads(source.read_text())
            if results[key]["status"] != "SUCCEEDED":
                missing.append(key)
            provenance[key] = {"result": str(source), "sha256": digest(source),
                "log_sha256": digest(source.parent / "gamess.log"), "recovery": False}
    rows = []
    pcm_audits = {}
    for name in frozen["names"]:
        state_keys = {"complex": f"pcm/{name}/complex", "receptor": "pcm/shared/receptor", "ligand": f"pcm/{name}/ligand"}
        pcm_audits[name] = require_consistent_pcm_cycle(
            {state: (root / key / "input.inp").read_text() for state, key in state_keys.items()},
            {state: (root / key / "gamess.log").read_text(errors="replace") for state, key in state_keys.items()})
        delta = {}
        for model in ("gas", "pcm"):
            keys = (f"{model}/{name}/complex", f"{model}/shared/receptor", f"{model}/{name}/ligand")
            if any(key not in results or results[key]["status"] != "SUCCEEDED" for key in keys):
                delta[model] = None
                continue
            complex_energy = results[f"{model}/{name}/complex"]["energy_hartree"]
            receptor_energy = results[f"{model}/shared/receptor"]["energy_hartree"]
            ligand_energy = results[f"{model}/{name}/ligand"]["energy_hartree"]
            delta[model] = (complex_energy - receptor_energy - ligand_energy) * HARTREE_TO_KCAL
        optimized = results[f"optimized/{name}/ligand"]
        bound = results[f"pcm/{name}/ligand"]["energy_hartree"]
        strain = (bound - optimized["energy_hartree"]) * HARTREE_TO_KCAL
        bound_deck = (root / f"pcm/{name}/ligand/input.inp").read_text()
        initial = np.array([list(map(float, line.split()[2:])) for line in group(bound_deck, "FMOXYZ")[1].strip().splitlines()])
        final = np.array([list(map(float, line.split()[2:])) for line in optimized["optimized_coordinates"]])
        if initial.shape != final.shape:
            raise ValueError("atom count changed during optimization")
        initial_center, final_center = initial.mean(axis=0), final.mean(axis=0)
        left, _, right = np.linalg.svd((initial - initial_center).T @ (final - final_center))
        correction = np.eye(3)
        correction[-1, -1] = np.linalg.det(left @ right)
        aligned = (initial - initial_center) @ left @ correction @ right + final_center
        available = {"tpie_pcm": tpie[name], "three_state_pcm_bound": delta["pcm"],
            "three_state_pcm_free_min": delta["pcm"] + strain,
            "tpie_plus_ligand_strain": tpie[name] + strain}
        if delta["gas"] is not None:
            available = scores(tpie[name], delta["gas"], delta["pcm"], strain)
        rows.append({"ligand": name, "experimental_dg": float(molecules[name].GetProp("r_exp_dg")),
            **available,
            "strain_kcal_mol": strain, "solvation_increment_kcal_mol": delta["pcm"] - delta["gas"] if delta["gas"] is not None else None,
            "optimized_aligned_allatom_rmsd_angstrom": float(np.sqrt(np.mean(np.sum((aligned - final)**2, axis=1))))})
    truth = [row["experimental_dg"] for row in rows]
    all_metrics = {arm: metrics(truth, [row[arm] for row in rows]) for arm in ARMS if all(arm in row for row in rows)}
    exploratory_names = ("three_state_electrostatic_only", "three_state_electrostatic_free",
        "tpie_plus_delta_electrostatic_solv", "tpie_plus_cavity", "tpie_plus_es_cavity")
    for row in rows:
        name = row["ligand"]
        inventories = pcm_audits[name]["states"]
        cavity_delta = inventories["complex"]["cavitation"] - inventories["receptor"]["cavitation"] - inventories["ligand"]["cavitation"]
        electrostatic_delta = inventories["complex"]["electrostatic"] - inventories["receptor"]["electrostatic"] - inventories["ligand"]["electrostatic"]
        optimized_inventory = pcm_energy_inventory((root / f"optimized/{name}/ligand/gamess.log").read_text())
        electronic_strain = row["strain_kcal_mol"] - inventories["ligand"]["cavitation"] + optimized_inventory["cavitation"]
        row.update(cavity_delta=cavity_delta, electrostatic_solvation_delta=electrostatic_delta,
            three_state_electrostatic_only=row["three_state_pcm_bound"] - cavity_delta,
            three_state_electrostatic_free=row["three_state_pcm_bound"] - cavity_delta + electronic_strain,
            tpie_plus_delta_electrostatic_solv=row["tpie_pcm"] + electrostatic_delta,
            tpie_plus_cavity=row["tpie_pcm"] + cavity_delta,
            tpie_plus_es_cavity=row["tpie_pcm"] + electrostatic_delta + cavity_delta)
    exploratory_metrics = {arm: metrics(truth, [row[arm] for row in rows]) for arm in exploratory_names}
    native_path = root / "native/execution.json"
    native = json.loads(native_path.read_text()) if native_path.exists() else {}
    if set(native) == set(frozen["names"]) and all(value["status"] == "SUCCEEDED" for value in native.values()):
        for row in rows:
            terms = native[row["ligand"]]["terms"]
            row["native_interaction"] = terms["interaction_kcal_mol"]
            row["native_interaction_plus_desolvation"] = terms["interaction_kcal_mol"] + terms["desolvation_kcal_mol"]
            row["native_binding"] = terms["binding_kcal_mol"]
            row["tpie_plus_native_desolvation"] = row["tpie_pcm"] + terms["desolvation_kcal_mol"]
            row["native_terms"] = terms
        for arm in ("native_interaction", "native_interaction_plus_desolvation", "native_binding", "tpie_plus_native_desolvation"):
            all_metrics[arm] = metrics(truth, [row[arm] for row in rows])
    pairs = []
    for first, second in combinations(rows, 2):
        experimental = second["experimental_dg"] - first["experimental_dg"]
        pairs.append({"first": first["ligand"], "second": second["ligand"], "experimental_delta": experimental,
            "arms": {arm: {"score_delta": second[arm] - first[arm],
                "correct_direction": (second[arm] - first[arm]) * experimental > 0} for arm in all_metrics}})
    report = {"schema": "cdk8-solution-cycle-complete-v1", "n": len(rows), "n_jobs": len(results),
        "n_recovered": sum(row["recovery"] for row in provenance.values()), "metrics": all_metrics,
        "failed_jobs": missing,
        "pcm_actual_energy_inventory_audit": pcm_audits,
        "posthoc_component_diagnostic_metrics": exploratory_metrics,
        "posthoc_warning": "Component subtraction from PCM-polarized states is diagnostic, not an independently optimized electrostatic-only model or free energy; no target-wise method selection",
        "native_complete_reconciled": len(native) == 10 and all(value["status"] == "SUCCEEDED" for value in native.values()),
        "primary_arm": frozen["primary_arm"], "primary_target_reached": all_metrics[frozen["primary_arm"]]["spearman"] >= 0.6,
        "target_reached": {arm: value["spearman"] >= 0.6 for arm, value in all_metrics.items()},
        "rows": rows, "pairs": pairs, "provenance": provenance,
        "source_sha256": {str(baseline): digest(baseline), str(sdf): digest(sdf)},
        "limitations": ["Ten-ligand development cohort; no held-out validation", "Single local free-ligand minimum, no entropy or state populations", "Fixed pocket and structural waters, not a full protein free-energy calculation", "Hybrid tPIE corrections are not a closed thermodynamic cycle", "No fitted weights or dropped cases"]}
    write_json(root / "final_evaluation.json", report)
    write_json(root / "final_status.json", {"status": "COMPLETE_10_LIGANDS_EVALUATED" if not missing else "PCM_COMPLETE_GAS_INCOMPLETE", "n": 10,
        "primary_target_reached": report["primary_target_reached"], "evaluation_sha256": digest(root / "final_evaluation.json")})
    lines = ["# CDK8：一致三态与溶液态参照", "", "全部10例；同一结构制备与分片；各项系数1，无拟合。", "",
        "|方法|Pearson|Spearman|Kendall|达到S≥0.6|", "|---|---:|---:|---:|---|"]
    lines.extend(f"|{arm}|{value['pearson']:.3f}|{value['spearman']:.3f}|{value['kendall']:.3f}|{report['target_reached'][arm]}|" for arm, value in all_metrics.items())
    lines.extend(["", "## 逐配体机制量", "", "|配体|实验DG|tPIE|气相三态|PCM三态|应变代价|游离优化三态|优化RMSD Å|", "|---|---:|---:|---:|---:|---:|---:|---:|"])
    lines.extend(f"|{row['ligand']}|{row['experimental_dg']:.2f}|{row['tpie_pcm']:.2f}|{row.get('three_state_gas_bound', float('nan')):.2f}|{row['three_state_pcm_bound']:.2f}|{row['strain_kcal_mol']:.2f}|{row['three_state_pcm_free_min']:.2f}|{row['optimized_aligned_allatom_rmsd_angstrom']:.3f}|" for row in rows)
    lines.extend(["", "## 限制", ""] + ["- " + value for value in report["limitations"]])
    lines.extend(["", "## 事后成分诊断（非预注册主要候选）", "", "|组合|Pearson|Spearman|Kendall|", "|---|---:|---:|---:|"])
    lines.extend(f"|{arm}|{value['pearson']:.3f}|{value['spearman']:.3f}|{value['kendall']:.3f}|" for arm, value in exploratory_metrics.items())
    lines.extend(["", "未通过作业：" + ", ".join(missing), "原生分解完整且与独立三态闭合：" + str(report["native_complete_reconciled"])])
    (root / "FINAL_REPORT.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({"metrics": all_metrics, "primary_target_reached": report["primary_target_reached"]}, indent=2))


if __name__ == "__main__":
    main()
