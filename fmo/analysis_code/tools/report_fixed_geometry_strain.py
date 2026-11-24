"""Independent complete-cohort verification for radial and graph NH2 controls."""
from __future__ import annotations

import argparse
import io
import json
from pathlib import Path
import sys

import numpy as np
from rdkit import Chem

from audit_cycle_ligand_chemistry import geometry_screen, stereochemistry
from fmotpiefinal.core.conjugated_amine_geometry import induced_conjugated_nh2, reconstruct_flagged_nh2
from fmotpiefinal.core.fmo_convergence import final_fmo_convergence
from fmotpiefinal.core.parser import parse_pie
from fmotpiefinal.core.runner import charge_aware_pair_outliers
from fmotpiefinal.core.solution_cycle import HARTREE_TO_KCAL_MOL, require_consistent_pcm_cycle
from report_bace_extensions import paired_rank_summary
from report_local_sidechain_transfer import load_truth
from run_cdk8_bound_qm import replace_coordinates
from run_cdk8_solution_cycle import digest, metrics, write_json
from run_direction_preserving_hydrogen_control import radial_coordinates
from run_hydrogen_normalized_strain import fixed_heavy_displacement, transplant_ligand, xyz_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--prior", type=Path, required=True)
    parser.add_argument("--cycle-protocol", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--hydrogen-run", type=Path)
    parser.add_argument("--registry", type=Path)
    args = parser.parse_args()
    sys.path.insert(0, str(args.workspace.resolve()/"src"))
    from sophosqm_baseline.relative_total_energy import parse_full_energy
    frozen = json.loads((args.run/"frozen_protocol.json").read_text())
    evaluation = json.loads((args.run/"evaluation.json").read_text())
    prior = json.loads(args.prior.read_text())
    cycle = json.loads(args.cycle_protocol.read_text())
    graph = frozen["schema"] == "graph-induced-conjugated-nh2-repair-control-v1"
    names = sorted(frozen["names"])
    if (evaluation["status"] != "COMPLETE_VERIFIED" or set(evaluation["rows"]) != set(names)
            or set(cycle["names"]) != set(names) or digest(args.prior) != frozen["prior_sha256"]):
        raise ValueError("complete matched frozen references required")
    driver = "run_conjugated_amine_control.py" if graph else "run_direction_preserving_hydrogen_control.py"
    if digest(Path(__file__).with_name(driver)) != frozen["driver_sha256"]:
        raise ValueError("control driver changed")
    if graph:
        if args.registry is None or digest(args.registry) != frozen["registry_sha256"] or digest(args.cycle_protocol) != frozen["source_protocol_sha256"]:
            raise ValueError("graph source fingerprints mismatch")
        core = Path(__file__).parents[1]/"src/fmotpiefinal/core/conjugated_amine_geometry.py"
        if digest(core) != frozen["graph_core_sha256"] or digest(Path(__file__).with_name("run_direction_preserving_hydrogen_control.py")) != frozen["executor_sha256"]:
            raise ValueError("graph core/executor changed")
        registry = json.loads(args.registry.read_text())["cohorts"][frozen["target"]]
        molecules = {molecule.GetProp("_Name").strip(): molecule for molecule in Chem.SDMolSupplier(registry["sdf"], removeHs=False) if molecule is not None}
        hydrogen = None
    else:
        if args.hydrogen_run is None or digest(args.hydrogen_run/"evaluation.json") != frozen["hydrogen_evaluation_sha256"] or digest(args.hydrogen_run/"frozen_protocol.json") != frozen["hydrogen_protocol_sha256"]:
            raise ValueError("constrained H source fingerprints mismatch")
        for name, fingerprint in frozen["dependencies"].items():
            if digest(Path(__file__).with_name(name)) != fingerprint:
                raise ValueError("control dependency changed")
        hydrogen = json.loads((args.hydrogen_run/"evaluation.json").read_text())
        molecules = {name: Chem.MolFromMolBlock(frozen["tasks"][name]["molblock"], removeHs=False) for name in names}
    scores, geometry = {}, {}
    for name in names:
        task, saved = frozen["tasks"][name], evaluation["rows"][name]
        source_keys = cycle["case_tasks"][name]
        ligand_reference = Path(cycle["tasks"][source_keys["ligand"]]["input"]).read_text()
        template = Path(cycle["tasks"][source_keys["complex"]]["input"]).read_text()
        initial_rows = xyz_rows(ligand_reference)
        before = np.asarray([row.split()[2:] for row in initial_rows], float)
        molecule = molecules[name]
        if graph:
            diagnostics = induced_conjugated_nh2(molecule, molecule.GetConformer().GetPositions(), before)
            points = reconstruct_flagged_nh2(molecule, before, diagnostics)
            expected_rows = [f"{row.split()[0]} {row.split()[1]} {point[0]:.10f} {point[1]:.10f} {point[2]:.10f}" for row, point in zip(initial_rows, points, strict=True)]
        else:
            diagnostics = []
            expected_rows = radial_coordinates(molecule, initial_rows, hydrogen["rows"][name]["hydrogen"]["optimized_coordinates"])
        expected = {"ligand": replace_coordinates(ligand_reference, expected_rows),
                    "complex": transplant_ligand(template, ligand_reference, expected_rows)}
        energies = {}
        decks, logs = {}, {}
        for state in ("ligand", "complex"):
            path = Path(task[state+"_input"])
            log = path.parent/"gamess.log"
            if (path.read_text() != expected[state] or digest(path) != saved[state]["input_sha256"]
                    or digest(log) != saved[state]["log_sha256"]
                    or digest(Path(task["runtime"]["gamess_dir"])/"gamess.00.x") != saved[state]["executable_sha256"]):
                raise ValueError("input/log/executable or reconstructed geometry mismatch")
            decks[state], logs[state] = path.read_text(), log.read_text()
            energies[state] = parse_full_energy(logs[state], input_text=decks[state])
            if abs(energies[state]-saved[state]["energy_hartree"]) > 1.e-9:
                raise ValueError("independent full energy extraction mismatch")
        receptor = Path(cycle["tasks"][source_keys["receptor"]]["input"])
        receptor_row = prior["rows"][source_keys["receptor"]]
        if digest(receptor) != receptor_row["input_sha256"] or digest(receptor.parent/"gamess.log") != receptor_row["log_sha256"]:
            raise ValueError("fixed receptor evidence changed")
        decks["receptor"], logs["receptor"] = receptor.read_text(), (receptor.parent/"gamess.log").read_text()
        audit = require_consistent_pcm_cycle(decks, logs)
        final_fmo_convergence(decks["complex"], io.StringIO(logs["complex"]))
        points = np.asarray([row.split()[2:] for row in expected_rows], float)
        shift = fixed_heavy_displacement(initial_rows, expected_rows)
        screen = geometry_screen(molecule, points, stereochemistry(molecule, before))
        if not screen["passed_screen"]:
            raise ValueError("independent chemical screen failed")
        fragment = task["runtime"]["ligand_fragment"]
        pairs = [pair for pair in parse_pie(Path(task["complex_input"]).parent/"gamess.log") if fragment in (pair.fragment_i, pair.fragment_j)]
        outliers, _ = charge_aware_pair_outliers(pairs, decks["complex"], 100.)
        tpie = sum(pair.value_kcal_mol for pair in pairs)
        if len(pairs) != fragment-1 or outliers or abs(tpie-saved["complex"]["tpie"]) > 1.e-8:
            raise ValueError("independent PIEDA inventory/gate/sum check failed")
        old = prior["scores"][name]
        bound = prior["rows"][source_keys["ligand"]]["energy_hartree"]
        change = (bound-energies["ligand"])*HARTREE_TO_KCAL_MOL
        residual = old["ligand_strain"]-change
        if graph:
            scores[name] = {"matched_tpie": old["matched_tpie"], "original_total_strain_score": old["tpie_plus_ligand_strain"],
                           "graph_repaired_tpie": tpie, "graph_repaired_local_strain": residual,
                           "graph_repaired_tpie_plus_local_strain": tpie+residual}
        else:
            scores[name] = {"matched_tpie": old["matched_tpie"], "original_total_strain_score": old["tpie_plus_ligand_strain"],
                "radial_component": change, "radial_residual": residual, "radial_tpie": tpie,
                "radial_tpie_plus_residual_strain": tpie+residual,
                "original_tpie_plus_radial_component": old["matched_tpie"]+change,
                "original_tpie_plus_radial_residual": old["matched_tpie"]+residual}
        if any(abs(value-evaluation["scores"][name][field]) > 1.e-8 for field, value in scores[name].items()):
            raise ValueError("independent score recomputation mismatch")
        geometry[name] = {"chemical_screen": screen, "max_heavy_displacement_A": shift,
                          "graph_diagnostics": diagnostics, "pcm_consistent": audit["consistent"]}
    truth = load_truth(args.workspace.resolve(), frozen["target"])
    if set(truth) != set(names):
        raise ValueError("complete truth inventory required")
    observed = [truth[name] for name in names]
    result = {"status": "COMPLETE_INDEPENDENTLY_VERIFIED", "target": frozen["target"], "n": len(names),
        "metrics": {arm: metrics(observed, [-scores[name][arm] for name in names]) for arm in frozen["arms"]},
        "scores": scores, "geometry": geometry, "production_default_changed": False,
        "driver_sha256": digest(__file__), "frozen_protocol_sha256": digest(args.run/"frozen_protocol.json")}
    pure = "graph_repaired_tpie" if graph else "radial_tpie"
    combined = "graph_repaired_tpie_plus_local_strain" if graph else "radial_tpie_plus_residual_strain"
    result["paired_comparisons"] = {label: paired_rank_summary(observed,
        [-scores[name][before] for name in names], [-scores[name][after] for name in names])
        for label, before, after in (("pure_vs_matched_tpie", "matched_tpie", pure), ("combined_vs_original_strain", "original_total_strain_score", combined))}
    write_json(args.output.with_suffix(".json"), result)
    lines = [f"# {frozen['target']}：{'图触发共轭NH2修复' if graph else '保留X-H方向的键长对照'}", "",
        f"完整{len(names)}例独立核验；几何、输入/日志/二进制哈希、严格收敛、完整PIEDA和一致PCM三态均复核。", "",
        "| 固定方法 | 原始Pearson | Spearman | Kendall |", "|---|---:|---:|---:|"]
    for arm, values in result["metrics"].items():
        lines.append(f"| {arm} | {values['pearson']:.4f} | {values['spearman']:.4f} | {values['kendall']:.4f} |")
    lines += ["", "仅几何候选，不是QM优化驻点或完整自由能；不得按target选较好方法，不自动替换生产默认。",
              "开发集配对bootstrap不能消除重复试验选择偏差，完整逐例证据与不确定性见JSON。"]
    args.output.with_suffix(".md").write_text("\n".join(lines)+"\n")


if __name__ == "__main__":
    main()
