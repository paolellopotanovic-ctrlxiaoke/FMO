"""Transfer the frozen CDK8 local-sidechain control without target-wise tuning."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import csv
import io
import json
from pathlib import Path
import re
import time

import numpy as np
from scipy.spatial import cKDTree

from build_model_comparison import COHORTS
from run_cdk8_bound_qm import replace_coordinates
from run_cdk8_solution_cycle import digest, group, metrics, write_json
from run_cdk8_structure_controls import coordinate_array, topology_signature


def directory_key(name):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", name)


def render_geometry(protein, ligands, name, keys):
    from fmotpiefinal.core.structure import Protein, residue_key
    from fmotpiefinal.core.fragmentation import make_facio_like_fragments
    from fmotpiefinal.core.fmo_input import render_fmo2_pieda
    pocket = Protein(tuple(atom for atom in protein.atoms if residue_key(atom) in keys))
    capped, fragments, metadata = make_facio_like_fragments(pocket, ligands[name],
        reference_protein=protein, merge_disulfides=True, repair_peptide_numbering=True)
    return render_fmo2_pieda(capped, ligands[name], fragments, dftb_parameter_dir="../params",
        fragment_bonds=metadata["hop_pairs_gamess_indices"], full_energy=True)


def register_cohort(baseline, target):
    from fmotpiefinal.core.structure import read_pdb, residue_key
    from fmotpiefinal.core.bound_geometry import load_bound_ligand, shared_pocket
    root = baseline / target
    verification_path = root / "independent_verification.json"
    verification = json.loads(verification_path.read_text())
    if verification["status"] != "COMPLETE_INDEPENDENTLY_VERIFIED":
        raise ValueError("baseline not independently verified")
    records = {row["ligand"]: row for row in verification["records"]}
    if len(records) != COHORTS[target][0]:
        raise ValueError("incomplete reference cohort")
    candidates = {}
    for path in (root / "run").rglob("input.inp"):
        if path.parent.parent.name == "qm":
            candidates.setdefault(digest(path), []).append(path)
    sources = {}
    for name, record in records.items():
        matches = candidates.get(record["input_sha256"], [])
        matches = [path for path in matches if path.parent.name == directory_key(name)]
        if not matches:
            raise ValueError(f"verified source deck missing: {name}")
        sources[name] = min(matches, key=lambda path: len(str(path)))
        if digest(sources[name].parent / "run/gamess.log") != record["log_sha256"]:
            raise ValueError(f"verified source log changed: {name}")
    source_roots = {path.parents[2] for path in sources.values()}
    if len(source_roots) != 1:
        raise ValueError("baseline does not use one registered target protocol")
    source_root = source_roots.pop()
    pdbs = list(source_root.glob("*_prepared.pdb"))
    if len(pdbs) != 1:
        raise ValueError("expected one registered prepared receptor")
    protein_path = pdbs[0]
    protein = read_pdb(protein_path)
    sdf = source_root / "ligands_no_properties.sdf"
    from rdkit import Chem
    molecules = list(Chem.SDMolSupplier(str(sdf), removeHs=False))
    if any(molecule is None or len(molecule.GetPropNames()) for molecule in molecules):
        raise ValueError("property-free valid source molecules required")
    minimized_root = root / "run/qm/full_target_fmo/minimized"
    if (source_root / "minimized").exists():
        minimized_root = source_root / "minimized"
    names = list(records)
    complexes = {name: minimized_root / name / "minimized_complex.pdb" for name in names}
    ligands = {name: load_bound_ligand(sdf, name, complexes[name]) for name in names}
    pocket = shared_pocket(protein, list(ligands.values()), 5.0)
    keys = {residue_key(atom) for atom in pocket.atoms}
    heavy = np.asarray([atom.xyz for ligand in ligands.values() for atom in ligand.atoms
                        if atom.element.upper() != "H"])
    tree = cKDTree(heavy)
    selected = sorted({(atom.chain, atom.residue_number, atom.residue_name)
        for atom in protein.atoms if atom.record == "ATOM" and atom.element.upper() != "H"
        and float(tree.query(atom.xyz)[0]) <= 4.0})
    if not selected:
        raise ValueError("empty local sidechain region")
    runtime = {}
    for name in names:
        source = sources[name].read_text()
        regenerated = render_geometry(protein, ligands, name, keys)
        if topology_signature(source) != topology_signature(regenerated):
            raise ValueError(f"baseline fragment/cap topology mismatch: {name}")
        symbols, positions = coordinate_array(source)
        other_symbols, other_positions = coordinate_array(regenerated)
        if symbols != other_symbols or not np.array_equal(positions, other_positions):
            raise ValueError(f"baseline cap/atom coordinates mismatch: {name}")
        result_path = sources[name].parent / "run/result.json"
        runtime[name] = json.loads(result_path.read_text())
        if runtime[name]["status"] != "SUCCEEDED":
            raise ValueError("baseline runtime not successful")
    parameters = (next(iter(sources.values())).parent / "params").resolve(strict=True)
    hashes = {str(path): digest(path) for path in
              [verification_path, protein_path, sdf, *sources.values(), *complexes.values()]}
    hashes.update({str(path): digest(path) for path in protein_path.parent.glob(protein_path.stem + ".*.json")})
    return {"target": target, "names": names, "protein": str(protein_path), "sdf": str(sdf),
        "sources": {name: str(path) for name, path in sources.items()},
        "complexes": {name: str(path) for name, path in complexes.items()},
        "local_residues": selected, "pocket_keys": sorted(keys), "parameters": str(parameters),
        "runtime": runtime, "input_sha256": hashes, "baseline": records,
        "baseline_regeneration": "exact atom coordinates, caps, topology, charge and HOP match"}


def prepare_job(payload):
    from fmotpiefinal.core.minimization import minimize_complex
    cohort, name, output = payload
    directory = Path(output) / cohort["target"] / directory_key(name)
    started = time.monotonic()
    try:
        result = minimize_complex(cohort["protein"], cohort["sdf"], name, directory / "minimization",
            restraint_k_kj_mol_a2=10.0, move_protein=False,
            local_flexible_residues=tuple(map(tuple, cohort["local_residues"])),
            receptor_restraint_k_kj_mol_a2=2.5)
        result = {"status": "SUCCEEDED", "complex": str(result.complex_pdb)}
    except Exception as error:
        result = {"status": "FAILED", "error": f"{type(error).__name__}: {error}"}
    result["wall_seconds"] = time.monotonic() - started
    write_json(directory / "preparation_status.json", result)
    return cohort["target"], name, result


def score_job(cohort, name, output, ncores):
    from fmotpiefinal.core.structure import read_pdb
    from fmotpiefinal.core.bound_geometry import load_bound_ligand
    from fmotpiefinal.pipeline.target_stages import receptor_at_minimized_geometry, verify_fmo_deck_heavy_geometry
    from fmotpiefinal.core.runner import run_gamess
    from fmotpiefinal.core.fmo_convergence import final_fmo_convergence
    from fmotpiefinal.core.parser import parse_pie
    parent = output / cohort["target"] / directory_key(name)
    try:
        prepared = json.loads((parent / "preparation_status.json").read_text())
        if prepared["status"] != "SUCCEEDED":
            raise ValueError("preparation failed; not silently skipped")
        protein = read_pdb(cohort["protein"])
        receptor = receptor_at_minimized_geometry(protein, Path(prepared["complex"]))
        ligand = load_bound_ligand(cohort["sdf"], name, Path(prepared["complex"]))
        regenerated = render_geometry(receptor, {name: ligand}, name, set(cohort["pocket_keys"]))
        source = Path(cohort["sources"][name]).read_text()
        if topology_signature(source) != topology_signature(regenerated):
            raise ValueError("local optimization changed fragment topology")
        deck = replace_coordinates(source, group(regenerated, "FMOXYZ")[1].strip().splitlines())
        verify_fmo_deck_heavy_geometry(deck)
        path = parent / "input.inp"
        if path.exists() and path.read_text() != deck:
            raise ValueError("refusing changed existing input")
        path.write_text(deck)
        count = int(re.search(r"NFRAG=(\d+)", deck)[1])
        runtime = cohort["runtime"][name]
        result_path = parent / "run/result.json"
        if result_path.exists():
            record = json.loads(result_path.read_text())
            if record["input_sha256"] != digest(path):
                raise ValueError("cached input mismatch")
        else:
            record = run_gamess(path, parent / "run", gamess_dir=runtime["gamess_dir"],
                library_dir=runtime["library_dir"], ligand_fragment=count,
                ncores=ncores, timeout_seconds=1800, parameter_dir=cohort["parameters"])
        if record["status"] != "SUCCEEDED":
            raise ValueError(f"QM status: {record['status']}")
        log = parent / "run/gamess.log"
        convergence = final_fmo_convergence(deck, io.StringIO(log.read_text()))
        pairs = [pair for pair in parse_pie(log) if count in (pair.fragment_i, pair.fragment_j)]
        if len(pairs) != count - 1:
            raise ValueError("incomplete ligand PIE inventory")
        score = sum(pair.value_kcal_mol for pair in pairs)
        result = {"status": "SUCCEEDED", "tpie": score, "baseline_tpie": cohort["baseline"][name]["tpie"],
            "pair_count": len(pairs), "convergence": convergence, "input_sha256": digest(path),
            "log_sha256": digest(log), "qm_wall_seconds": record.get("wall_seconds"),
            "pairs": [{"i": pair.fragment_i, "j": pair.fragment_j, "pie": pair.value_kcal_mol} for pair in pairs]}
    except Exception as error:
        result = {"status": "FAILED", "error": f"{type(error).__name__}: {error}"}
    write_json(parent / "summary.json", result)
    return cohort["target"], name, result


def evaluate(workspace, baseline, results, output):
    from rdkit import Chem
    originals = list(csv.DictReader((workspace / "data/fmo87_cohort.csv").open()))
    groups = {}
    for target, rows in results.items():
        expected = COHORTS[target][0]
        group_result = {"n_expected": expected, "n_succeeded": sum(row["status"] == "SUCCEEDED" for row in rows.values()),
                        "failures": {name: row for name, row in rows.items() if row["status"] != "SUCCEEDED"}}
        if len(rows) == expected and group_result["n_succeeded"] == expected:
            if target in {"CDK2", "JNK1", "P38", "TYK2"}:
                truth = {row["ligand_name"]: float(row["experimental_dg_kcal_mol"])
                         for row in originals if row["target"] == target}
            else:
                sdf = workspace / "vendor/public_binding_free_energy_benchmark/fep_benchmark_inputs/structure_inputs" / (COHORTS[target][1] + "_ligands.sdf")
                truth = {molecule.GetProp("_Name").strip(): float(molecule.GetProp("r_exp_dg"))
                         for molecule in Chem.SDMolSupplier(str(sdf), removeHs=False) if molecule is not None}
                if target == "Galectin3":
                    labels = workspace / "FMOtpiefinal_HIV1lab/datasets/evaluation_only/Galectin3.csv"
                    truth = {row["ligand"]: -float(row["experimental_score_minus_dG"])
                             for row in csv.DictReader(labels.open())}
            if set(truth) != set(rows):
                raise ValueError(f"complete identity mismatch: {target}")
            names = sorted(rows)
            observed = [truth[name] for name in names]
            group_result["local_sidechains"] = metrics(observed, [rows[name]["tpie"] for name in names])
            group_result["baseline"] = metrics(observed, [rows[name]["baseline_tpie"] for name in names])
        groups[target] = group_result
    result = {"groups": groups, "production_default_changed": False, "partial_cohort_metrics_forbidden": True}
    write_json(output / "evaluation.json", result)
    lines = ["# 匹配结构对照：跨集合完整验证", "", "具体处理、实际执行路径与新 QM/复用数量见本目录 frozen_protocol.json 和 progress.json；不按测试集调参。", "",
             "| 集合 | 完成 | 基线 P / S / K | 当前对照 P / S / K |", "|---|---:|---|---|"]
    for target, values in groups.items():
        display = lambda key: " / ".join(f"{values[key][field]:.3f}" for field in ("pearson", "spearman", "kendall")) if key in values else "未完成，不报告子集指标"
        lines.append(f"| {target} | {values['n_succeeded']}/{values['n_expected']} | {display('baseline')} | {display('local_sidechains')} |")
    lines += ["", "基线为独立验收的冻结结果；无触发路径的复用必须在实际执行记录中明确标记。",
              "保留每个集合实际已验收的 QM 数值设置；不是将 CDK8 的收敛器强加到其它集合。",
              "这是一项候选流程回归；不会按实验标签自动挑选处理分支。"]
    (output / "RESULTS.md").write_text("\n".join(lines) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", required=True, type=Path)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--targets", nargs="+", default=[target for target in COHORTS if target != "CDK8"])
    parser.add_argument("--prepare-workers", type=int, default=6)
    parser.add_argument("--qm-workers", type=int, default=12)
    parser.add_argument("--ncores", type=int, default=2)
    parser.add_argument("--register-only", action="store_true")
    args = parser.parse_args()
    if min(args.prepare_workers, args.qm_workers, args.ncores) < 1:
        parser.error("positive workers/cores required")
    if len(set(args.targets)) != len(args.targets) or set(args.targets) - set(COHORTS):
        parser.error("unique registered targets required")
    workspace, baseline, output = args.workspace.resolve(), args.baseline.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    plan_path = output / "frozen_protocol.json"
    if plan_path.exists():
        plan = json.loads(plan_path.read_text())
        if plan["targets"] != args.targets:
            raise ValueError("frozen target list changed")
    else:
        cohorts, errors = {}, {}
        for target in args.targets:
            try:
                cohorts[target] = register_cohort(baseline, target)
                print("REGISTER", target, len(cohorts[target]["names"]), flush=True)
            except Exception as error:
                errors[target] = f"{type(error).__name__}: {error}"
                print("REGISTER_FAILED", target, errors[target], flush=True)
        write_json(output / "registration.json", {"cohorts": cohorts, "errors": errors})
        if errors:
            raise ValueError(f"registration failed; no requested target silently excluded: {errors}")
        plan = {"targets": args.targets, "cohorts": cohorts, "labels_used_for_generation": False,
            "driver_sha256": digest(__file__), "baseline": str(baseline),
            "local_radius_angstrom": 4.0, "local_receptor_k": 2.5, "ligand_k": 10.0,
            "production_default_changed": False, "baseline_policy": "reuse independently verified complete cohort",
            "code_sha256": {str(path): digest(path) for path in [Path(__file__),
                Path(__file__).parents[1] / "src/fmotpiefinal/core/minimization.py"]}}
        write_json(plan_path, plan)
    for path, expected in plan["code_sha256"].items():
        if digest(path) != expected:
            raise ValueError("frozen implementation changed")
    for cohort in plan["cohorts"].values():
        for path, expected in cohort["input_sha256"].items():
            if digest(path) != expected:
                raise ValueError("frozen source changed")
    if args.register_only:
        return
    started = time.monotonic()
    results = {target: {} for target in args.targets}
    with ProcessPoolExecutor(max_workers=args.prepare_workers) as preparation, ThreadPoolExecutor(max_workers=args.qm_workers) as quantum:
        futures = {}
        for cohort in plan["cohorts"].values():
            for name in cohort["names"]:
                path = output / cohort["target"] / directory_key(name) / "preparation_status.json"
                if path.exists() and json.loads(path.read_text())["status"] == "SUCCEEDED":
                    futures[quantum.submit(score_job, cohort, name, output, args.ncores)] = "qm"
                else:
                    futures[preparation.submit(prepare_job, (cohort, name, str(output)))] = "prepare"
        while futures:
            future = next(as_completed(futures))
            stage = futures.pop(future)
            target, name, record = future.result()
            print(stage.upper(), target, name, record["status"], record.get("tpie", record.get("error", "")), flush=True)
            if stage == "prepare" and record["status"] == "SUCCEEDED":
                futures[quantum.submit(score_job, plan["cohorts"][target], name, output, args.ncores)] = "qm"
            else:
                results[target][name] = record
                write_json(output / "progress.json", results)
                if len(results[target]) == len(plan["cohorts"][target]["names"]):
                    evaluate(workspace, baseline, results, output)
    evaluation = evaluate(workspace, baseline, results, output)
    complete = all(group["n_succeeded"] == group["n_expected"] for group in evaluation["groups"].values())
    write_json(output / "final_status.json", {"status": "COMPLETE" if complete else "INCOMPLETE_WITH_FAILURES",
        "wall_seconds": time.monotonic() - started, "n_requested": sum(len(cohort["names"]) for cohort in plan["cohorts"].values()),
        "evaluation_sha256": digest(output / "evaluation.json")})


if __name__ == "__main__":
    main()
