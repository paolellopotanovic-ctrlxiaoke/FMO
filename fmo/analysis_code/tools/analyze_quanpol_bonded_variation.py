"""Attribute MM bonded fluctuations on all frozen frames without fitting affinities."""
from __future__ import annotations

import argparse
import bz2
import json
import math
from pathlib import Path
import statistics

import numpy as np

from analyze_endpoint_internal_energy import evaluate, parameterized_context, physical_atom_maps, read_topology
from report_quanpol_multiframe import describe
from run_matched_openfe import digest, write_json


def physical_constraint_pairs(system, saved_particle_indices, ligand_saved_indices):
    saved = np.asarray(saved_particle_indices, dtype=int)
    ligand = np.asarray(ligand_saved_indices, dtype=int)
    if saved.ndim != 1 or len(set(saved)) != len(saved) or np.any(saved < 0) or np.any(saved >= system.getNumParticles()):
        raise ValueError("invalid saved-to-system particle mapping")
    if np.any(ligand < 0) or np.any(ligand >= len(saved)):
        raise ValueError("invalid ligand saved-particle indices")
    reverse = {int(saved[original]): index for index, original in enumerate(ligand)}
    constrained = set()
    for offset in range(system.getNumConstraints()):
        first, second, _ = system.getConstraintParameters(offset)
        if first in reverse and second in reverse:
            constrained.add(tuple(sorted((reverse[first], reverse[second]))))
    return constrained


def bonded_observables(system, molecule, coordinates, constrained_pairs):
    import openmm
    from openmm import unit

    terms = []
    for force in system.getForces():
        if isinstance(force, openmm.HarmonicBondForce):
            for offset in range(force.getNumBonds()):
                first, second, equilibrium, constant = force.getBondParameters(offset)
                distance = float(np.linalg.norm(coordinates[first] - coordinates[second]))
                equilibrium = equilibrium.value_in_unit(unit.angstrom)
                coefficient = .5 * constant.value_in_unit(unit.kilocalorie_per_mole / unit.angstrom**2)
                terms.append({"id": f"bond_{offset}", "kind": "bond", "atoms": [first, second],
                              "elements": [molecule.GetAtomWithIdx(index).GetSymbol() for index in (first, second)],
                              "has_hydrogen": any(molecule.GetAtomWithIdx(index).GetAtomicNum() == 1 for index in (first, second)),
                              "constrained_in_source": tuple(sorted((first, second))) in constrained_pairs,
                              "observed": distance, "equilibrium": equilibrium, "unit": "angstrom",
                              "energy_kcal_mol": coefficient * (distance - equilibrium)**2})
        elif isinstance(force, openmm.HarmonicAngleForce):
            for offset in range(force.getNumAngles()):
                first, middle, last, equilibrium, constant = force.getAngleParameters(offset)
                left, right = coordinates[first] - coordinates[middle], coordinates[last] - coordinates[middle]
                cosine = np.dot(left, right) / (np.linalg.norm(left) * np.linalg.norm(right))
                angle = math.acos(float(np.clip(cosine, -1, 1)))
                equilibrium = equilibrium.value_in_unit(unit.radian)
                coefficient = .5 * constant.value_in_unit(unit.kilocalorie_per_mole / unit.radian**2)
                terms.append({"id": f"angle_{offset}", "kind": "angle", "atoms": [first, middle, last],
                              "elements": [molecule.GetAtomWithIdx(index).GetSymbol() for index in (first, middle, last)],
                              "has_hydrogen": any(molecule.GetAtomWithIdx(index).GetAtomicNum() == 1 for index in (first, middle, last)),
                              "observed": math.degrees(angle), "equilibrium": math.degrees(equilibrium), "unit": "degree",
                              "energy_kcal_mol": coefficient * (angle - equilibrium)**2})
    if not all(math.isfinite(term["energy_kcal_mol"]) for term in terms):
        raise ValueError("nonfinite bonded energy")
    return terms


def describe_correlation(values, corrections):
    summary = describe(values)
    deviation = statistics.stdev(values) if len(values) > 1 else 0
    correction_sd = statistics.stdev(corrections) if len(corrections) > 1 else 0
    summary["pearson_with_internal_delta_descriptive"] = (
        float(np.corrcoef(values, corrections)[0, 1]) if deviation > 1e-10 and correction_sd > 1e-10 else None)
    return summary


def summarize(rows):
    states = {}
    for leg in ("bound", "solution"):
        for endpoint in (0, 10):
            selected = [row for row in rows if row["leg"] == leg and row["endpoint"] == endpoint]
            if not selected:
                raise ValueError("missing state")
            corrections = [row["internal_correction_kcal_mol"] for row in selected]
            components = {}
            for kind in ("bond", "angle"):
                for has_hydrogen in (False, True):
                    name = f"{kind}_{'H' if has_hydrogen else 'heavy'}"
                    values = [sum(term["energy_kcal_mol"] for term in row["terms"]
                                  if term["kind"] == kind and term["has_hydrogen"] == has_hydrogen) for row in selected]
                    components[name] = describe_correlation(values, corrections)
            for name in selected[0]["mm_components"]:
                components[name] = describe_correlation([row["mm_components"][name] for row in selected], corrections)
            by_term = []
            for reference in selected[0]["terms"]:
                entries = [next(term for term in row["terms"] if term["id"] == reference["id"]) for row in selected]
                metadata = {key: reference[key] for key in ("id", "kind", "atoms", "elements", "has_hydrogen", "unit")}
                if reference["kind"] == "bond":
                    metadata["constrained_in_source"] = reference["constrained_in_source"]
                by_term.append({**metadata, "energy": describe_correlation([term["energy_kcal_mol"] for term in entries], corrections),
                                "geometry": describe(term["observed"] for term in entries),
                                "equilibrium": reference["equilibrium"]})
            states[f"{leg}_{endpoint}"] = {"n_frames": len(selected), "components": components, "all_bonded_terms": by_term}
    return states


def render(report):
    lines = ["# CDK8 固定帧配体 MM 键合项审计", "",
             "全部 60 帧；参数来自同一冻结缓存。H 表示含氢的键/角，heavy 表示全部为重原子的键/角。",
             "下表为 MM 项本身的波动，不是 QM 修正的物理贡献比例，也不是新亲和力分数。", "",
             "| 状态 | 键 heavy SD | 键 H SD | 角 heavy SD | 角 H SD | 扭转 SD | 内部非键 SD |",
             "|---|---:|---:|---:|---:|---:|---:|"]
    for state, details in report["states"].items():
        values = details["components"]
        names = ("bond_heavy", "bond_H", "angle_heavy", "angle_H", "PeriodicTorsionForce", "NonbondedForce")
        cells = [f"{values[name]['sample_sd']:.4f}" for name in names]
        lines.append(f"| {state} | " + " | ".join(cells) + " |")
    lines.extend(["", "单位：kcal/mol。所有逐键/逐角的几何分布、平衡值、能量和描述性相关完整保留于 JSON。",
                  "相关仅在各单态内部计算，不跨不同分子的绝对 QM 能量混算。MM 项在 QM−MM 差值里带负号，机械相关不能当作因果归因。", "",
                  f"- 最大 MM bond/angle 分项闭合误差：{report['maximum_bonded_closure_kcal_mol']:.3g} kcal/mol。",
                  "- 本报告不删除约束键、不修改力场，也不根据相关性拟合亲和力。", ""])
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    from gufe import LigandAtomMapping
    from netCDF4 import Dataset
    import openmm

    source = json.loads(args.source.read_text())
    if source.get("status") != "COMPLETE" or source.get("labels_used") is not False:
        raise ValueError("complete label-free source required")
    parent_path = Path(source["source_status"])
    if digest(parent_path) != source["source_status_sha256"]:
        raise ValueError("source status changed")
    parent = json.loads(parent_path.read_text())
    lookup = {tuple(row[field] for field in ("repeat", "iteration", "leg", "endpoint")): row for row in parent["results"]}
    if len(source["rows"]) != source["expected_points"] or len(lookup) != len(source["rows"]):
        raise ValueError("source frame inventory mismatch")
    contexts = {}
    rows = []
    maximum_closure = 0
    try:
        for item in source["rows"]:
            key = tuple(item[field] for field in ("repeat", "iteration", "leg", "endpoint"))
            directory = Path(lookup[key]["qm"]["log"]).parent
            manifest = json.loads((directory / "input_manifest.json").read_text())
            topology = Path(manifest["topology"])
            cache_key = str(topology), item["endpoint"]
            if cache_key not in contexts:
                mapping = LigandAtomMapping.from_json(topology.parents[2] / "mapping.json")
                maps, molecules = physical_atom_maps(read_topology(topology), mapping)
                molecule = molecules[item["endpoint"]]
                database = topology.with_name("A_db.json" if item["endpoint"] == 0 else "B_db.json")
                context, integrator, groups, metadata = parameterized_context(database, molecule)
                hybrid = openmm.XmlSerializer.deserialize(bz2.decompress(topology.with_name("hybrid_system.xml.bz2").read_bytes()).decode())
                with Dataset(manifest["trajectory"]) as trajectory:
                    saved = np.asarray(trajectory["analysis_particle_indices"][:], dtype=int)
                constrained = physical_constraint_pairs(hybrid, saved, maps[item["endpoint"]])
                contexts[cache_key] = context, integrator, groups, metadata, molecule, constrained
            context, _, groups, metadata, molecule, constrained = contexts[cache_key]
            with np.load(directory / "physical_coordinates.npz") as data:
                coordinates = data["ligand_xyz"]
            components = evaluate(context, groups, coordinates)
            if abs(components["total"] - item["vacuum_mm_kcal_mol"]) > 1e-6:
                raise ValueError("MM reference energy changed")
            terms = bonded_observables(context.getSystem(), molecule, coordinates, constrained)
            for kind, label in (("bond", "HarmonicBondForce"), ("angle", "HarmonicAngleForce")):
                closure = abs(sum(term["energy_kcal_mol"] for term in terms if term["kind"] == kind) - components[label])
                maximum_closure = max(maximum_closure, closure)
                if closure > 1e-6:
                    raise ValueError("bonded energy bookkeeping does not close")
            rows.append({**{field: item[field] for field in ("repeat", "iteration", "leg", "endpoint")},
                         "internal_correction_kcal_mol": item["internal_correction_kcal_mol"],
                         "mm_components": components, "terms": terms, "mm_metadata": metadata})
    finally:
        contexts.clear()
    report = {"status": "COMPLETE", "labels_used": False, "source_sha256": digest(args.source),
              "driver_sha256": digest(Path(__file__)), "frame_count": len(rows),
              "maximum_bonded_closure_kcal_mol": maximum_closure, "rows": rows, "states": summarize(rows)}
    args.report.parent.mkdir(parents=True, exist_ok=True)
    write_json(args.report.with_suffix(".json"), report)
    args.report.write_text(render(report))


if __name__ == "__main__":
    main()
