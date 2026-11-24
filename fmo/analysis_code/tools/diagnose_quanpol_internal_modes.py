"""Run label-free hydrogen-only geometry controls on frozen QM/MM diagnostic frames."""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import re
import statistics
import subprocess
import sys
import tempfile
import time

import numpy as np

from analyze_endpoint_internal_energy import evaluate, parameterized_context
from decompose_quanpol_local_correction import parse_vacuum
from report_quanpol_multiframe import describe, single_state_weight_ess
from run_matched_openfe import digest, write_json


def rebind_parameters(text, alias):
    parameters = {}

    def substitute(match):
        first, second, old_path = match.groups()
        source = alias / Path(old_path).name
        if not source.is_file():
            raise ValueError(f"missing DFTB parameter: {source}")
        if source.stem.upper() != f"{first}-{second}".upper():
            raise ValueError("DFTB species/parameter filename mismatch")
        parameters[source.name] = digest(source)
        return f' {first} {second} "{source}"'

    updated, count = re.subn(r'^\s*(\w+)\s+(\w+)\s+"([^"\n]+\.skf)"\s*$', substitute, text, flags=re.M)
    if not count:
        raise ValueError("no DFTB parameter paths found")
    return updated, parameters


def replace_vacuum_coordinates(text, coordinates):
    lines = text.splitlines()
    starts = [index for index, line in enumerate(lines) if line.strip().upper() == "$DATA"]
    if len(starts) != 1 or "$FFDATA" in text or "$QUANPO" in text:
        raise ValueError("one isolated ligand DATA group required")
    first = starts[0] + 3
    last = next(index for index in range(first, len(lines)) if lines[index].strip().upper() == "$END")
    if last - first != len(coordinates):
        raise ValueError("ligand atom inventory mismatch")
    for offset, position in enumerate(coordinates):
        if len(position) != 3 or not np.isfinite(position).all():
            raise ValueError("nonfinite or invalid ligand coordinates")
        fields = lines[first + offset].split()
        if len(fields) != 5:
            raise ValueError("invalid ligand coordinate record")
        lines[first + offset] = f" {fields[0]} {fields[1]} " + " ".join(f"{value:.10f}" for value in position)
    return "\n".join(lines) + "\n"


def minimize_hydrogens(context, groups, coordinates, hydrogen_indices):
    from openmm import unit
    from scipy.optimize import minimize

    coordinates = np.asarray(coordinates, dtype=float)
    indices = np.asarray(hydrogen_indices, dtype=int)
    if not len(indices) or len(set(indices)) != len(indices) or min(indices) < 0 or max(indices) >= len(coordinates):
        raise ValueError("unique valid hydrogen indices required")
    mask = sum(1 << group for group in groups)

    def objective(values):
        current = coordinates.copy()
        current[indices] = values.reshape((-1, 3))
        context.setPositions(current * .1 * unit.nanometer)
        state = context.getState(getEnergy=True, getForces=True, groups=mask)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        forces = state.getForces(asNumpy=True).value_in_unit(unit.kilocalorie_per_mole / unit.angstrom)
        gradient = -np.asarray(forces)[indices].reshape(-1)
        if not math.isfinite(energy) or not np.isfinite(gradient).all():
            raise ValueError("nonfinite energy or gradient during hydrogen minimization")
        return energy, gradient

    before, _ = objective(coordinates[indices].reshape(-1))
    result = minimize(objective, coordinates[indices].reshape(-1), method="L-BFGS-B", jac=True,
                      options={"gtol": .005, "ftol": 1e-13, "maxiter": 1000, "maxls": 40})
    after, gradient = objective(result.x)
    if not result.success or max(abs(gradient)) > .02 or after > before + 1e-7:
        raise ValueError(f"hydrogen minimization failed validation: {result.message}; gradient={max(abs(gradient))}")
    final = coordinates.copy()
    final[indices] = result.x.reshape((-1, 3))
    return final, {"iterations": int(result.nit), "function_evaluations": int(result.nfev),
                   "max_abs_hydrogen_gradient_kcal_mol_angstrom": float(max(abs(gradient))),
                   "mm_energy_change_kcal_mol": float(after - before),
                   "maximum_hydrogen_displacement_angstrom": float(np.linalg.norm(final[indices] - coordinates[indices], axis=1).max()),
                   "heavy_atoms_displacement_angstrom": 0.0,
                   "force_mask_excludes_duplicate_lj_group31": 31 not in groups}


def summarize(rows):
    lookup = {(row["repeat"], row["iteration"], row["leg"], row["endpoint"]): row for row in rows}
    if len(lookup) != len(rows) or not rows:
        raise ValueError("nonempty unique frame matrix required")
    states = {}
    for leg in ("bound", "solution"):
        for endpoint in (0, 10):
            selected = [row for row in rows if row["leg"] == leg and row["endpoint"] == endpoint]
            if not selected:
                raise ValueError("missing thermodynamic state")
            controls = {}
            for control in ("original", "hydrogen_relaxed"):
                values = [row[control]["internal_correction_kcal_mol"] for row in selected]
                if not all(math.isfinite(value) for value in values):
                    raise ValueError("nonfinite internal correction")
                ess, _, _ = single_state_weight_ess(values)
                controls[control] = {**describe(values), "weight_ess_uncorrected": ess}
            states[f"{leg}_{endpoint}"] = controls
    cycles = []
    for repeat, iteration in sorted({key[:2] for key in lookup}):
        row = {"repeat": repeat, "iteration": iteration}
        for control in ("original", "hydrogen_relaxed"):
            values = {}
            for leg in ("bound", "solution"):
                for endpoint in (0, 10):
                    entry = lookup.get((repeat, iteration, leg, endpoint))
                    if entry is None:
                        raise ValueError("incomplete four-state matrix")
                    values[leg, endpoint] = entry[control]["internal_correction_kcal_mol"]
            row[control] = values["bound", 10] - values["bound", 0] - values["solution", 10] + values["solution", 0]
        cycles.append(row)
    return {"by_state": states, "cycles": cycles,
            "cycle_statistics": {control: describe(row[control] for row in cycles)
                                 for control in ("original", "hydrogen_relaxed")},
            "max_heavy_atom_displacement_angstrom": max(row["minimization"]["heavy_atoms_displacement_angstrom"] for row in rows),
            "max_hydrogen_displacement_angstrom": max(row["minimization"]["maximum_hydrogen_displacement_angstrom"] for row in rows)}


def render(report):
    lines = ["# CDK8 配体内部模式：固定重原子、仅氢松弛对照", "",
             "同一批 60 个冻结帧；重原子严格固定，仅在孤立配体 MM 势能上松弛氢，然后以相同 DFTB3 模型重新做孤立配体 QM 单点。",
             "这不是溶液态重采样，不是 QM 氢优化，也不是可直接替换原轨迹的自由能修正。", "",
             "| 四态内部能量修正双差 | 帧数 | 均值 (kcal/mol) | SD (kcal/mol) |", "|---|---:|---:|---:|"]
    for control, values in report["cycle_statistics"].items():
        deviation = values["sample_sd"]
        display = "NA" if deviation is None else f"{deviation:.4f}"
        lines.append(f"| {control} | {values['n']} | {values['mean']:.4f} | {display} |")
    lines.extend(["", "| 单态内部修正 | 原始 SD / ESS | 仅氢松弛 SD / ESS |", "|---|---:|---:|"])
    for state, controls in report["by_state"].items():
        cells = []
        for control in ("original", "hydrogen_relaxed"):
            values = controls[control]
            deviation = values["sample_sd"]
            display = "NA" if deviation is None else f"{deviation:.3f}"
            cells.append(f"{display} / {values['weight_ess_uncorrected']:.2f}")
        lines.append(f"| {state} | " + " | ".join(cells) + " |")
    lines.extend(["", f"- 最大重原子位移：{report['max_heavy_atom_displacement_angstrom']:.6g} Å。",
                  f"- 最大氢原子位移：{report['max_hydrogen_displacement_angstrom']:.6g} Å。", "",
                  "## 方法边界", "",
                  "- 内部修正仍为同坐标 QM−MM 差；不拿新几何的内部修正去拼接旧几何的环境修正。",
                  "- 松弛是确定性的几何映射，会改变参考分布；没有对应统计修正，不能将映射后的 ESS 当作原目标自由能收敛证据。",
                  "- MM 力常数/键角与能量差的相关仅作描述；因为 MM 能量本身就在差值中，相关不等于因果。",
                  "- 不读取实验标签、不删帧、不改生产 tPIE；本轮不生成全队列排序指标。", ""])
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--gamess", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.source.read_text())
    if source.get("status") != "COMPLETE" or source.get("labels_used") is not False:
        raise ValueError("complete label-free source required")
    parent_path = Path(source["source_status"])
    if digest(parent_path) != source["source_status_sha256"]:
        raise ValueError("parent status hash changed")
    parent = json.loads(parent_path.read_text())
    parent_rows = {(row["repeat"], row["iteration"], row["leg"], row["endpoint"]): row for row in parent["results"]}
    if len(parent_rows) != len(source["rows"]) or len(source["rows"]) != source["expected_points"]:
        raise ValueError("source coverage mismatch")
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    parameter_alias = Path(tempfile.mkdtemp(prefix="qph-")) / "sk"
    parameter_alias.symlink_to((args.gamess / "auxdata/DFTB/3OB-3-1").resolve(), target_is_directory=True)
    protocol = {"status": "RUNNING", "schema": "quanpol-internal-hydrogen-control-v1", "labels_used": False,
                "expected_points": len(source["rows"]), "source_sha256": digest(args.source),
                "driver_sha256": digest(Path(__file__)), "selection": "ALL_FROZEN_SOURCE_FRAMES",
                "operation": "isolated_MM_hydrogen_only_minimization_fixed_heavy_atoms",
                "target_unchanged": False, "production_pipeline_modified": False,
                "parameter_alias": str(parameter_alias),
                "gamess_binary_sha256": digest(args.gamess / "gamess.00.x")}
    write_json(args.output / "protocol.json", protocol)
    from gufe import LigandAtomMapping
    from openmm import unit

    environment = os.environ.copy()
    environment.update(LD_LIBRARY_PATH=str(Path(sys.executable).resolve().parents[1] / "lib") + ":" + environment.get("LD_LIBRARY_PATH", ""),
                       OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1",
                       OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    rows = []
    contexts = {}
    try:
        for item in source["rows"]:
            key = tuple(item[field] for field in ("repeat", "iteration", "leg", "endpoint"))
            original = parent_rows[key]
            directory = Path(original["qm"]["log"]).parent
            manifest = json.loads((directory / "input_manifest.json").read_text())
            topology = Path(manifest["topology"])
            context_key = str(topology), int(item["endpoint"])
            if context_key not in contexts:
                mapping = LigandAtomMapping.from_json(topology.parents[2] / "mapping.json")
                component = mapping.componentA if item["endpoint"] == 0 else mapping.componentB
                molecule = component.to_rdkit()
                database = topology.with_name("A_db.json" if item["endpoint"] == 0 else "B_db.json")
                contexts[context_key] = (*parameterized_context(database, molecule), molecule)
            context, _, groups, metadata, molecule = contexts[context_key]
            with np.load(directory / "physical_coordinates.npz") as data:
                coordinates = data["ligand_xyz"]
            before = evaluate(context, groups, coordinates)
            if abs(before["total"] - item["vacuum_mm_kcal_mol"]) > 1e-6:
                raise ValueError("original cached ligand MM energy does not reproduce")
            hydrogen = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() == 1]
            heavy = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() != 1]
            relaxed, minimization = minimize_hydrogens(context, groups, coordinates, hydrogen)
            if not np.array_equal(relaxed[heavy], coordinates[heavy]):
                raise ValueError("fixed heavy atom moved")
            after = evaluate(context, groups, relaxed)
            mask = sum(1 << group for group in groups)
            context.setPositions(coordinates * .1 * unit.nanometer)
            forces = context.getState(getForces=True, groups=mask).getForces(asNumpy=True).value_in_unit(unit.kilocalorie_per_mole / unit.angstrom)
            name = f"hcontrol_r{item['repeat']}_i{item['iteration']}_{item['leg']}_{item['endpoint']}"
            deck = args.output / f"{name}.inp"
            log = deck.with_suffix(".stdout")
            vacuum_source = Path(item["log"]).with_suffix(".inp")
            if digest(vacuum_source) != item["deck_sha256"]:
                raise ValueError("original vacuum deck hash changed")
            rebound, parameter_hashes = rebind_parameters(vacuum_source.read_text(), parameter_alias)
            baseline_deck = args.output / f"baseline_{name}.inp"
            baseline_log = baseline_deck.with_suffix(".stdout")
            baseline_deck.write_text(rebound)
            baseline_started = time.monotonic()
            with baseline_log.open("w") as handle:
                baseline_process = subprocess.run(["timeout", "180", "tcsh", "-f", str(args.gamess.resolve() / "rungms"), str(baseline_deck), "00", "2"],
                                                  cwd=args.output, env=environment, stdout=handle, stderr=subprocess.STDOUT)
            if baseline_process.returncode:
                raise ValueError(f"baseline GAMESS process failed: {baseline_log}")
            baseline_energy = parse_vacuum(baseline_log.read_text(errors="replace"))
            if abs(baseline_energy - item["vacuum_qm_kcal_mol"]) > 1e-6:
                raise ValueError("historical QM energy does not reproduce in current runtime")
            baseline_wall = time.monotonic() - baseline_started
            deck.write_text(replace_vacuum_coordinates(rebound, relaxed))
            started = time.monotonic()
            with log.open("w") as handle:
                process = subprocess.run(["timeout", "180", "tcsh", "-f", str(args.gamess.resolve() / "rungms"), str(deck), "00", "2"],
                                         cwd=args.output, env=environment, stdout=handle, stderr=subprocess.STDOUT)
            if process.returncode:
                raise ValueError(f"GAMESS process failed: {log}")
            qm_energy = parse_vacuum(log.read_text(errors="replace"))
            np.savez(args.output / f"{name}.npz", original=coordinates, hydrogen_relaxed=relaxed,
                     hydrogen_indices=hydrogen, heavy_indices=heavy)
            rows.append({**{field: item[field] for field in ("repeat", "iteration", "leg", "endpoint")},
                         "original": {"mm": before, "qm_kcal_mol": item["vacuum_qm_kcal_mol"],
                                      "internal_correction_kcal_mol": item["internal_correction_kcal_mol"]},
                         "hydrogen_relaxed": {"mm": after, "qm_kcal_mol": qm_energy,
                                              "internal_correction_kcal_mol": qm_energy - after["total"]},
                         "mm_hydrogen_force_rms_kcal_mol_angstrom": float(np.sqrt(np.mean(np.asarray(forces)[hydrogen] ** 2))),
                         "mm_heavy_force_rms_kcal_mol_angstrom": float(np.sqrt(np.mean(np.asarray(forces)[heavy] ** 2))),
                         "minimization": minimization, "mm_metadata": metadata,
                         "baseline_log": str(baseline_log), "baseline_energy_kcal_mol": baseline_energy,
                         "baseline_deck_sha256": digest(baseline_deck), "baseline_wall_seconds": baseline_wall,
                         "parameter_sha256": parameter_hashes,
                         "log": str(log), "deck_sha256": digest(deck), "wall_seconds": time.monotonic() - started})
            write_json(args.output / "status.json", {**protocol, "rows": rows})
            print(f"{len(rows)}/{len(source['rows'])}: {name}", flush=True)
    except BaseException as error:
        write_json(args.output / "status.json", {**protocol, "status": "FAILED", "rows": rows,
                                                "error": f"{type(error).__name__}: {error}"})
        raise
    finally:
        contexts.clear()
    report = {**protocol, "status": "COMPLETE", "rows": rows, **summarize(rows)}
    write_json(args.output / "status.json", report)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    write_json(args.report.with_suffix(".json"), report)
    args.report.write_text(render(report))


if __name__ == "__main__":
    main()
