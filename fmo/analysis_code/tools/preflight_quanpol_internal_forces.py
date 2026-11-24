"""Check isolated DFTB3 gradients and internal force mismatch on four fixed states."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys

import numpy as np

from analyze_endpoint_internal_energy import parameterized_context
from decompose_quanpol_local_correction import HARTREE_KCAL_MOL, parse_vacuum
from diagnose_quanpol_internal_modes import replace_vacuum_coordinates
from run_matched_openfe import digest, write_json

BOHR_ANGSTROM = .529177210903


def parse_gradient(text, symbols):
    parse_vacuum(text)
    blocks = text.split("GRADIENT OF THE ENERGY")
    if len(blocks) != 2 or "UNITS ARE HARTREE/BOHR" not in blocks[1]:
        raise ValueError("one Hartree/bohr energy gradient required")
    entries = re.findall(r"^\s*(\d+)\s+([A-Za-z]+)\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s*$", blocks[1], flags=re.M)
    if len(entries) != len(symbols):
        raise ValueError("gradient atom inventory mismatch")
    for index, entry in enumerate(entries):
        if int(entry[0]) != index + 1 or entry[1].upper() != symbols[index].upper():
            raise ValueError("gradient atom order differs")
    gradient = np.array([[float(value.replace("D", "E")) for value in entry[2:]] for entry in entries])
    if not np.isfinite(gradient).all():
        raise ValueError("nonfinite gradient")
    return gradient * HARTREE_KCAL_MOL / BOHR_ANGSTROM


def tangent_projection(gradient, coordinates, constrained_pairs):
    jacobian = []
    for first, second in constrained_pairs:
        direction = coordinates[first] - coordinates[second]
        length = np.linalg.norm(direction)
        if length < 1e-8:
            raise ValueError("degenerate constraint bond")
        row = np.zeros_like(coordinates)
        row[first] = direction / length
        row[second] = -direction / length
        jacobian.append(row.reshape(-1))
    gradient = np.asarray(gradient, dtype=float)
    if not jacobian:
        return gradient.copy(), np.zeros_like(gradient)
    jacobian = np.asarray(jacobian)
    normal = jacobian.T @ np.linalg.pinv(jacobian @ jacobian.T) @ jacobian @ gradient.reshape(-1)
    return gradient - normal.reshape(gradient.shape), normal.reshape(gradient.shape)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--parent-status", type=Path, required=True)
    parser.add_argument("--bonded-report", type=Path, required=True)
    parser.add_argument("--gamess", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    from gufe import LigandAtomMapping
    from openmm import unit

    source = json.loads(args.source.read_text())
    bonded = json.loads(args.bonded_report.read_text())
    if source.get("status") not in ("RUNNING", "COMPLETE") or bonded.get("status") != "COMPLETE":
        raise ValueError("validated source results required")
    parent = json.loads(args.parent_status.read_text())
    if parent.get("status") != "COMPLETE" or parent.get("labels_used") is not False:
        raise ValueError("complete label-free parent required")
    parent_lookup = {tuple(row[field] for field in ("repeat", "iteration", "leg", "endpoint")): row for row in parent["results"]}
    earliest_repeat = min(row["repeat"] for row in source["rows"])
    earliest_time = min(row["iteration"] for row in source["rows"])
    selected = [row for row in source["rows"] if row["repeat"] == earliest_repeat and row["iteration"] == earliest_time]
    if {(row["leg"], row["endpoint"]) for row in selected} != {(leg, endpoint) for leg in ("bound", "solution") for endpoint in (0, 10)}:
        raise ValueError("four earliest fixed states required")
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment.update(LD_LIBRARY_PATH=str(Path(sys.executable).resolve().parents[1] / "lib") + ":" + environment.get("LD_LIBRARY_PATH", ""),
                       OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1",
                       OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")

    def execute(deck):
        log = deck.with_suffix(".stdout")
        with log.open("w") as handle:
            process = subprocess.run(["timeout", "180", "tcsh", "-f", str(args.gamess.resolve() / "rungms"), str(deck), "00", "2"],
                                     cwd=args.output, env=environment, stdout=handle, stderr=subprocess.STDOUT)
        if process.returncode:
            raise ValueError(f"GAMESS launch failed: {log}")
        text = log.read_text(errors="replace")
        return parse_vacuum(text), text, str(log)

    rows = []
    for index, item in enumerate(selected):
        state = f"{item['leg']}_{item['endpoint']}"
        original = next(row for row in bonded["rows"] if all(row[field] == item[field] for field in ("repeat", "iteration", "leg", "endpoint")))
        baseline_deck = Path(item["baseline_log"]).with_suffix(".inp")
        if digest(baseline_deck) != item["baseline_deck_sha256"]:
            raise ValueError("baseline deck changed")
        text = baseline_deck.read_text()
        block = text.split(" $DATA\n", 1)[1].split(" $END", 1)[0].splitlines()[2:]
        symbols = [line.split()[0] for line in block if line.strip()]
        coordinates = np.array([[float(value) for value in line.split()[2:]] for line in block if line.strip()])
        parent_row = parent_lookup[tuple(item[field] for field in ("repeat", "iteration", "leg", "endpoint"))]
        directory = Path(parent_row["qm"]["log"]).parent
        manifest = json.loads((directory / "input_manifest.json").read_text())
        topology = Path(manifest["topology"])
        mapping = LigandAtomMapping.from_json(topology.parents[2] / "mapping.json")
        molecule = (mapping.componentA if item["endpoint"] == 0 else mapping.componentB).to_rdkit()
        context, integrator, groups, _ = parameterized_context(topology.with_name("A_db.json" if item["endpoint"] == 0 else "B_db.json"), molecule)
        try:
            context.setPositions(coordinates * .1 * unit.nanometer)
            force = context.getState(getForces=True, groups=sum(1 << group for group in groups)).getForces(asNumpy=True).value_in_unit(unit.kilocalorie_per_mole / unit.angstrom)
            mm_gradient = -np.asarray(force)
        finally:
            del context, integrator
        deck = args.output / f"internal_grad_{state}.inp"
        deck.write_text(text.replace("RUNTYP=ENERGY", "RUNTYP=GRADIENT"))
        energy, output, log = execute(deck)
        if abs(energy - item["baseline_energy_kcal_mol"]) > 1e-6:
            raise ValueError("gradient and energy jobs disagree")
        qm_gradient = parse_gradient(output, symbols)
        difference = qm_gradient - mm_gradient
        constrained = [term["atoms"] for term in original["terms"] if term["kind"] == "bond" and term["constrained_in_source"]]
        tangent, normal = tangent_projection(difference, coordinates, constrained)
        finite_difference = None
        if index == 0:
            step = .0002
            energies = []
            for sign in (-1, 1):
                shifted = coordinates.copy()
                shifted[0, 0] += sign * step
                finite_deck = args.output / f"internal_fd_{'plus' if sign > 0 else 'minus'}.inp"
                finite_deck.write_text(replace_vacuum_coordinates(text, shifted))
                energies.append(execute(finite_deck)[0])
            numerical = (energies[1] - energies[0]) / (2 * step)
            error = abs(numerical - qm_gradient[0, 0])
            finite_difference = {"atom_index": 0, "axis": "x", "step_angstrom": step,
                                 "analytic_kcal_mol_angstrom": float(qm_gradient[0, 0]),
                                 "numerical_kcal_mol_angstrom": numerical, "absolute_error": float(error)}
            if error > .01:
                raise ValueError("analytic gradient fails finite difference check")
        rows.append({"state": state, "source_frame": {field: item[field] for field in ("repeat", "iteration")},
                     "gradient_difference_norm_kcal_mol_angstrom": float(np.linalg.norm(difference)),
                     "constraint_tangent_difference_norm_kcal_mol_angstrom": float(np.linalg.norm(tangent)),
                     "constraint_normal_difference_norm_kcal_mol_angstrom": float(np.linalg.norm(normal)),
                     "qm_gradient": qm_gradient.tolist(), "mm_gradient": mm_gradient.tolist(),
                     "constraint_count": len(constrained), "finite_difference": finite_difference,
                     "deck_sha256": digest(deck), "log": log})
        print(f"{len(rows)}/4: {state}", flush=True)
    report = {"status": "COMPLETE_ISOLATED_FORCE_PREFLIGHT_ONLY", "labels_used": False,
              "selection": "earliest_repeat_and_earliest_frozen_iteration_all_four_states",
              "driver_sha256": digest(Path(__file__)), "bonded_report_sha256": digest(args.bonded_report),
              "parent_status_sha256": digest(args.parent_status),
              "source_sha256": digest(args.source), "rows": rows,
              "full_qmmm_md_ready": False,
              "limitations": ["Only isolated internal gradients were checked, not full embedded periodic QM/MM forces.",
                              "Projection is Euclidean geometry diagnostic, not a mass-weighted integrator or force implementation.",
                              "Four fixed configurations cannot establish ensemble overlap or ranking accuracy."]}
    write_json(args.output / "status.json", report)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    write_json(args.report.with_suffix(".json"), report)
    lines = ["# CDK8 孤立配体 QM−MM 梯度前置核查", "",
             "事先按最早重复/时点取四个状态；不是按能量异常或实验标签选帧。单位：kcal/mol/Å。", "",
             "| 状态 | 梯度差范数 | 约束切向部分 | 约束法向部分 |", "|---|---:|---:|---:|"]
    for row in rows:
        lines.append(f"| {row['state']} | {row['gradient_difference_norm_kcal_mol_angstrom']:.3f} | {row['constraint_tangent_difference_norm_kcal_mol_angstrom']:.3f} | {row['constraint_normal_difference_norm_kcal_mol_angstrom']:.3f} |")
    check = rows[0]["finite_difference"]
    lines.extend(["", f"第一个固定构型的 atom 0 x 方向解析/中心差分误差为 {check['absolute_error']:.6g} kcal/mol/Å（步长 ±{check['step_angstrom']} Å）。仅是一个坐标方向的符号/单位局部核查，不是全维梯度验证。",
                  "约束投影使用实际轨迹粒子映射下的 H-bond 约束；这里只是欧氏几何诊断，不是质量加权约束动力学。",
                  "切向梯度差非零表明模型差异不全是约束键长法向上的差异，但不能据四个构型证明其系综后果或把范数当作结合能。",
                  "本轮仅核对孤立配体内部力，未验证嵌入环境/周期目标的完整力，尚不能启动受验证的 QM/MM 桥接 MD。", ""])
    args.report.write_text("\n".join(lines))


if __name__ == "__main__":
    main()
