"""Diagnose frozen-frame QM/MM corrections using matched isolated ligand energies."""
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
import time

from report_quanpol_multiframe import describe, single_state_weight_ess, summarize
from run_matched_openfe import digest, write_json

HARTREE_KCAL_MOL = 627.5094740631


def isolated_deck(text):
    if len(re.findall(r"\$FFDATA\b", text, re.I)) != 1:
        raise ValueError("exactly one FFDATA group required")
    result = text[:text.index(" $FFDATA")]
    result, count = re.subn(r" \$QUANPO\b.*?\$END", "", result, flags=re.S | re.I)
    if count != 1 or "$DATA" not in result or "$DFTB" not in result:
        raise ValueError("expected DFTB QM/MM input groups")
    return result


def parse_vacuum(text):
    if "EXECUTION OF GAMESS TERMINATED NORMALLY" not in text:
        raise ValueError("isolated ligand GAMESS job did not terminate normally")
    if "SCF IS UNCONVERGED" in text:
        raise ValueError("isolated ligand SCF unconverged")
    matches = re.findall(r"FINAL RDFTB ENERGY IS\s+([-+\d.EeDd]+)", text)
    if len(matches) != 1:
        raise ValueError("one final isolated RDFTB energy required")
    energy = float(matches[0].replace("D", "E")) * HARTREE_KCAL_MOL
    if not math.isfinite(energy):
        raise ValueError("nonfinite isolated QM energy")
    return energy


def decompose(rows):
    by_key = {}
    for row in rows:
        key = (row["repeat"], row["iteration"], row["leg"], row["endpoint"])
        if key in by_key:
            raise ValueError("duplicate decomposition frame")
        components = [float(row[f"{component}_correction_kcal_mol"])
                      for component in ("total", "internal", "environment")]
        if not all(math.isfinite(value) for value in components):
            raise ValueError("nonfinite correction")
        if abs(components[0] - components[1] - components[2]) > 1e-7:
            raise ValueError("single-state decomposition does not close")
        by_key[key] = row
    cycles = []
    for repeat, iteration in sorted({key[:2] for key in by_key}):
        cycle = {"repeat": repeat, "iteration": iteration}
        for component in ("total", "internal", "environment"):
            values = {}
            for leg in ("bound", "solution"):
                for endpoint in (0, 10):
                    row = by_key.get((repeat, iteration, leg, endpoint))
                    if row is None:
                        raise ValueError("incomplete four-state decomposition")
                    value = float(row[f"{component}_correction_kcal_mol"])
                    if not math.isfinite(value):
                        raise ValueError("nonfinite correction")
                    values[leg, endpoint] = value
            cycle[component] = values["bound", 10] - values["bound", 0] - values["solution", 10] + values["solution", 0]
        cycle["closure"] = cycle["total"] - cycle["internal"] - cycle["environment"]
        if abs(cycle["closure"]) > 1e-7:
            raise ValueError("decomposition does not close")
        cycles.append(cycle)
    if not cycles:
        raise ValueError("no decomposition frames")
    by_state = {}
    for leg in ("bound", "solution"):
        for endpoint in (0, 10):
            components = {}
            for component in ("total", "internal", "environment"):
                values = [row[f"{component}_correction_kcal_mol"] for row in rows
                          if row["leg"] == leg and row["endpoint"] == endpoint]
                ess, _, _ = single_state_weight_ess(values)
                components[component] = {**describe(values), "weight_ess_uncorrected": ess}
            by_state[f"{leg}_{endpoint}"] = components
    internal = [cycle["internal"] for cycle in cycles]
    environment = [cycle["environment"] for cycle in cycles]
    covariance = statistics.covariance(internal, environment) if len(cycles) > 1 else None
    return {
        "cycles": cycles,
        "cycle_statistics": {component: describe(cycle[component] for cycle in cycles)
                             for component in ("total", "internal", "environment")},
        "internal_environment_sample_covariance": covariance,
        "max_abs_closure_kcal_mol": max(abs(cycle["closure"]) for cycle in cycles),
        "by_state": by_state,
    }


def render(report):
    lines = [f"# {report['target']} 局部 QM/MM 修正来源诊断", "",
             "所有原先固定的帧均保留；不优化几何、不读取实验亲和力、不按结果删帧。", "",
             "定义：内部修正 = 同坐标气相配体 QM − 同坐标气相配体 MM；环境修正 = 原 QM/MM−MM 修正 − 内部修正。",
             "环境修正含嵌入导致的 QM 电子响应等，不是独立去溶剂化自由能，也不是 tPIE。", "",
             "| 固定帧四态双差 (kcal/mol) | n | 均值 | SD |", "|---|---:|---:|---:|"]
    for component in ("total", "internal", "environment"):
        stats = report["cycle_statistics"][component]
        deviation = stats["sample_sd"]
        display = f"{deviation:.4f}" if deviation is not None else "NA"
        lines.append(f"| {component} | {stats['n']} | {stats['mean']:.4f} | {display} |")
    lines.extend(["", f"- 最大闭合误差：{report['max_abs_closure_kcal_mol']:.3g} kcal/mol。",
                  f"- 内部与环境双差的样本协方差：{report['internal_environment_sample_covariance']}。",
                  "- 两部分有协方差，不能用各自方差除以总方差当作独立因果贡献比例。", "",
                  "| 状态 | 总修正 SD / ESS | 内部修正 SD / ESS | 环境修正 SD / ESS |",
                  "|---|---:|---:|---:|"])
    for state, components in report["by_state"].items():
        cells = []
        for component in ("total", "internal", "environment"):
            stats = components[component]
            deviation = stats["sample_sd"]
            display = f"{deviation:.3f}" if deviation is not None else "NA"
            cells.append(f"{display} / {stats['weight_ess_uncorrected']:.2f}")
        lines.append(f"| {state} | " + " | ".join(cells) + " |")
    lines.extend(["", "## 解释边界", "", *[f"- {text}" for text in report["limitations"]], ""])
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--status", type=Path, required=True)
    parser.add_argument("--gamess", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    validation = summarize([args.status])
    if validation["qm_center_nonzero_single_points"]:
        raise ValueError("QM-center contaminated source cannot be decomposed")
    source = json.loads(args.status.read_text())
    args.output.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment.update(LD_LIBRARY_PATH=str(Path(sys.executable).resolve().parents[1] / "lib") + ":" + environment.get("LD_LIBRARY_PATH", ""),
                       OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1",
                       OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    from analyze_endpoint_internal_energy import evaluate, parameterized_context
    from gufe import LigandAtomMapping
    import numpy as np

    contexts = {}
    rows = []
    try:
        for original in source["results"]:
            directory = Path(original["qm"]["log"]).parent
            manifest = json.loads((directory / "input_manifest.json").read_text())
            topology = Path(manifest["topology"])
            endpoint = int(original["endpoint"])
            context_key = (str(topology), endpoint)
            if context_key not in contexts:
                mapping = LigandAtomMapping.from_json(topology.parents[2] / "mapping.json")
                component = mapping.componentA if endpoint == 0 else mapping.componentB
                contexts[context_key] = parameterized_context(topology.with_name("A_db.json" if endpoint == 0 else "B_db.json"), component.to_rdkit())
            context, _, groups, metadata = contexts[context_key]
            with np.load(directory / "physical_coordinates.npz") as coordinates:
                ligand_xyz = coordinates["ligand_xyz"]
            mm_energy = evaluate(context, groups, ligand_xyz)["total"]
            candidates = [path for path in directory.glob("qm_*.inp") if digest(path) == original["qm"]["deck_sha256"]]
            if len(candidates) != 1:
                raise ValueError("unique original QM deck with matching hash required")
            name = f"vac_r{original['repeat']}_i{original['iteration']}_{original['leg']}_{endpoint}"
            deck = args.output.resolve() / f"{name}.inp"
            log = deck.with_suffix(".stdout")
            deck.write_text(isolated_deck(candidates[0].read_text()))
            started = time.monotonic()
            with log.open("w") as handle:
                process = subprocess.run(["timeout", "180", "tcsh", "-f", str(args.gamess.resolve() / "rungms"), str(deck), "00", "2"],
                                         cwd=deck.parent, env=environment, stdout=handle, stderr=subprocess.STDOUT)
            if process.returncode:
                raise ValueError(f"vacuum QM process failed: {log}")
            qm_energy = parse_vacuum(log.read_text(errors="replace"))
            internal = qm_energy - mm_energy
            total = float(original["delta_u_kcal_mol"])
            rows.append({**{key: original[key] for key in ("repeat", "iteration", "leg", "endpoint")},
                         "total_correction_kcal_mol": total, "internal_correction_kcal_mol": internal,
                         "environment_correction_kcal_mol": total - internal,
                         "vacuum_qm_kcal_mol": qm_energy, "vacuum_mm_kcal_mol": mm_energy,
                         "deck_sha256": digest(deck), "source_deck_sha256": digest(candidates[0]),
                         "coordinates_sha256": digest(directory / "physical_coordinates.npz"),
                         "mm_metadata": metadata, "log": str(log), "wall_seconds": time.monotonic() - started})
            write_json(args.output / "status.json", {"status": "RUNNING", "expected_points": len(source["results"]), "rows": rows})
            print(f"{len(rows)}/{len(source['results'])}: {name}", flush=True)
    except BaseException as error:
        write_json(args.output / "status.json", {
            "status": "FAILED", "expected_points": len(source["results"]),
            "source_status_sha256": digest(args.status), "rows": rows,
            "error": f"{type(error).__name__}: {error}",
        })
        raise
    finally:
        contexts.clear()
    report = {"status": "COMPLETE", "target": source["target"], "labels_used": False,
              "source_status_sha256": digest(args.status), "source_status": str(args.status.resolve()),
              "driver_sha256": digest(Path(__file__)), "expected_points": len(source["results"]),
              **decompose(rows), "rows": rows,
              "limitations": [
                  "单条配体变换边、稀疏冻结帧的机制诊断，不是亲和力排序或自由能结果。",
                  "内部/环境拆分依赖这里的气相参照定义；不是唯一物理分解，也不是 PIEDA。",
                  "所有 ESS 仅为 298.15 K 单向权重诊断，未校正自相关；不能用环境部分 ESS 提高来证明总修正可重加权。",
                  "不能删除内部修正并称之为相同 QM/MM 模型；这样做会改变目标 Hamiltonian。",
                  "如使用 H_target = H_MM,PME + δ_local，需明确它是局部替换模型，不是完整周期 QM/MM，并单独验证重叠与收敛。",
              ]}
    write_json(args.output / "status.json", report)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    write_json(args.report.with_suffix(".json"), report)
    args.report.write_text(render(report))


if __name__ == "__main__":
    main()
