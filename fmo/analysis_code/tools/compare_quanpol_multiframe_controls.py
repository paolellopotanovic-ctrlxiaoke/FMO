"""Compare matched finite-shell runs with and without the QM-center penalty."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import statistics

from quanpol_constants import QM_CENTER_RESTRAINT_THRESHOLD_ANGSTROM


def read_cycles(path):
    report = json.loads(path.read_text())
    if report.get("status") != "COMPLETE" or report.get("labels_used"):
        raise ValueError(f"input must be a complete label-independent run: {path}")
    rows = {}
    for row in report["results"]:
        key = (int(row["repeat"]), int(row["iteration"]), row["leg"], int(row["endpoint"]))
        if key in rows:
            raise ValueError(f"duplicate result {key} in {path}")
        rows[key] = row
    repeats = sorted({key[0] for key in rows})
    iterations = sorted({key[1] for key in rows})
    expected = {(repeat, iteration, leg, endpoint) for repeat in repeats for iteration in iterations
                for leg in ("bound", "solution") for endpoint in (0, 10)}
    if set(rows) != expected:
        raise ValueError(f"incomplete paired matrix in {path}")
    for key, row in rows.items():
        for model in ("qm", "mm"):
            result = row[model]
            if result.get("status") != "SUCCEEDED" or result.get("scf_unconverged"):
                raise ValueError(f"failed or unconverged {model} state {key} in {path}")
    return report, rows


def build_comparison(with_center_path, without_center_path):
    with_center, original = read_cycles(with_center_path)
    without_center, corrected = read_cycles(without_center_path)
    if with_center.get("target") != without_center.get("target"):
        raise ValueError("runs must use the same target")
    if float(with_center.get("swrb_angstrom", 1001.0)) == float(without_center.get("swrb_angstrom", 1001.0)):
        raise ValueError("comparison runs must use distinct SWRB settings")
    if set(original) != set(corrected):
        raise ValueError("runs must use identical repeats, frames, and states")

    def sha256(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    def matched_artifacts(row):
        frame_dir = Path(row["qm"]["log"]).parent
        stem = (f"r{row['repeat']}_i{row['iteration']}_"
                f"{row['leg']}_{row['endpoint']}")
        coordinate_path = frame_dir / "physical_coordinates.npz"
        environment_path = frame_dir / "environment.pdb"
        qm_deck_path = frame_dir / f"qm_{stem}.inp"
        mm_deck_path = frame_dir / f"mm_{stem}.inp"
        required = (coordinate_path, environment_path, qm_deck_path, mm_deck_path)
        if not all(path.is_file() for path in required):
            raise ValueError(f"missing frozen comparison artifact for {stem}")
        qm_text = qm_deck_path.read_text()
        qm_text = re.sub(r"/tmp/qp-endpoint-[^/]+", "/tmp/qp-endpoint-<TMP>", qm_text)
        qm_text, swrb_count = re.subn(r"\bSWRB\s*=\s*[-+0-9.Ee]+", "SWRB=<VALUE>", qm_text)
        if swrb_count != 1:
            raise ValueError(f"expected exactly one SWRB field in {qm_deck_path}")
        coordinate_block = qm_text.split("COORDINATES\n", 1)[1].split("\nSTOP", 1)[0]
        coordinates = []
        for line in coordinate_block.splitlines():
            fields = line.split()
            if len(fields) == 5:
                try:
                    coordinates.append([float(value) for value in fields[-3:]])
                except ValueError:
                    continue
        if not coordinates:
            raise ValueError(f"missing FFDATA coordinates in {qm_deck_path}")
        axes = list(zip(*coordinates, strict=True))
        coordinate_box_diagonal_bound = sum(
            (max(axis) - min(axis)) ** 2 for axis in axes
        ) ** 0.5
        return {
            "coordinates": sha256(coordinate_path),
            "environment": sha256(environment_path),
            "qm_deck_without_swrb": qm_text,
            "mm_deck": sha256(mm_deck_path),
            "trajectory": row["trajectory"],
            "replica": row["replica"],
            "coordinate_box_diagonal_bound_angstrom": coordinate_box_diagonal_bound,
        }

    max_coordinate_box_diagonal_bound = 0.0
    for key in original:
        before = matched_artifacts(original[key])
        after = matched_artifacts(corrected[key])
        if before != after:
            raise ValueError(f"comparison changed coordinates or non-SWRB input for {key}")
        max_coordinate_box_diagonal_bound = max(
            max_coordinate_box_diagonal_bound,
            before["coordinate_box_diagonal_bound_angstrom"],
        )

    with_center_swrb = float(with_center.get("swrb_angstrom", 1001.0))
    without_center_swrb = float(without_center["swrb_angstrom"])
    if max_coordinate_box_diagonal_bound >= min(with_center_swrb, without_center_swrb):
        raise ValueError("a finite-shell pair may be cut off by one of the compared SWRB values")

    def cycle(rows, repeat, iteration, field):
        def value(leg, endpoint):
            return float(rows[(repeat, iteration, leg, endpoint)][field])
        return value("bound", 10) - value("bound", 0) - value("solution", 10) + value("solution", 0)

    def qm_center(row):
        path = Path(row["qm"]["log"])
        matches = re.findall(r"QM CENTER\s+ENERGY =\s*([-+0-9.Ee]+)\s+KCAL/MOL",
                             path.read_text(errors="replace"))
        if not matches:
            raise ValueError(f"missing QM-center kcal/mol record in {path}")
        return float(matches[-1])

    repeats = sorted({key[0] for key in original})
    iterations = sorted({key[1] for key in original})
    rows = []
    for repeat in repeats:
        for iteration in iterations:
            old_cycle = cycle(original, repeat, iteration, "delta_u_kcal_mol")
            new_cycle = cycle(corrected, repeat, iteration, "delta_u_kcal_mol")
            center_cycle = (qm_center(original[(repeat, iteration, "bound", 10)])
                            - qm_center(original[(repeat, iteration, "bound", 0)])
                            - qm_center(original[(repeat, iteration, "solution", 10)])
                            + qm_center(original[(repeat, iteration, "solution", 0)]))
            mm_differences = [
                float(original[key]["mm"]["potential_kcal_mol"])
                - float(corrected[key]["mm"]["potential_kcal_mol"])
                for key in original if key[:2] == (repeat, iteration)
            ]
            rows.append({"repeat": repeat, "iteration": iteration,
                         "with_center_kcal_mol": old_cycle,
                         "without_center_kcal_mol": new_cycle,
                         "qm_center_component_kcal_mol": center_cycle,
                         "observed_cycle_shift_kcal_mol": new_cycle - old_cycle,
                         "mm_max_abs_change_kcal_mol": max(abs(value) for value in mm_differences)})
    old_values = [row["with_center_kcal_mol"] for row in rows]
    new_values = [row["without_center_kcal_mol"] for row in rows]
    center_values = [row["qm_center_component_kcal_mol"] for row in rows]
    qm_single_point_closure_errors = []
    mm_single_point_changes = []
    openmm_single_point_changes = []
    for key in original:
        old_row = original[key]
        new_row = corrected[key]
        center = qm_center(old_row)
        qm_shift = (float(old_row["qm"]["potential_kcal_mol"])
                    - float(new_row["qm"]["potential_kcal_mol"]))
        qm_single_point_closure_errors.append(abs(qm_shift - center))
        mm_single_point_changes.append(
            float(old_row["mm"]["potential_kcal_mol"])
            - float(new_row["mm"]["potential_kcal_mol"])
        )
        openmm_single_point_changes.append(
            float(old_row["openmm_ffdata"]["potential_kcal_mol"])
            - float(new_row["openmm_ffdata"]["potential_kcal_mol"])
        )
    max_center_single_point = max(
        abs(qm_center(row)) for row in original.values()
    )
    max_without_center_single_point = max(
        abs(qm_center(row)) for row in corrected.values()
    )
    max_single_point_closure = max(qm_single_point_closure_errors)
    max_mm_change = max(abs(value) for value in mm_single_point_changes)
    max_openmm_change = max(abs(value) for value in openmm_single_point_changes)
    if max_single_point_closure > 1e-6:
        raise ValueError("QM single-point energy shift does not close against the QM-center term")
    if max_without_center_single_point > 1e-8:
        raise ValueError("the comparison run still contains a nonzero QM-center term")
    if max_mm_change > 1e-8 or max_openmm_change > 1e-8:
        raise ValueError("MM or OpenMM reference energies changed between controls")
    return {
        "schema": "quanpol-qm-center-factorial-comparison-v1",
        "status": "COMPLETE",
        "target": with_center["target"],
        "labels_used": False,
        "with_center_swrb_angstrom": with_center_swrb,
        "without_center_swrb_angstrom": without_center_swrb,
        "max_coordinate_box_diagonal_bound_angstrom": max_coordinate_box_diagonal_bound,
        "qm_center_restraint_branch_threshold_angstrom": QM_CENTER_RESTRAINT_THRESHOLD_ANGSTROM,
        "fixed_frame_count": len(rows),
        "with_center_cycle_mean_kcal_mol": statistics.mean(old_values),
        "without_center_cycle_mean_kcal_mol": statistics.mean(new_values),
        "with_center_cycle_sample_sd_kcal_mol": statistics.stdev(old_values) if len(old_values) > 1 else None,
        "without_center_cycle_sample_sd_kcal_mol": statistics.stdev(new_values) if len(new_values) > 1 else None,
        "qm_center_cycle_mean_kcal_mol": statistics.mean(center_values),
        "qm_center_cycle_max_abs_kcal_mol": max(abs(value) for value in center_values),
        "qm_center_single_point_max_abs_kcal_mol": max_center_single_point,
        "qm_center_nonzero_single_point_count": sum(
            abs(qm_center(row)) > 1e-8 for row in original.values()
        ),
        "without_center_single_point_max_abs_kcal_mol": max_without_center_single_point,
        "without_center_nonzero_single_point_count": sum(
            abs(qm_center(row)) > 1e-8 for row in corrected.values()
        ),
        "max_abs_closure_error_kcal_mol": max(
            abs(row["observed_cycle_shift_kcal_mol"] + row["qm_center_component_kcal_mol"])
            for row in rows
        ),
        "max_abs_single_point_closure_error_kcal_mol": max_single_point_closure,
        "mm_max_abs_change_kcal_mol": max(row["mm_max_abs_change_kcal_mol"] for row in rows),
        "mm_single_point_max_abs_change_kcal_mol": max_mm_change,
        "openmm_ffdata_single_point_max_abs_change_kcal_mol": max_openmm_change,
        "cycles": rows,
        "interpretation_limits": [
            "This is a single Hamiltonian-factor sensitivity control on sparse fixed frames, not an affinity estimate or free energy.",
            "The QM-center term is a GAMESS finite-shell algorithm restraint; removing it does not validate the remaining finite-shell model against PME.",
            "All selected frames are retained; no frame is deleted based on observed energy or affinity.",
        ],
    }


def render_markdown(report):
    lines = [
        f"# {report['target']} QM-center 罚能敏感性对照",
        "",
        "对完全相同的坐标、配体端态、FFDATA MM 参考与采样时点，仅改变 GAMESS `SWRB`，比较触发/关闭 QM-center 人工约束时的 QM-minus-MM 四态双差值。",
        "",
        f"- `SWRB={report['with_center_swrb_angstrom']:g} Å`（中心罚能可触发）：平均双差值 `{report['with_center_cycle_mean_kcal_mol']:.3f}`，SD `{report['with_center_cycle_sample_sd_kcal_mol']:.3f}` kcal/mol。",
        f"- `SWRB={report['without_center_swrb_angstrom']:g} Å`（源码条件关闭中心罚能）：平均双差值 `{report['without_center_cycle_mean_kcal_mol']:.3f}`，SD `{report['without_center_cycle_sample_sd_kcal_mol']:.3f}` kcal/mol。",
        f"- 中心项双差均值 `{report['qm_center_cycle_mean_kcal_mol']:.3f}`，最大绝对值 `{report['qm_center_cycle_max_abs_kcal_mol']:.3f}` kcal/mol。",
        f"- 旧设置下 QM-center 项非零的单点 `{report['qm_center_nonzero_single_point_count']}/{4 * report['fixed_frame_count']}`；修正设置下 `{report['without_center_nonzero_single_point_count']}/{4 * report['fixed_frame_count']}`。",
        f"- 全部坐标包围盒对角线上界 `{report['max_coordinate_box_diagonal_bound_angstrom']:.3f}` Å，小于较小的 `SWRB`，因此该几何对中两组有限壳层距离截断都未触及。",
        f"- GAMESS QM-center 分支门槛为 `SWRB² < 1e8 bohr²`（约 `{report['qm_center_restraint_branch_threshold_angstrom']:.2f}` Å）；中心最近 QM 原子距几何中心超过 2 Å 时才增加罚能。",
        f"- 逐帧变化与 QM-center 项的最大闭合误差 `{report['max_abs_closure_error_kcal_mol']:.3g}` kcal/mol；MM 参考最大变化 `{report['mm_max_abs_change_kcal_mol']:.3g}` kcal/mol。",
        f"- 单点 QM 能量变化与日志罚能的最大闭合误差 `{report['max_abs_single_point_closure_error_kcal_mol']:.3g}` kcal/mol；MM-only 与 OpenMM FFDATA 参考变化分别为 `{report['mm_single_point_max_abs_change_kcal_mol']:.3g}` / `{report['openmm_ffdata_single_point_max_abs_change_kcal_mol']:.3g}` kcal/mol。",
        "",
        "| repeat | iteration | 有中心罚能 | 无中心罚能 | QM-center 项 | 新旧差值 |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    lines.extend(
        f"| {row['repeat']} | {row['iteration']} | {row['with_center_kcal_mol']:.3f} | {row['without_center_kcal_mol']:.3f} | {row['qm_center_component_kcal_mol']:.3f} | {row['observed_cycle_shift_kcal_mol']:.3f} |"
        for row in report["cycles"]
    )
    lines.extend([
        "",
        "结论：若表中闭合良好且 MM 未变化，原双差值中的对应变化可直接归因于 GAMESS QM-center 项，而不是 tPIE、PME 或实验亲和力。该项是有限壳层 QM/MM 计算的算法性中心约束，不应未经验证混入结合自由能修正。关闭后结果仍是固定帧势能诊断，尚非 ΔG/FEP。",
        "",
    ])
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--with-center", type=Path, required=True)
    parser.add_argument("--without-center", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path, required=True)
    args = parser.parse_args()
    report = build_comparison(args.with_center, args.without_center)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(report, indent=2) + "\n")
    args.output_md.parent.mkdir(parents=True, exist_ok=True)
    args.output_md.write_text(render_markdown(report))
    print(json.dumps({key: report[key] for key in (
        "status", "fixed_frame_count", "with_center_cycle_mean_kcal_mol",
        "without_center_cycle_mean_kcal_mol", "qm_center_cycle_mean_kcal_mol",
        "max_abs_closure_error_kcal_mol")}, indent=2))


if __name__ == "__main__":
    main()
