"""Compare matched QM/MM corrections across explicit-solvent shell radii."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import statistics

import numpy as np


def load_run(path):
    report = json.loads(path.read_text())
    if report.get("status") != "COMPLETE" or report.get("labels_used"):
        raise ValueError(f"run must be complete and label-independent: {path}")
    if report.get("qm_center_restraint_allowed"):
        raise ValueError(f"QM-center restraint control is not a production shell run: {path}")
    rows = {}
    for row in report["results"]:
        key = (int(row["repeat"]), int(row["iteration"]), row["leg"], int(row["endpoint"]))
        if key in rows:
            raise ValueError(f"duplicate frame state {key} in {path}")
        for model in ("qm", "mm"):
            calculation = row[model]
            if calculation.get("status") != "SUCCEEDED" or calculation.get("scf_unconverged"):
                raise ValueError(f"failed or unconverged {model} result for {key} in {path}")
        if abs(float(row["mm_parity_difference_kcal_mol"])) > 0.1:
            raise ValueError(f"matched MM parity failed for {key} in {path}")
        rows[key] = row
    if not rows:
        raise ValueError(f"run contains no calculation rows: {path}")
    return report, rows


def read_frame_manifest(row):
    frame_dir = Path(row["qm"]["log"]).parent
    manifest_path = frame_dir / "input_manifest.json"
    coordinates_path = frame_dir / "physical_coordinates.npz"
    if not manifest_path.is_file() or not coordinates_path.is_file():
        raise ValueError(f"missing frozen input artifacts in {frame_dir}")
    manifest = json.loads(manifest_path.read_text())
    with np.load(coordinates_path) as coordinates:
        ligand_coordinates = np.array(coordinates["ligand_xyz"], copy=True)
        environment_coordinates = np.array(coordinates["environment_xyz"], copy=True)
        environment_indices = np.array(coordinates["environment_original_indices"], dtype=int, copy=True)
    solvent_names = {"HOH", "WAT", "NA", "CL"}
    protein_indices = {int(index) for entry in manifest["residues"]
                       if entry["name"].upper() not in solvent_names for index in entry["indices"]}
    protein_coordinates = {int(index): coordinate for index, coordinate in zip(
        environment_indices, environment_coordinates, strict=True) if int(index) in protein_indices}
    water_residues = sum(entry["name"].upper() in ("HOH", "WAT") for entry in manifest["residues"])
    ion_residues = sum(entry["name"].upper() in ("NA", "CL") for entry in manifest["residues"])
    return manifest, ligand_coordinates, protein_coordinates, water_residues, ion_residues


def summarize(named_paths):
    grouped = {}
    for shell, path in named_paths:
        report, rows = load_run(path)
        if not np.isclose(float(report.get("shell_angstrom", np.nan)), shell):
            raise ValueError(f"declared radius {shell:g} Å does not match run metadata in {path}")
        if float(report.get("swrb_angstrom", 1001.0)) <= 5291.77210903:
            raise ValueError(f"run may activate the GAMESS QM-center penalty: {path}")
        group = grouped.setdefault(shell, {"shell": shell, "paths": [], "report": report, "rows": {}})
        if (report.get("target") != group["report"].get("target")
                or float(report.get("swrb_angstrom")) != float(group["report"].get("swrb_angstrom"))):
            raise ValueError(f"partial runs at {shell:g} Å use incompatible settings")
        duplicate_keys = set(group["rows"]) & set(rows)
        if duplicate_keys:
            raise ValueError(f"overlapping rows for radius {shell:g} Å: {sorted(duplicate_keys)}")
        group["paths"].append(path)
        group["rows"].update(rows)
    runs = [grouped[shell] for shell in sorted(grouped)]
    if len(runs) < 2:
        raise ValueError("at least two distinct shell radii are required")

    target_set = {run["report"].get("target") for run in runs}
    key_sets = {tuple(sorted(run["rows"])) for run in runs}
    if len(target_set) != 1 or len(key_sets) != 1:
        raise ValueError("all shell runs must use the same target and complete frame/state matrix")
    reference_keys = set(runs[0]["rows"])
    keys = sorted(reference_keys)
    expected = {(repeat, iteration, leg, endpoint)
                for repeat, iteration in {(key[0], key[1]) for key in keys}
                for leg in ("bound", "solution") for endpoint in (0, 10)}
    if reference_keys != expected:
        raise ValueError("shell sensitivity requires complete four-state cycles")

    frame_records = []
    for key in keys:
        reference_row = runs[0]["rows"][key]
        reference_manifest, reference_ligand, reference_protein, _, _ = read_frame_manifest(reference_row)
        for run in runs[1:]:
            row = run["rows"][key]
            manifest, ligand_coordinates, protein_coordinates, _, _ = read_frame_manifest(row)
            if (row.get("trajectory") != reference_row.get("trajectory")
                    or row.get("replica") != reference_row.get("replica")
                    or manifest.get("trajectory") != reference_manifest.get("trajectory")
                    or manifest.get("replica") != reference_manifest.get("replica")
                    or not np.allclose(ligand_coordinates, reference_ligand, atol=1e-7, rtol=0)
                    or set(protein_coordinates) != set(reference_protein)
                    or any(not np.allclose(protein_coordinates[index], reference_protein[index],
                                           atol=1e-7, rtol=0) for index in reference_protein)):
                raise ValueError(f"shell comparison changed the sampled frame or QM ligand geometry: {key}")

        by_shell = []
        for run in runs:
            row = run["rows"][key]
            manifest, _, _, water_count, ion_count = read_frame_manifest(row)
            by_shell.append({"shell_angstrom": run["shell"],
                             "qm_atoms": int(manifest["qm_atoms"]),
                             "mm_atoms": int(manifest["mm_atoms"]),
                             "water_residues": water_count,
                             "ion_residues": ion_count,
                             "delta_u_kcal_mol": float(row["delta_u_kcal_mol"]),
                             "mm_parity_difference_kcal_mol": float(row["mm_parity_difference_kcal_mol"])})
        frame_records.append({"repeat": key[0], "iteration": key[1], "leg": key[2],
                              "endpoint": key[3], "trajectory": reference_row["trajectory"],
                              "replica": reference_row["replica"], "shell_results": by_shell})

    cycles = []
    for repeat, iteration in sorted({(key[0], key[1]) for key in keys}):
        shell_cycles = []
        for run in runs:
            rows = run["rows"]
            def delta(leg, endpoint):
                return float(rows[(repeat, iteration, leg, endpoint)]["delta_u_kcal_mol"])
            correction = delta("bound", 10) - delta("bound", 0) - delta("solution", 10) + delta("solution", 0)
            shell_cycles.append({"shell_angstrom": run["shell"], "correction_kcal_mol": correction})
        cycles.append({"repeat": repeat, "iteration": iteration, "shell_corrections": shell_cycles})

    summaries = {}
    for run in runs:
        values = [result["correction_kcal_mol"] for cycle in cycles
                  for result in cycle["shell_corrections"]
                  if result["shell_angstrom"] == run["shell"]]
        summaries[str(run["shell"])] = {
            "n": len(values), "mean_kcal_mol": statistics.mean(values),
            "sample_sd_kcal_mol": statistics.stdev(values) if len(values) > 1 else None,
            "median_kcal_mol": statistics.median(values),
            "min_kcal_mol": min(values), "max_kcal_mol": max(values),
        }
    reference_shell = 12.0 if 12.0 in grouped else runs[len(runs) // 2]["shell"]
    reference_values = {
        (cycle["repeat"], cycle["iteration"]): next(
            result["correction_kcal_mol"] for result in cycle["shell_corrections"]
            if result["shell_angstrom"] == reference_shell)
        for cycle in cycles
    }
    differences_vs_reference = {}
    for run in runs:
        if run["shell"] == reference_shell:
            continue
        differences = []
        for cycle in cycles:
            shell_value = next(result["correction_kcal_mol"] for result in cycle["shell_corrections"]
                               if result["shell_angstrom"] == run["shell"])
            differences.append(shell_value - reference_values[(cycle["repeat"], cycle["iteration"])])
        differences_vs_reference[str(run["shell"])] = {
            "reference_shell_angstrom": reference_shell, "n": len(differences),
            "mean_difference_kcal_mol": statistics.mean(differences),
            "sample_sd_kcal_mol": statistics.stdev(differences) if len(differences) > 1 else None,
            "min_kcal_mol": min(differences), "max_kcal_mol": max(differences),
        }

    return {"schema": "quanpol-shell-sensitivity-v1", "status": "COMPLETE",
            "target": next(iter(target_set)), "labels_used": False,
            "swrb_angstrom": [run["report"]["swrb_angstrom"] for run in runs],
            "shell_radii_angstrom": [run["shell"] for run in runs],
            "cycle_summary_by_shell_kcal_mol": summaries,
            "paired_difference_vs_reference_shell_kcal_mol": differences_vs_reference,
            "reference_shell_angstrom": reference_shell,
            "fixed_four_state_cycles": len(cycles), "frame_states": frame_records,
            "cycles": cycles,
            "interpretation_limits": [
                "This is a finite explicit-water/ion shell sensitivity test; protein residues and QM ligand coordinates are held fixed.",
                "It does not restore PME, periodic long-range electrostatics, mutual polarization, or missing conformational sampling.",
                "Sparse fixed-frame results are diagnostic only and cannot establish affinity ranking or a converged free-energy correction.",
            ]}


def render_markdown(report):
    lines = [
        f"# {report['target']} QM/MM 水/离子壳层敏感性",
        "",
        "同一条 OpenFE 轨迹的同一帧、相同配体 QM 坐标和固定蛋白残基；仅改变纳入的水/离子半径。所有计算使用关闭 QM-center 罚能的 GAMESS 设置。",
        "",
        f"- 固定四态双差数量：{report['fixed_four_state_cycles']}；壳层半径：`{report['shell_radii_angstrom']}` Å。",
        "- 蛋白残基全部保留；变动部分仅为按配体距离筛选的水/离子。每个壳层下仍要求匹配 MM 与 OpenMM FFDATA 能量误差小于 0.1 kcal/mol。",
        "",
        "| repeat | iteration | shell (Å) | cycle correction (kcal/mol) |",
        "|---:|---:|---:|---:|",
    ]
    for cycle in report["cycles"]:
        for result in cycle["shell_corrections"]:
            lines.append(f"| {cycle['repeat']} | {cycle['iteration']} | {result['shell_angstrom']:g} | {result['correction_kcal_mol']:.4f} |")
    lines.extend([
        "",
        "## 队列统计",
        "",
        "| shell (Å) | n | mean | SD | median | range |",
        "|---:|---:|---:|---:|---:|---:|",
    ])
    for shell, values in report["cycle_summary_by_shell_kcal_mol"].items():
        sd = "NA" if values["sample_sd_kcal_mol"] is None else f"{values['sample_sd_kcal_mol']:.4f}"
        lines.append(f"| {float(shell):g} | {values['n']} | {values['mean_kcal_mol']:.4f} | {sd} | {values['median_kcal_mol']:.4f} | {values['min_kcal_mol']:.4f} to {values['max_kcal_mol']:.4f} |")
    lines.extend([
        "",
        f"配对比较以 {report['reference_shell_angstrom']:g} Å 为参照：",
        "",
        "| shell (Å) | n | mean difference | SD | range |",
        "|---:|---:|---:|---:|---:|",
    ])
    for shell, values in report["paired_difference_vs_reference_shell_kcal_mol"].items():
        sd = "NA" if values["sample_sd_kcal_mol"] is None else f"{values['sample_sd_kcal_mol']:.4f}"
        lines.append(f"| {float(shell):g} | {values['n']} | {values['mean_difference_kcal_mol']:.4f} | {sd} | {values['min_kcal_mol']:.4f} to {values['max_kcal_mol']:.4f} |")
    lines.extend([
        "",
        "## 四态输入规模",
        "",
        "| repeat | iteration | state | shell (Å) | QM atoms | MM atoms | waters | ions |",
        "|---:|---:|---|---:|---:|---:|---:|---:|",
    ])
    for frame in report["frame_states"]:
        for result in frame["shell_results"]:
            lines.append(f"| {frame['repeat']} | {frame['iteration']} | {frame['leg']}_{frame['endpoint']} | {result['shell_angstrom']:g} | {result['qm_atoms']} | {result['mm_atoms']} | {result['water_residues']} | {result['ion_residues']} |")
    lines.extend([
        "",
        "## 判读边界",
        "",
        "若双差随半径明显变化，说明该局部 QM/MM 结果对水/离子截断敏感，不能把某一个半径下的值当作稳定结合修正。若变化较小，也只说明这些冻结帧上的显式壳层半径影响有限；它不证明周期 PME 长程静电或系综平均已经处理正确。本轮只覆盖一条配体变换边，不能用于亲和力排序。",
        "",
    ])
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="append", required=True, metavar="RADIUS=STATUS_JSON",
                        help="repeat for each matched shell radius, e.g. 8=run8/status.json")
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path, required=True)
    args = parser.parse_args()
    named_paths = []
    for value in args.run:
        radius, separator, raw_path = value.partition("=")
        if not separator:
            parser.error("each --run must be RADIUS=STATUS_JSON")
        try:
            shell = float(radius)
        except ValueError:
            parser.error(f"invalid shell radius: {radius}")
        if not np.isfinite(shell) or shell <= 0:
            parser.error("shell radii must be finite and positive")
        named_paths.append((shell, Path(raw_path)))
    named_paths.sort(key=lambda row: row[0])
    report = summarize(named_paths)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(report, indent=2) + "\n")
    args.output_md.parent.mkdir(parents=True, exist_ok=True)
    args.output_md.write_text(render_markdown(report))
    print(json.dumps({"status": report["status"], "fixed_four_state_cycles": report["fixed_four_state_cycles"],
                      "shell_radii_angstrom": report["shell_radii_angstrom"],
                      "cycle_summary_by_shell_kcal_mol": report["cycle_summary_by_shell_kcal_mol"],
                      "paired_difference_vs_reference_shell_kcal_mol": report["paired_difference_vs_reference_shell_kcal_mol"]}, indent=2))


if __name__ == "__main__":
    main()
