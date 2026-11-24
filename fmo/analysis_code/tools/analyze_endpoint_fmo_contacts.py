"""Join matched structural contact graphs and FMO pairs without score fitting."""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from run_cdk8_solution_cycle import array, group
from run_matched_openfe import digest, write_json


def fragment_label(contact_key):
    residue, identity = contact_key.split(":", 1)
    fields = identity.split()
    if len(fields) != 2:
        raise ValueError("ambiguous contact residue identity")
    return fields[0] + residue + fields[1]


def descriptive_correlation(first, second):
    if len(first) < 3 or np.ptp(first) < 1e-10 or np.ptp(second) < 1e-10:
        return None
    return float(spearmanr(first, second).statistic)


def rotamer_fractions(angles):
    centers = np.array([-60., 60., 180.])
    delta = (np.asarray(angles)[:, None] - centers[None] + 180) % 360 - 180
    assignment = np.abs(delta).argmin(axis=1)
    return {name: float(np.mean(assignment == index))
            for index, name in enumerate(("gauche_minus", "gauche_plus", "trans"))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--structure", type=Path, required=True)
    args = parser.parse_args()
    root = args.run.resolve()
    plan_path, progress_path = root / "frozen_protocol.json", root / "progress.json"
    plan = json.loads(plan_path.read_text())
    rows = json.loads(progress_path.read_text())
    report = {"status": "INCOMPLETE", "targets": {}, "labels_used": False, "main_pipeline_changed": False,
        "source_sha256": {str(path): digest(path) for path in (plan_path, progress_path, args.structure)},
        "driver_sha256": digest(__file__), "limitations": [
            "Within-repeat correlations are descriptive; adjacent frames are correlated and no significance claim is made.",
            "Whole-residue contacts and FACIO fragment atom memberships are not identical.",
            "Distance contacts and three-bin rotamer summaries do not prove a hydrogen bond or causal affinity mechanism.",
            "No regression coefficient, score weight or frame selection is learned from experimental affinity."]}
    expected = {task["directory"] for task in plan["tasks"]}
    if len(rows) != len(expected) or {row["directory"] for row in rows} != expected or any(row["result"]["status"] != "SUCCEEDED" for row in rows):
        write_json(root / "contact_energy_analysis.json", report)
        return
    structure = json.loads(args.structure.read_text())
    for target, data in structure["targets"].items():
        structural = {}
        for repeat in data["bound"]["complete_repeats"]:
            for frame in repeat["frames"]:
                structural[repeat["trajectory"], repeat["endpoint"], frame["iteration"]] = frame
        selected = [row for row in rows if row["target"] == target]
        grouped = defaultdict(list)
        for row in selected:
            frame = structural[row["trajectory"], row["endpoint"], row["iteration"]]
            if frame["replica"] != row["replica"]:
                raise ValueError("structural endpoint replica differs from FMO")
            grouped[row["name"], row["trajectory"]].append((row, frame))
        records = []
        for (name, trajectory), frames in sorted(grouped.items()):
            frames.sort(key=lambda pair: pair[0]["iteration"])
            if [row["iteration"] for row, frame in frames] != plan["fixed_iterations"]:
                raise ValueError("complete fixed graph schedule required")
            source = Path(frames[0][0]["source"])
            labels = array(group(source.read_text(), "FMO")[1], "FRGNAM", "INDAT")[:-1]
            energy_by_fragment = {}
            for label_index, label in enumerate(labels, 1):
                values = []
                for row, frame in frames:
                    ligand = row["ligand_fragment"]
                    pair = [pair for pair in row["result"]["pairs"]
                            if {pair["fragment_i"], pair["fragment_j"]} == {ligand, label_index}]
                    if len(pair) != 1:
                        raise ValueError("unique ligand pair required")
                    values.append(pair[0]["value_kcal_mol"])
                energy_by_fragment[label] = values
            contacts = {}
            for key in frames[0][1]["residue_contacts"]:
                label = fragment_label(key)
                if label not in energy_by_fragment:
                    continue
                distance = [frame["residue_contacts"][key] for row, frame in frames]
                energy = energy_by_fragment[label]
                record = {"mean_contact_distance_angstrom": float(np.mean(distance)),
                    "fraction_below_3p5_angstrom": float(np.mean(np.asarray(distance) < 3.5)),
                    "mean_fragment_pie_kcal_mol": float(np.mean(energy)),
                    "within_repeat_pie_std_kcal_mol": float(np.std(energy, ddof=1)),
                    "within_repeat_distance_pie_spearman_descriptive": descriptive_correlation(distance, energy)}
                if key in frames[0][1]["sidechain_chi1_degrees"]:
                    record["chi1_rotamer_fractions"] = rotamer_fractions([frame["sidechain_chi1_degrees"][key] for row, frame in frames])
                contacts[label] = record
            tpie = [row["result"]["tpie_kcal_mol"] for row, frame in frames]
            descriptors = {}
            for field in ("local_receptor_fit_ligand_pose_rmsd_angstrom", "ligand_self_fit_rmsd_angstrom",
                          "water_oxygen_within_3p5A_of_ligand_heteroatom", "distance_only_water_bridge_count"):
                values = [frame[field] for row, frame in frames]
                descriptors[field] = {"mean": float(np.mean(values)),
                    "within_repeat_tpie_spearman_descriptive": descriptive_correlation(values, tpie)}
            records.append({"name": name, "trajectory": trajectory, "n_frames": len(frames),
                            "residue_contacts": contacts, "descriptors": descriptors})
        report["targets"][target] = records
    report["status"] = "COMPLETE"
    write_json(root / "contact_energy_analysis.json", report)


if __name__ == "__main__":
    main()
