"""Report paired repeats from complete, frozen same-frame FMO controls."""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np

from run_matched_openfe import digest, write_json


COMPONENTS = {"ES": "electrostatic_kcal_mol", "EX": "exchange_repulsion_kcal_mol",
              "CT": "charge_transfer_kcal_mol", "DI": "dispersion_kcal_mol",
              "SOLV": "solvation_kcal_mol"}


def summarize(plan, rows, endpoint_names, baselines):
    expected = {task["directory"] for task in plan["tasks"]}
    actual = [row["directory"] for row in rows]
    if len(actual) != len(set(actual)) or not set(actual) <= expected:
        raise ValueError("duplicate or unknown result")
    tasks = {task["directory"]: task for task in plan["tasks"]}
    for row in rows:
        if any(row[key] != tasks[row["directory"]][key]
               for key in ("target", "name", "trajectory", "iteration", "ligand_fragment")):
            raise ValueError("result identity differs from frozen task")
    failed = [{"directory": row["directory"], "status": row["result"]["status"],
               "error": row["result"].get("error")} for row in rows if row["result"]["status"] != "SUCCEEDED"]
    report = {"status": "COMPLETE" if set(actual) == expected and not failed else "INCOMPLETE",
        "n_expected": len(expected), "n_received": len(rows), "failed": failed,
        "missing": sorted(expected - set(actual)), "targets": {}, "labels_used": False,
        "cohort_ranking_available": False, "binding_free_energy_claim_allowed": False,
        "scope": plan["scope"], "main_pipeline_changed": False}
    if report["status"] != "COMPLETE":
        return report
    iterations = plan.get("fixed_iterations", [400, 1200, 2000])
    if not iterations or sorted(set(iterations)) != list(iterations):
        raise ValueError("unique increasing fixed frame schedule required")
    for target, names in endpoint_names.items():
        first, second = names
        selected = [row for row in rows if row["target"] == target]
        if len(selected) != 6 * len(iterations) or Counter(row["name"] for row in selected) != {first: 3 * len(iterations), second: 3 * len(iterations)}:
            raise ValueError("complete two-endpoint three-repeat schedule required")
        repeats = {}
        residue_differences = []
        for trajectory in sorted({row["trajectory"] for row in selected}):
            values = {}
            for name in names:
                frames = [row for row in selected if row["trajectory"] == trajectory and row["name"] == name]
                if sorted(row["iteration"] for row in frames) != list(iterations):
                    raise ValueError("missing fixed frame schedule")
                scores = []
                for row in frames:
                    result = row["result"]
                    pairs = [pair for pair in result["pairs"] if row["ligand_fragment"] in (pair["fragment_i"], pair["fragment_j"])]
                    if len(pairs) != row["ligand_fragment"] - 1:
                        raise ValueError("incomplete ligand pair inventory")
                    other_fragments = [pair["fragment_j"] if pair["fragment_i"] == row["ligand_fragment"]
                                       else pair["fragment_i"] for pair in pairs]
                    if sorted(other_fragments) != list(range(1, row["ligand_fragment"])):
                        raise ValueError("duplicate or unknown receptor pair")
                    score = {"tPIE": result["tpie_kcal_mol"]}
                    score.update({name: sum(pair[field] for pair in pairs) for name, field in COMPONENTS.items()})
                    if not np.isfinite(list(score.values())).all():
                        raise ValueError("nonfinite score")
                    scores.append(score)
                values[name] = {key: float(np.mean([score[key] for score in scores])) for key in scores[0]}
            residue_means = {}
            for name in names:
                frames = [row for row in selected if row["trajectory"] == trajectory and row["name"] == name]
                frame_pairs = []
                for row in frames:
                    fragment = row["ligand_fragment"]
                    frame_pairs.append({pair["fragment_j"] if pair["fragment_i"] == fragment else pair["fragment_i"]:
                        pair["value_kcal_mol"] for pair in row["result"]["pairs"]
                        if fragment in (pair["fragment_i"], pair["fragment_j"])})
                if any(set(pairs) != set(frame_pairs[0]) for pairs in frame_pairs):
                    raise ValueError("receptor pair identities changed across frames")
                residue_means[name] = {fragment: float(np.mean([pairs[fragment] for pairs in frame_pairs]))
                                       for fragment in frame_pairs[0]}
            if set(residue_means[first]) != set(residue_means[second]):
                raise ValueError("receptor fragments differ across ligand endpoints")
            residue_differences.append({fragment: residue_means[second][fragment] - residue_means[first][fragment]
                                        for fragment in residue_means[first]})
            repeats[trajectory] = {"ligand_means_kcal_mol": values,
                "B_minus_A_kcal_mol": {key: values[second][key] - values[first][key] for key in values[first]}}
        if len(repeats) != 3:
            raise ValueError("three independent repeats required")
        differences = [repeat["B_minus_A_kcal_mol"] for repeat in repeats.values()]
        averaged = {key: float(np.mean([value[key] for value in differences])) for key in differences[0]}
        deviations = {key: float(np.std([value[key] for value in differences], ddof=1)) for key in differences[0]}
        report["targets"][target] = {"first": first, "second": second, "n_independent_repeats": 3,
            "repeats": repeats, "mean_B_minus_A_kcal_mol": averaged,
            "repeat_std_B_minus_A_kcal_mol": deviations,
            "baseline_tpie_B_minus_A_kcal_mol": baselines[target][second] - baselines[target][first],
            "receptor_fragment_B_minus_A_kcal_mol": [{"fragment_id": fragment,
                "mean": float(np.mean([values[fragment] for values in residue_differences])),
                "repeat_std": float(np.std([values[fragment] for values in residue_differences], ddof=1))}
                for fragment in sorted(residue_differences[0])],
            "ligands": {name: {"baseline_tpie_kcal_mol": baselines[target][name],
                "mean_tpie_kcal_mol": float(np.mean([repeat["ligand_means_kcal_mol"][name]["tPIE"] for repeat in repeats.values()])),
                "repeat_std_tpie_kcal_mol": float(np.std([repeat["ligand_means_kcal_mol"][name]["tPIE"] for repeat in repeats.values()], ddof=1))}
                for name in names}, "convergence_proven": False}
        time_slices = {}
        for iteration in iterations:
            differences_at_time = []
            for trajectory in repeats:
                frames = {row["name"]: row for row in selected
                          if row["iteration"] == iteration and row["trajectory"] == trajectory}
                differences_at_time.append(frames[second]["result"]["tpie_kcal_mol"] - frames[first]["result"]["tpie_kcal_mol"])
            time_slices[str(iteration)] = {"production_ns": iteration * .0025,
                "mean_B_minus_A_kcal_mol": float(np.mean(differences_at_time)),
                "repeat_std_B_minus_A_kcal_mol": float(np.std(differences_at_time, ddof=1)),
                "n_independent_repeats": 3}
        report["targets"][target]["fixed_time_slice_sensitivity"] = time_slices
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--openfe", type=Path, required=True)
    args = parser.parse_args()
    root = args.run.resolve()
    plan = json.loads((root / "frozen_protocol.json").read_text())
    rows = json.loads((root / "progress.json").read_text())
    controls = json.loads((args.baseline / "progress.json").read_text())
    baselines = {}
    for row in controls:
        if row["kind"] == "baseline_replay":
            if row["result"]["status"] != "SUCCEEDED":
                raise ValueError("verified baseline replay required")
            baselines.setdefault(row["target"], {})[row["name"]] = row["result"]["tpie_kcal_mol"]
    openfe = json.loads((args.openfe / "free_energy_snapshot.json").read_text())
    names = {target: [data["first"], data["second"]] for target, data in openfe["targets"].items()}
    report = summarize(plan, rows, names, baselines)
    if report["status"] == "COMPLETE":
        from run_cdk8_solution_cycle import array, group

        for target, data in report["targets"].items():
            sources = {task["name"]: task["source"] for task in plan["tasks"] if task["target"] == target}
            fragment_names = {}
            for name, source in sources.items():
                text = Path(source).read_text()
                fragment_names[name] = array(group(text, "FMO")[1], "FRGNAM", "INDAT")[:-1]
            first, second = names[target]
            if fragment_names[first] != fragment_names[second]:
                raise ValueError("baseline receptor fragment labels differ across endpoints")
            for row in data["receptor_fragment_B_minus_A_kcal_mol"]:
                row["fragment_name"] = fragment_names[first][row["fragment_id"] - 1]
            baseline_pairs = {}
            for control in controls:
                if control["target"] == target and control["kind"] == "baseline_replay":
                    ligand_fragment = control["result"]["ligand_fragment"]
                    baseline_pairs[control["name"]] = {pair["fragment_j"] if pair["fragment_i"] == ligand_fragment else pair["fragment_i"]:
                        pair["value_kcal_mol"] for pair in control["result"]["pairs"]
                        if ligand_fragment in (pair["fragment_i"], pair["fragment_j"])}
            for row in data["receptor_fragment_B_minus_A_kcal_mol"]:
                fragment = row["fragment_id"]
                row["baseline_B_minus_A"] = baseline_pairs[second][fragment] - baseline_pairs[first][fragment]
                row["ensemble_shift_from_baseline"] = row["mean"] - row["baseline_B_minus_A"]
    report["report_driver_sha256"] = digest(__file__)
    report["source_sha256"] = {str(path): digest(path) for path in (root / "frozen_protocol.json",
        root / "progress.json", args.baseline / "progress.json", args.openfe / "free_energy_snapshot.json")}
    write_json(root / "evaluation.json", report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
