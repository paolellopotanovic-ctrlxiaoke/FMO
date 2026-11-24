"""Evaluate only a complete, audited frozen-receptor QM geometry cohort."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from run_cdk8_solution_cycle import HARTREE_TO_KCAL, metrics, write_json
from fmotpiefinal.core.solution_cycle import require_consistent_pcm_cycle


def geometry_scores(complex_energy, receptor_energy, bound_energy, free_energy, tpie):
    bound = (complex_energy - receptor_energy - bound_energy) * HARTREE_TO_KCAL
    strain = (bound_energy - free_energy) * HARTREE_TO_KCAL
    return {"qm_bound_tpie": tpie, "qm_bound_three_state": bound,
        "qm_free_three_state": bound + strain,
        "qm_bound_strain_hybrid": tpie + strain, "qm_strain_kcal_mol": strain}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    args = parser.parse_args()
    reference, run = args.reference.resolve(), args.run.resolve()
    protocol = json.loads((run / "frozen_protocol.json").read_text())
    names = json.loads((reference / "frozen_protocol.json").read_text())["names"]
    results = {}
    for name in names:
        summary = run / name / "summary.json"
        if summary.exists():
            results[name] = json.loads(summary.read_text())
    missing = [name for name in names if results.get(name, {}).get("status") != "SUCCEEDED"]
    if set(protocol["names"]) != set(names):
        raise ValueError("partial pilot cohort cannot publish full-cohort metrics")
    if missing:
        write_json(run / "evaluation.json", {"status": "INCOMPLETE", "n_expected": len(names),
            "n_succeeded": len(names) - len(missing), "missing_or_failed": missing, "metrics": {}})
        print("INCOMPLETE", missing)
        return
    previous = json.loads((reference / "final_evaluation.json").read_text())
    sampling = json.loads((reference / "solution_sampling4/final_evaluation.json").read_text())
    if sampling.get("failed_starts"):
        raise ValueError("free-state sampling is incomplete")
    baseline = {row["ligand"]: row for row in previous["rows"]}
    sampled = {row["ligand"]: row for row in sampling["rows"]}
    receptor_path = reference / "pcm/shared/receptor"
    receptor_energy = json.loads((receptor_path / "result.json").read_text())["energy_hartree"]
    rows, audits = [], {}
    for name in names:
        state_paths = {"complex": run / name / "complex", "ligand": run / name / "ligand",
                       "receptor": receptor_path}
        audits[name] = require_consistent_pcm_cycle(
            {state: (path / "input.inp").read_text() for state, path in state_paths.items()},
            {state: (path / "gamess.log").read_text() for state, path in state_paths.items()})
        bound_reference = json.loads((reference / f"pcm/{name}/ligand/result.json").read_text())["energy_hartree"]
        free_energy = bound_reference - sampled[name]["solution_min_strain"] / HARTREE_TO_KCAL
        result = results[name]
        row = {"ligand": name, "experimental_dg": baseline[name]["experimental_dg"],
            "original_tpie": baseline[name]["tpie_pcm"], **geometry_scores(
                result["energies_hartree"]["complex"], receptor_energy,
                result["energies_hartree"]["ligand"], free_energy, result["tpie_kcal_mol"])}
        rows.append(row)
    arms = ("original_tpie", "qm_bound_tpie", "qm_bound_three_state", "qm_bound_strain_hybrid")
    truth = [row["experimental_dg"] for row in rows]
    evaluated = {arm: metrics(truth, [row[arm] for row in rows]) for arm in arms}
    write_json(run / "evaluation.json", {"status": "COMPLETE", "n": len(rows), "rows": rows,
        "metrics": evaluated, "pcm_audits": audits, "production_changed": False,
        "warning": "Development-set geometry control, not a free energy or external validation"})
    print(json.dumps(evaluated, indent=2))


if __name__ == "__main__":
    main()
