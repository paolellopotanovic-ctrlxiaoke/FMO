"""Full-cohort conditional H reconstruction control; no affinity-driven selection."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np
from rdkit import Chem

from fmotpiefinal.core.bound_geometry import load_bound_ligand
from fmotpiefinal.core.ligand_hydrogen_geometry import audit_hch_angles, rebuild_flagged_carbon_hydrogens
from run_cdk8_solution_cycle import digest, write_json
from run_local_sidechain_transfer import directory_key, evaluate, score_job


TARGETS = ("CHK1", "BACE", "Galectin3", "CDK8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source, output = args.source.resolve(), args.output.resolve()
    parent_plan = json.loads((source / "frozen_protocol.json").read_text())
    source_rows = json.loads((source / "progress.json").read_text())
    cohorts = {target: parent_plan["cohorts"][target] for target in TARGETS}
    decisions = {}
    hashes = {}
    for target, cohort in cohorts.items():
        molecules = {molecule.GetProp("_Name").strip(): molecule for molecule in Chem.SDMolSupplier(cohort["sdf"], removeHs=False) if molecule is not None}
        if set(molecules) != set(cohort["names"]):
            raise ValueError("full molecular inventory required")
        decisions[target] = {}
        for name in cohort["names"]:
            row = source_rows[target][name]
            directory = source / target / directory_key(name)
            path = directory / "minimization/minimized_complex.pdb"
            if row["status"] != "SUCCEEDED" or digest(directory / "input.inp") != row["input_sha256"] or digest(directory / "run/gamess.log") != row["log_sha256"]:
                raise ValueError("verified complete source QM required")
            positions = np.asarray([atom.xyz for atom in load_bound_ligand(cohort["sdf"], name, path).atoms])
            repaired, changed = rebuild_flagged_carbon_hydrogens(molecules[name], positions)
            decisions[target][name] = {"audit_before": audit_hch_angles(molecules[name], positions),
                "changed_hydrogen_indices": changed, "new_coordinates_A": repaired.tolist() if changed else None,
                "audit_after": audit_hch_angles(molecules[name], repaired)}
            hashes[str(path)] = digest(path)
            hashes[str(directory / "input.inp")] = digest(directory / "input.inp")
            hashes[str(directory / "run/gamess.log")] = digest(directory / "run/gamess.log")
    files = [Path(__file__), Path(__file__).parents[1] / "src/fmotpiefinal/core/ligand_hydrogen_geometry.py"]
    frozen = {"schema": "flagged-carbon-hydrogen-reconstruction-v1", "cohorts": cohorts,
        "decisions": decisions, "source_sha256": hashes, "code_sha256": {str(path.resolve()): digest(path) for path in files},
        "labels_used_for_selection": False, "scores_used_for_selection": False,
        "control": "only gross SP3 H-C-H angles outside 80..150 degrees trigger RDKit topology-based H coordinate reconstruction; no deletion, no heavy atom movement, no energy-based selection",
        "limitations": "direct reconstructed-H single point; not a complete H optimization or production repair",
        "production_default_changed": False}
    path = output / "frozen_protocol.json"
    if path.exists() and json.loads(path.read_text()) != frozen:
        raise ValueError("refusing changed frozen H reconstruction")
    write_json(path, frozen)
    rows = {target: {} for target in TARGETS}
    started = time.monotonic()
    for target, cohort in cohorts.items():
        for name in cohort["names"]:
            decision = decisions[target][name]
            source_directory = source / target / directory_key(name)
            destination = output / target / directory_key(name)
            if not decision["changed_hydrogen_indices"]:
                result = dict(source_rows[target][name], execution_path="UNCHANGED_VERIFIED_SIDECHAIN_SOURCE", new_qm_executed=False)
            else:
                lines = (source_directory / "minimization/minimized_complex.pdb").read_text().splitlines()
                ligand_index = 0
                replaced = 0
                for index, line in enumerate(lines):
                    if line.startswith("HETATM") and line[17:20] == "UNK":
                        if ligand_index in decision["changed_hydrogen_indices"]:
                            point = decision["new_coordinates_A"][ligand_index]
                            lines[index] = line[:30]+"".join(f"{value:8.3f}" for value in point)+line[54:]
                            replaced += 1
                        ligand_index += 1
                if replaced != len(decision["changed_hydrogen_indices"]):
                    raise ValueError("hydrogen identity mapping incomplete")
                complex_path = destination / "minimization/minimized_complex.pdb"
                complex_path.parent.mkdir(parents=True, exist_ok=True)
                complex_path.write_text("\n".join(lines)+"\n")
                write_json(destination / "preparation_status.json", {"status": "SUCCEEDED", "complex": str(complex_path),
                    "control": "topology-based flagged-carbon H reconstruction; all other coordinates unchanged"})
                _, _, result = score_job(cohort, name, output, 2)
                result.update(execution_path="FLAGGED_H_RECONSTRUCTION_AND_NEW_QM", new_qm_executed=True)
            rows[target][name] = result
            write_json(destination / "summary.json", result)
            write_json(output / "progress.json", rows)
            print(target, name, result["status"], result["execution_path"], result.get("tpie"), flush=True)
    evaluate(args.workspace.resolve(), args.workspace.resolve() / "FMOtpiefinal_HIV1lab/runs/all_accelerated_robust_20261020", rows, output)
    failures = sum(row["status"] != "SUCCEEDED" for group in rows.values() for row in group.values())
    write_json(output / "final_status.json", {"status": "COMPLETE" if not failures else "INCOMPLETE_WITH_FAILURES",
        "n_requested": sum(len(group) for group in rows.values()), "n_failed": failures,
        "n_new_qm": sum(row["new_qm_executed"] for group in rows.values() for row in group.values()),
        "wall_seconds": time.monotonic()-started, "evaluation_sha256": digest(output / "evaluation.json")})


if __name__ == "__main__":
    main()
