"""Finite-difference check of bound ligand analytic PCM/DFTB gradients."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import io
import json
from pathlib import Path
import re
import sys

from run_cdk8_bound_qm import replace_coordinates
from run_cdk8_solution_cycle import digest, group, write_json


def first_gradient(path, atom_number):
    with Path(path).open() as handle:
        for line in handle:
            if "GRADIENT (HARTREE/BOHR)" not in line:
                continue
            for row in handle:
                fields = row.split()
                if len(fields) == 6 and fields[0].isdigit() and int(fields[0]) == atom_number:
                    return list(map(float, fields[-3:]))
            break
    raise ValueError("missing requested initial analytic gradient")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--bound-run", type=Path, required=True)
    args = parser.parse_args()
    workspace, reference, bound_run = args.workspace.resolve(), args.reference.resolve(), args.bound_run.resolve()
    sys.path.insert(0, str(workspace / "src"))
    from sophosqm_baseline.relative_total_energy import execute_deck
    from fmotpiefinal.core.fmo_convergence import final_fmo_convergence
    root = bound_run / "gradient_check"
    root.mkdir(exist_ok=True)
    parameters = (reference / "pcm/13/params").resolve(strict=True)
    tasks, definitions = [], {}
    for name in ("13", "17"):
        deck = (reference / f"pcm/{name}/complex/input.inp").read_text()
        rows = group(deck, "FMOXYZ")[1].strip().splitlines()
        ligand_size = len(group((reference / f"pcm/{name}/ligand/input.inp").read_text(), "FMOXYZ")[1].strip().splitlines())
        atom_number = len(rows) - ligand_size + 10
        if rows[atom_number-1].split()[0] != "O":
            raise ValueError("expected atom-mapped carbonyl oxygen")
        gradient = first_gradient(bound_run / name / "optimization/gamess.log", atom_number)
        axis = max(range(3), key=lambda index: abs(gradient[index]))
        definitions[name] = {"atom_number": atom_number, "axis": axis, "analytic_hartree_bohr": gradient[axis]}
        for step in (0.001, 0.0005):
            for sign in (-1, 1):
                changed = rows.copy()
                fields = changed[atom_number-1].split()
                fields[axis+2] = f"{float(fields[axis+2]) + sign * step:.10f}"
                changed[atom_number-1] = " ".join(fields)
                tasks.append((name, step, sign, replace_coordinates(deck, changed)))
    write_json(root / "plan.json", {"driver_sha256": digest(__file__), "definitions": definitions,
        "steps_angstrom": [0.001, 0.0005], "labels_used": False,
        "selection": "Primary amide 13 and anomalously twisted cyclopropyl amide 17; diagnostic, no scoring selection"})

    def execute(task):
        name, step, sign, deck = task
        parent = root / name
        parent.mkdir(exist_ok=True)
        link = parent / "params"
        if not link.exists():
            link.symlink_to(parameters, target_is_directory=True)
        directory = parent / f"h{step}_{sign:+d}"
        directory.mkdir(exist_ok=True)
        (directory / "input.inp").write_text(deck)
        runtime = json.loads((reference / f"pcm/{name}/complex/result.json").read_text())
        result = execute_deck(deck, directory, gamess_dir=runtime["gamess_dir"], library_dir=runtime["library_dir"], ncores=2, timeout=1800)
        if result["status"] == "SUCCEEDED":
            final_fmo_convergence(deck, io.StringIO((directory / "gamess.log").read_text()))
        return name, step, sign, result

    with ThreadPoolExecutor(max_workers=4) as executor:
        results = list(executor.map(execute, tasks))
    differences = []
    for name, definition in definitions.items():
        for step in (0.001, 0.0005):
            selected = {sign: value for key, size, sign, value in results if key == name and size == step}
            if any(value["status"] != "SUCCEEDED" for value in selected.values()):
                differences.append({"ligand": name, "step_angstrom": step, "status": "FAILED"})
                continue
            numerical = (selected[1]["energy_hartree"] - selected[-1]["energy_hartree"]) / (2 * step) * 0.529177210903
            differences.append({"ligand": name, "step_angstrom": step, "status": "COMPLETE",
                **definition, "numerical_hartree_bohr": numerical,
                "difference_hartree_bohr": numerical - definition["analytic_hartree_bohr"]})
    write_json(root / "evaluation.json", {"rows": differences,
        "warning": "Two atoms/axes only, not exhaustive gradient validation; PCM tessellation and energy precision affect finite differences"})
    print(json.dumps(differences, indent=2))


if __name__ == "__main__":
    main()
