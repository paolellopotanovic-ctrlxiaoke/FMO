"""Gate embedded QM/MM forces against same-model Cartesian energy differences."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

import numpy as np

from audit_quanpol_mm_parity import evaluate_export, read_sections
from diagnose_quanpol_internal_modes import rebind_parameters
from run_matched_openfe import digest, write_json

HARTREE_KCAL = 627.509469
BOHR_ANGSTROM = .52917724924


def parse_embedded(text):
    if "EXECUTION OF GAMESS TERMINATED NORMALLY" not in text or "SCF IS UNCONVERGED" in text:
        raise ValueError("normal SCF-converged termination required")
    blocks = text.split("QMMM_FORCE_EXPORT_BEGIN")
    if len(blocks) != 2 or "QMMM_FORCE_EXPORT_END" not in blocks[1]:
        raise ValueError("exactly one unmasked embedded-gradient block required")
    block = blocks[1].split("QMMM_FORCE_EXPORT_END")[0]
    header = block.splitlines()[0].split()
    nat, nff = map(int, header[:2])
    energies = re.findall(r"POTENTIAL\s+ENERGY =\s+([-+0-9.EeDd]+)\s+KCAL/MOL", text)
    if len(energies) != 1:
        raise ValueError("one zero-step QuanPol potential energy required")
    energy = float(energies[0].replace("D", "E"))
    qm, mm, mapping = [], [], []
    for line in block.splitlines()[1:]:
        words = line.split()
        if not words:
            continue
        if words[0] == "QMG":
            if int(words[1]) != len(qm) + 1:
                raise ValueError("QM gradient order mismatch")
            qm.append([float(value.replace("D", "E")) for value in words[2:]])
        elif words[0] == "MMG":
            if int(words[1]) != len(mm) + 1:
                raise ValueError("MM gradient order mismatch")
            mapping.append(int(words[2]))
            mm.append([float(value.replace("D", "E")) for value in words[3:]])
    qm, mm, mapping = np.asarray(qm), np.asarray(mm), np.asarray(mapping)
    if qm.shape != (nat, 3) or mm.shape != (nff, 3) or not np.isfinite(mm).all() or not np.isfinite(qm).all() or not np.isfinite(energy):
        raise ValueError("incomplete or nonfinite embedded gradients")
    if not np.array_equal(mapping[:nat], np.arange(1, nat + 1)) or np.any(mapping[nat:] != 0):
        raise ValueError("requires one-to-one first-NAT ligand mapping, no link atoms")
    if not np.allclose(qm, mm[:nat], atol=1e-12, rtol=0):
        raise ValueError("duplicated QM/MM atom gradients disagree")
    return {"energy_kcal_mol": energy, "gradient_kcal_mol_angstrom": mm * HARTREE_KCAL / BOHR_ANGSTROM,
            "n_qm": nat, "n_physical": nff}


def replace_coordinates(text, coordinates, n_qm):
    coordinates = np.asarray(coordinates, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3 or not np.isfinite(coordinates).all():
        raise ValueError("finite Nx3 coordinates required")
    lines = text.splitlines()
    data = next(index for index, line in enumerate(lines) if line.strip().upper() == "$DATA") + 3
    data_end = next(index for index in range(data, len(lines)) if lines[index].strip().upper() == "$END")
    ff = next(index for index, line in enumerate(lines) if line.strip().upper() == "$FFDATA")
    ff_start = next(index for index in range(ff + 1, len(lines)) if lines[index].strip().upper() == "COORDINATES") + 1
    ff_end = next(index for index in range(ff_start, len(lines)) if lines[index].strip().upper() == "STOP")
    if data_end - data != n_qm or ff_end - ff_start != len(coordinates):
        raise ValueError("coordinate inventories do not match deck")
    for start, count in ((data, n_qm), (ff_start, len(coordinates))):
        for offset in range(count):
            words = lines[start + offset].split()
            if len(words) != 5:
                raise ValueError("unexpected coordinate row")
            lines[start + offset] = " " + " ".join(words[:2]) + " " + " ".join(f"{value:.12f}" for value in coordinates[offset])
    return "\n".join(lines) + "\n"


def select_directions(sections, n_qm, gradient):
    xyz = np.array([row[2:5] for row in sections["COORDINATES"]], dtype=float)
    masses = np.array([float(row[1]) for row in sections["PARAMETERS"]])
    distance = np.min(np.linalg.norm(xyz[:, None] - xyz[None, :n_qm], axis=-1), axis=1)
    neighbors = {index: [] for index in range(len(xyz))}
    for row in sections.get("BOND", []):
        first, second = int(row[1]) - 1, int(row[2]) - 1
        neighbors[first].append(second)
        neighbors[second].append(first)
    water = [index for index in range(n_qm, len(xyz)) if 15 < masses[index] < 17 and
             len(neighbors[index]) == 2 and all(masses[other] < 2 for other in neighbors[index])]
    other_heavy = [index for index in range(n_qm, len(xyz)) if masses[index] > 2 and index not in water]
    choices = [("qm_heavy", next(index for index in range(n_qm) if masses[index] > 2)),
               ("qm_hydrogen", next(index for index in range(n_qm) if masses[index] < 2))]
    if water:
        oxygen = min(water, key=lambda index: distance[index])
        choices.append(("mm_water_oxygen", oxygen))
    if other_heavy:
        choices.append(("mm_nonwater_heavy", min(other_heavy, key=lambda index: distance[index])))
    elif water:
        choices.append(("mm_water_hydrogen", neighbors[oxygen][0]))
    return [{"kind": kind, "atom": atom, "axis": int(np.argmax(np.abs(gradient[atom])))} for kind, atom in choices]


def direction_vector(direction, shape, n_qm):
    vector = np.zeros(shape)
    if "atom" in direction:
        vector[direction["atom"], direction["axis"]] = 1
    elif direction["kind"] == "collective_qm_translation":
        vector[:n_qm, 0] = 1
    else:
        vector = np.random.default_rng(direction["seed"]).normal(size=shape)
        if direction["kind"] == "collective_environment":
            vector[:n_qm] = 0
    return vector / np.linalg.norm(vector)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent-status", type=Path, required=True)
    parser.add_argument("--force-build", type=Path, required=True)
    parser.add_argument("--parameters", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    parent = json.loads(args.parent_status.read_text())
    if parent.get("status") != "COMPLETE" or parent.get("labels_used") is not False:
        raise ValueError("complete label-free parent required")
    selected = [row for row in parent["results"] if row["repeat"] == 0 and row["iteration"] == 400]
    if {(row["leg"], row["endpoint"]) for row in selected} != {(leg, endpoint) for leg in ("bound", "solution") for endpoint in (0, 10)}:
        raise ValueError("four fixed representative states required")
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    alias = Path("/tmp/qpforce_sk_20261101")
    if not alias.exists():
        alias.symlink_to(args.parameters.resolve(), target_is_directory=True)
    if alias.resolve() != args.parameters.resolve():
        raise ValueError("parameter alias points to different source")
    environment = os.environ.copy()
    environment.update(LD_LIBRARY_PATH=str(Path(sys.executable).resolve().parents[1] / "lib") + ":" + environment.get("LD_LIBRARY_PATH", ""),
                       OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    report = {"status": "RUNNING", "labels_used": False, "production_pipeline_modified": False,
              "scope": "full embedded ligand-QM/environment-MM finite-shell Cartesian forces, four states; not full periodic QM/MM",
              "build_manifest": str(args.force_build.resolve() / "build_manifest.json"), "states": []}
    write_json(args.output / "status.json", report)

    def execute(deck):
        log = deck.with_suffix(".stdout")
        if log.exists():
            try:
                return parse_embedded(log.read_text())
            except ValueError:
                pass
        started = time.monotonic()
        with log.open("w") as handle:
            process = subprocess.run(["timeout", "180", "tcsh", "-f", str(args.force_build.resolve() / "rungms"), str(deck), "force", "2"],
                                     cwd=args.output, env=environment, stdout=handle, stderr=subprocess.STDOUT)
        if process.returncode:
            raise RuntimeError(f"GAMESS launcher failed: {deck}: {process.returncode}")
        result = parse_embedded(log.read_text())
        write_json(deck.with_suffix(".timing.json"), {"wall_seconds": time.monotonic() - started, "deck_sha256": digest(deck), "log_sha256": digest(log)})
        return result

    try:
        for row in selected:
            key = f"{row['leg'][0]}{row['endpoint']}"
            original = Path(row["qm"]["log"]).parent
            source = next(original.glob("qm_*.inp"))
            text, _ = rebind_parameters(source.read_text(), alias)
            text = re.sub(r"RUNTYP=\w+", "RUNTYP=GRADIENT", text, count=1)
            deck = args.output / f"efv_{key}_base.inp"
            deck.write_text(text)
            base = execute(deck)
            sections = read_sections(text)
            xyz = np.array([record[2:5] for record in sections["COORDINATES"]], dtype=float)
            gradients = base["gradient_kcal_mol_angstrom"]
            reference_mm_deck = next(original.glob("mm_*.inp"))
            reference_mm = evaluate_export(reference_mm_deck, include_forces=True)
            delta_gradient = gradients + np.array(reference_mm["forces_kcal_mol_angstrom"])
            directions = select_directions(sections, base["n_qm"], gradients)
            directions += [{"kind": "collective_all", "seed": 47103},
                           {"kind": "collective_environment", "seed": 47104},
                           {"kind": "collective_qm_translation"}]
            jobs = []
            for ordinal, direction in enumerate(directions):
                vector = direction_vector(direction, xyz.shape, base["n_qm"])
                for step in (.001, .0002):
                    for sign in (-1, 1):
                        perturbed = xyz + sign * step * vector
                        perturbed_deck = args.output / f"efv_{key}_d{ordinal}_h{step:g}_{'p' if sign > 0 else 'm'}.inp"
                        perturbed_deck.write_text(replace_coordinates(text, perturbed, base["n_qm"]))
                        jobs.append((ordinal, step, sign, perturbed_deck))
            with ThreadPoolExecutor(max_workers=args.workers) as pool:
                results = list(pool.map(lambda job: execute(job[3]), jobs))
            energies = {(ordinal, step, sign): result["energy_kcal_mol"] for (ordinal, step, sign, _), result in zip(jobs, results)}
            comparisons = []
            for ordinal, direction in enumerate(directions):
                analytic = float(np.sum(gradients * direction_vector(direction, xyz.shape, base["n_qm"])))
                for step in (.001, .0002):
                    numerical = (energies[ordinal, step, 1] - energies[ordinal, step, -1]) / (2 * step)
                    comparisons.append({**direction, "step_angstrom": step, "analytic_kcal_mol_angstrom": analytic,
                                        "numerical_kcal_mol_angstrom": numerical, "absolute_error": abs(analytic - numerical),
                                        "passed": abs(analytic - numerical) < .005 + abs(analytic) * 2e-5})
            net = float(np.linalg.norm(gradients.sum(axis=0)))
            parity = abs(base["energy_kcal_mol"] - row["qm"]["potential_kcal_mol"])
            np.savez_compressed(args.output / f"{key}_gradients.npz", qm_mm_gradient=gradients, delta_gradient=delta_gradient,
                                coordinates_angstrom=xyz, mm_forces=np.array(reference_mm["forces_kcal_mol_angstrom"]))
            record = {"state": key, "leg": row["leg"], "endpoint": row["endpoint"], "n_qm": base["n_qm"], "n_physical": base["n_physical"],
                      "instrumentation_energy_difference_kcal_mol": parity, "net_gradient_kcal_mol_angstrom": net,
                      "delta_gradient_norm_kcal_mol_angstrom": float(np.linalg.norm(delta_gradient)),
                      "finite_differences": comparisons,
                      "passed": parity < 1e-6 and net < .01 and all(item["passed"] for item in comparisons)}
            report["states"].append(record)
            write_json(args.output / "status.json", report)
            print(json.dumps({"state": key, "passed": record["passed"], "net_gradient": net, "max_error": max(item["absolute_error"] for item in comparisons)}), flush=True)
        report["status"] = "PASSED" if all(row["passed"] for row in report["states"]) else "FAILED_VALIDATION"
        report["bridge_sampling_permitted"] = report["status"] == "PASSED"
    except Exception as error:
        report.update(status="ERROR", bridge_sampling_permitted=False, error=repr(error))
        write_json(args.output / "status.json", report)
        raise
    write_json(args.output / "status.json", report)


if __name__ == "__main__":
    main()
