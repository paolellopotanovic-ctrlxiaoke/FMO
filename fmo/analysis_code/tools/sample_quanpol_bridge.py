"""Sample explicit local-QM/MM bridge states using Metropolized periodic-MM proposals.

The target is U_MM,PME + lambda * (U_QMMM,shell - U_MM,shell).
This is not periodic QM electrostatics and does not modify production tPIE.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
import resource
import shutil
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

from run_matched_openfe import digest, write_json
from validate_quanpol_embedded_forces import parse_embedded, replace_coordinates

KJ_PER_KCAL = 4.184
GAS_KJ = .00831446261815324


def log_acceptance(mm_energy_change_kj, kinetic_change_kj, delta_change_kcal, coupling, temperature):
    values = (mm_energy_change_kj, kinetic_change_kj, delta_change_kcal, coupling, temperature)
    if not all(math.isfinite(value) for value in values) or not 0 <= coupling <= 1 or temperature <= 0:
        raise ValueError("finite energies, valid coupling, and positive temperature required")
    return -(mm_energy_change_kj + kinetic_change_kj + coupling * delta_change_kcal * KJ_PER_KCAL) / (GAS_KJ * temperature)


def compose_gradient(full_mm_gradient, shell_delta_gradient, indices):
    full = np.asarray(full_mm_gradient, dtype=float).copy()
    local = np.asarray(shell_delta_gradient, dtype=float)
    indices = np.asarray(indices, dtype=int)
    if len(set(indices.tolist())) != len(indices) or local.shape != (len(indices), 3) or np.any(indices < 0) or np.any(indices >= len(full)):
        raise ValueError("invalid local/full force mapping")
    full[indices] += local
    return full


def validate_composed_target(context, shell_context, source_positions, shell_reference, indices, gradient_data, record, validation_root):
    from openmm import unit
    from validate_quanpol_embedded_forces import direction_vector
    from diagnose_quanpol_periodic_force_gate import cutoff_jump_diagnostics

    context.setPositions(source_positions * unit.nanometer)
    state = context.getState(getForces=True)
    mm_gradient = -np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)) / 41.84
    target_gradient = compose_gradient(mm_gradient, gradient_data["delta_gradient"], indices)
    comparisons = []
    cutoff_cache = {}
    box = np.asarray(context.getState().getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer))
    for ordinal, direction in enumerate(record["finite_differences"][::2]):
        vector = direction_vector(direction, shell_reference.shape, record["n_qm"])
        analytic = float(np.sum(target_gradient[indices] * vector))
        cutoff = cutoff_jump_diagnostics(context.getSystem(), context, source_positions, box, indices, vector, [.001, .0002], cache=cutoff_cache)
        for step, jump in zip((.001, .0002), cutoff):
            energies = []
            for sign in (-1, 1):
                coordinates = shell_reference + sign * step * vector
                perturbed = source_positions.copy()
                perturbed[indices] += sign * step * .1 * vector
                context.setPositions(perturbed * unit.nanometer)
                periodic = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
                shell_context.setPositions(coordinates * .1 * unit.nanometer)
                finite_mm = shell_context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
                path = validation_root / f"efv_{record['state']}_d{ordinal}_h{step:g}_{'p' if sign > 0 else 'm'}.stdout"
                qm_mm = parse_embedded(path.read_text())["energy_kcal_mol"]
                energies.append(periodic + qm_mm - finite_mm)
            numerical = (energies[1] - energies[0]) / (2 * step)
            accounted = numerical - jump["predicted_gradient_discontinuity"]
            comparisons.append({"kind": direction["kind"], "step_angstrom": step, "analytic": analytic,
                                "numerical": numerical, "raw_absolute_error": abs(numerical - analytic),
                                "reference_cutoff_jump": jump, "numerical_with_reference_cutoff_jump_removed": accounted,
                                "absolute_error": abs(accounted - analytic),
                                "passed": abs(accounted - analytic) < .02 + 2e-5 * abs(analytic)})
    context.setPositions(source_positions * unit.nanometer)
    return {"passed": all(row["passed"] for row in comparisons), "checks": comparisons,
            "reference_pme_net_gradient_kcal_mol_angstrom": float(np.linalg.norm(mm_gradient.sum(axis=0))),
            "reference_hamiltonian_modified": False, "sampling_energy_or_force_shifted_to_pass_gate": False,
            "cutoff_discontinuities": "analytic known pair potential jumps accounted ONLY for finite-difference verification; all unmodified energies used in Metropolis acceptance",
            "net_target_gradient_kcal_mol_angstrom": float(np.linalg.norm(target_gradient.sum(axis=0)))}


def physical_indices(task, directory, saved):
    from gufe import LigandAtomMapping
    from analyze_endpoint_internal_energy import physical_atom_maps
    from analyze_openfe_endpoints import read_topology

    topology = Path(task["topology"])
    mapping = LigandAtomMapping.from_json(topology.parents[2] / "mapping.json")
    ligand, _ = physical_atom_maps(read_topology(topology), mapping)
    data = np.load(directory / "physical_coordinates.npz")
    result = np.concatenate((saved[ligand[task["endpoint"]]], data["environment_original_indices"].astype(int)))
    if len(set(result.tolist())) != len(result):
        raise ValueError("duplicate physical-to-system indices")
    return result


def read_endpoint_state(reporter, endpoint):
    import copy
    from openmmtools.utils import deserialize

    def bottom(state):
        while "thermodynamic_state" in state:
            state = state["thermodynamic_state"]
        return state

    def serialization(key):
        state = copy.deepcopy(reporter.read_dict(key))
        inner = bottom(state)
        reference = inner.pop("_Reporter__compatible_state", None)
        if reference is not None:
            group, state_id = reference.rsplit("/", 1)
            if state_id.isdigit():
                reference = f"{group}/state{state_id}"
            inner["standard_system"] = bottom(serialization(reference))["standard_system"]
        return state

    return deserialize(serialization(f"thermodynamic_states/state{endpoint}"))


def run_state(args):
    import openmm
    from openmm import unit
    from netCDF4 import Dataset
    from openmmtools.integrators import VelocityVerletIntegrator
    from openmmtools.multistate import MultiStateReporter
    from audit_quanpol_mm_parity import build_export_system

    validation = json.loads(args.validation.read_text())
    if validation.get("status") != "PASSED" or validation.get("bridge_sampling_permitted") is not True:
        raise ValueError("embedded force validation must pass before sampling")
    build_manifest_path = Path(validation["build_manifest"])
    if build_manifest_path.parent.resolve() != args.force_build.resolve():
        raise ValueError("force validation refers to a different GAMESS build")
    build_manifest = json.loads(build_manifest_path.read_text())
    if digest(args.force_build / "gamess.force.x") != build_manifest["binary_sha256"]:
        raise ValueError("GAMESS binary changed after instrumented build")
    record = next(row for row in validation["states"] if row["state"] == args.state)
    if not record["passed"]:
        raise ValueError("this state's force validation failed")
    parent = json.loads(args.parent_status.read_text())
    row = next(row for row in parent["results"] if row["repeat"] == 0 and row["iteration"] == 400 and
               row["leg"] == record["leg"] and row["endpoint"] == record["endpoint"])
    source_dir = Path(row["qm"]["log"]).parent
    task = json.loads((source_dir / "input_manifest.json").read_text())
    reporter = MultiStateReporter(task["trajectory"], open_mode="r")
    try:
        thermodynamic = read_endpoint_state(reporter, task["endpoint"])
        system = thermodynamic.get_system(remove_thermostat=True, remove_barostat=True)
        temperature = thermodynamic.temperature.value_in_unit(unit.kelvin)
    finally:
        reporter.close()
    with Dataset(task["trajectory"]) as dataset:
        if int(dataset["states"][task["iteration"], task["replica"]]) != task["endpoint"]:
            raise ValueError("source replica/endpoint mismatch")
        saved = np.asarray(dataset["analysis_particle_indices"][:], dtype=int)
        source_positions = np.asarray(dataset["positions"][task["iteration"], task["replica"]], dtype=float)
        box = np.asarray(dataset["box_vectors"][task["iteration"], task["replica"]], dtype=float)
    if len(saved) != system.getNumParticles() or not np.array_equal(np.sort(saved), np.arange(len(saved))):
        raise ValueError("full saved-particle permutation required")
    positions = np.empty_like(source_positions)
    positions[saved] = source_positions
    indices = physical_indices(task, source_dir, saved)
    full_reference = positions.copy()
    exported = np.load(args.validation.parent / f"{args.state}_gradients.npz")
    shell_reference = exported["coordinates_angstrom"]
    if len(indices) != len(shell_reference):
        raise ValueError("embedded/full-system atom inventory mismatch")
    qm_count = record["n_qm"]
    shell_system, _, _ = build_export_system(next(source_dir.glob("mm_*.inp")))
    shell_integrator = openmm.VerletIntegrator(.001)
    shell_context = openmm.Context(shell_system, shell_integrator, openmm.Platform.getPlatformByName("Reference"))
    for index in reversed(range(system.getNumForces())):
        if isinstance(system.getForce(index), openmm.CMMotionRemover):
            system.removeForce(index)
    integrator = VelocityVerletIntegrator(args.timestep_fs * unit.femtosecond)
    integrator.setConstraintTolerance(1e-8)
    properties = {"DeviceIndex": str(args.device), "Precision": "double"} if args.platform == "CUDA" else {"Threads": "2"}
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName(args.platform), properties)
    context.setPeriodicBoxVectors(*(vector * unit.nanometer for vector in box))
    context.setPositions(positions * unit.nanometer)
    context.applyConstraints(1e-8)
    constrained = context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    initial_projection = float(np.max(np.linalg.norm(constrained - positions, axis=1)) * 10)
    if initial_projection > .001:
        raise ValueError(f"source constraint projection too large: {initial_projection} angstrom")
    positions = np.asarray(constrained).copy()
    masses = np.array([system.getParticleMass(index).value_in_unit(unit.dalton) for index in range(system.getNumParticles())])
    if np.any(masses <= 0):
        raise ValueError("virtual/massless particles require a separately validated proposal")
    output = args.output / args.state
    output.mkdir(parents=True, exist_ok=True)
    protocol_path = output / "protocol.json"
    if (output / "samples.jsonl").exists():
        raise ValueError("refusing to overwrite an existing chain; use a new output directory")
    if protocol_path.exists():
        shutil.copy2(protocol_path, output / f"previous_attempt_{time.time_ns()}.json")
    template = (args.validation.parent / f"efv_{args.state}_base.inp").read_text()
    protocol = {"status": "RUNNING", "state": args.state, "labels_used": False,
                "hamiltonian": "U_MM_PME + lambda*(U_DFTB3_QuanPol_fixed_shell - U_OpenMM_matched_finite_MM_shell)",
                "ensemble": "NVT at the fixed saved periodic box, not NPT",
                "temperature_kelvin": temperature, "box_nm": box.tolist(),
                "fixed_shell_atom_inventory": True, "continuous_fixed_reference_images": True,
                "source_constraints_and_HMR_masses_preserved": True,
                "hybrid_dummy_atoms_retained_in_MM_reference": True,
                "proposal": "full-system MM velocity-Verlet/RATTLE + refreshed constrained Maxwell velocities; exact target Metropolis correction",
                "qm_forces_used_for_propagation": False,
                "embedded_gradient_arrays_and_selected_directions_validated_before_sampling": True,
                "validation_sha256": digest(args.validation), "driver_sha256": digest(Path(__file__)),
                "gamess_binary_sha256": digest(args.force_build / "gamess.force.x"),
                "system_xml_sha256": hashlib.sha256(openmm.XmlSerializer.serialize(system).encode()).hexdigest(),
                "platform": args.platform, "device": args.device, "precision": properties.get("Precision", "CPU"),
                "cuda_nvrtc_preload": os.environ.get("LD_PRELOAD", ""), "shell_mm_platform": "Reference",
                "lambdas": args.lambdas, "moves_per_window": args.moves, "burnin_moves_per_window": args.burnin,
                "md_steps_per_proposal": args.md_steps, "timestep_fs": args.timestep_fs,
                "particles": len(positions), "shell_particles": len(indices), "constraints": system.getNumConstraints(),
                "initial_constraint_projection_max_angstrom": initial_projection,
                "production_tpie_unchanged": True, "finite_shell_not_periodic_qm": True,
                "free_energy_converged": False, "seed": args.seed}
    write_json(protocol_path, protocol)
    (output / "system.xml").write_text(openmm.XmlSerializer.serialize(system))
    np.savez_compressed(output / "initial.npz", positions_nm=positions, box_nm=box, shell_to_system=indices,
                        shell_reference_angstrom=shell_reference, full_reference_nm=full_reference, masses_dalton=masses)
    environment = os.environ.copy()
    environment.update(LD_LIBRARY_PATH=str(Path(sys.executable).resolve().parents[1] / "lib") + ":" + environment.get("LD_LIBRARY_PATH", ""),
                       OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    environment.pop("LD_PRELOAD", None)
    rng = np.random.default_rng(args.seed)
    calls = 0
    start = time.monotonic()

    def delta_energy(current):
        nonlocal calls
        coordinates = shell_reference + 10 * (current[indices] - full_reference[indices])
        deck = output / f"bridge_{args.state}_{calls:05d}.inp"
        deck.write_text(replace_coordinates(template, coordinates, qm_count))
        log = deck.with_suffix(".stdout")
        before = time.monotonic()
        usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)
        with log.open("w") as handle:
            completed = subprocess.run(["timeout", "180", "tcsh", "-f", str(args.force_build / "rungms"), str(deck), "force", "2"],
                                       cwd=output, env=environment, stdout=handle, stderr=subprocess.STDOUT)
        if completed.returncode:
            raise RuntimeError(f"QM bridge evaluation failed at call {calls}")
        result = parse_embedded(log.read_text())
        shell_context.setPositions(coordinates * .1 * unit.nanometer)
        mm_energy = shell_context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        value = result["energy_kcal_mol"] - mm_energy
        if not math.isfinite(value):
            raise ValueError("nonfinite target correction")
        timing = {"call": calls, "qm_kcal_mol": result["energy_kcal_mol"], "mm_shell_kcal_mol": mm_energy,
                  "delta_kcal_mol": value, "wall_seconds": time.monotonic() - before,
                  "deck_sha256": digest(deck), "log_sha256": digest(log)}
        usage_after = resource.getrusage(resource.RUSAGE_CHILDREN)
        timing["qm_child_cpu_seconds"] = usage_after.ru_utime + usage_after.ru_stime - usage_before.ru_utime - usage_before.ru_stime
        write_json(deck.with_suffix(".energy.json"), timing)
        calls += 1
        return value

    def energy_and_kinetic(state):
        velocity = state.getVelocities(asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond)
        kinetic = float(np.sum(.5 * masses[:, None] * np.asarray(velocity)**2))
        return state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole), kinetic

    try:
        composition = validate_composed_target(context, shell_context, full_reference, shell_reference, indices, exported, record, args.validation.parent)
        write_json(output / "target_force_validation.json", composition)
        protocol["composed_periodic_mm_local_qmmm_force_validation_passed"] = composition["passed"]
        if not composition["passed"]:
            raise ValueError("periodic MM plus local QM/MM target force composition failed")
        context.setPositions(positions * unit.nanometer)
        context.setVelocitiesToTemperature(temperature * unit.kelvin, args.seed)
        context.applyVelocityConstraints(1e-8)
        initial = context.getState(getPositions=True, getVelocities=True)
        initial_velocity = initial.getVelocities(asNumpy=True)
        integrator.step(5)
        forward = context.getState(getVelocities=True)
        context.setVelocities(-forward.getVelocities(asNumpy=True))
        integrator.step(5)
        reversed_position = context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        reversal = float(np.max(np.linalg.norm(np.asarray(reversed_position) - positions, axis=1)) * 10)
        protocol["proposal_reversibility_max_angstrom"] = reversal
        if reversal > 1e-4:
            raise ValueError(f"MM proposal reversibility check failed: {reversal} angstrom")
        context.setPositions(positions * unit.nanometer)
        context.setVelocities(initial_velocity)
        current_delta = delta_energy(positions)
        protocol["delta_offset_kcal_mol"] = current_delta
        write_json(protocol_path, protocol)
        for window, coupling in enumerate(args.lambdas):
            accepted_window = 0
            for move in range(args.moves + args.burnin):
                context.setPositions(positions * unit.nanometer)
                context.setVelocitiesToTemperature(temperature * unit.kelvin, int(rng.integers(1, 2**31 - 1)))
                context.applyVelocityConstraints(1e-8)
                before = context.getState(getEnergy=True, getVelocities=True)
                old_u, old_k = energy_and_kinetic(before)
                md_start = time.monotonic()
                integrator.step(args.md_steps)
                proposed = context.getState(getEnergy=True, getVelocities=True, getPositions=True)
                md_wall = time.monotonic() - md_start
                new_u, new_k = energy_and_kinetic(proposed)
                proposal_positions = np.asarray(proposed.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
                proposal_delta = delta_energy(proposal_positions)
                log_alpha = log_acceptance(new_u - old_u, new_k - old_k, proposal_delta - current_delta, coupling, temperature)
                accepted = math.log(rng.uniform()) < min(0., log_alpha)
                previous_delta = current_delta
                rms_move = float(np.sqrt(np.mean(np.sum((proposal_positions[indices[:qm_count]] - positions[indices[:qm_count]])**2, axis=1))) * 10)
                if accepted:
                    positions, current_delta = proposal_positions.copy(), proposal_delta
                    accepted_window += 1
                else:
                    context.setPositions(positions * unit.nanometer)
                sample = {"window": window, "lambda": coupling, "move": move, "burnin": move < args.burnin,
                          "accepted": accepted, "log_acceptance": log_alpha,
                          "mm_integration_energy_error_kj_mol": new_u - old_u + new_k - old_k,
                          "previous_delta_kcal_mol": previous_delta, "proposal_delta_kcal_mol": proposal_delta,
                          "current_delta_kcal_mol": current_delta, "md_wall_seconds": md_wall,
                          "proposal_ligand_rms_move_angstrom": rms_move,
                          "ligand_rms_displacement_from_start_angstrom": float(np.sqrt(np.mean(np.sum((positions[indices[:qm_count]] - full_reference[indices[:qm_count]])**2, axis=1))) * 10)}
                with (output / "samples.jsonl").open("a") as handle:
                    handle.write(json.dumps(sample) + "\n")
                if move % args.save_every == 0 or move == args.moves + args.burnin - 1:
                    np.savez_compressed(output / f"frame_w{window:02d}_m{move:04d}.npz", positions_nm=positions, box_nm=box,
                                        delta_kcal_mol=current_delta, coupling=coupling)
                    (output / "latest.chk").write_bytes(context.createCheckpoint())
                protocol.update(current_window=window, current_lambda=coupling, completed_move=move + 1,
                                qm_evaluations=calls, elapsed_seconds=time.monotonic() - start)
                write_json(protocol_path, protocol)
                print(json.dumps({"state": args.state, "lambda": coupling, "move": move, "accepted": accepted,
                                  "delta_kcal_mol": current_delta}), flush=True)
            protocol.setdefault("windows", []).append({"lambda": coupling, "accepted": accepted_window,
                                                        "attempted": args.moves + args.burnin,
                                                        "acceptance_fraction": accepted_window / (args.moves + args.burnin)})
        protocol["status"] = "PILOT_COMPLETE_NOT_CONVERGED_FREE_ENERGY"
        protocol["elapsed_seconds"] = time.monotonic() - start
    except Exception as error:
        protocol.update(status="ERROR", error=repr(error), qm_evaluations=calls, elapsed_seconds=time.monotonic() - start)
        write_json(protocol_path, protocol)
        raise
    finally:
        del context, integrator, shell_context, shell_integrator
    write_json(protocol_path, protocol)
    return protocol


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation", type=Path, required=True)
    parser.add_argument("--parent-status", type=Path, required=True)
    parser.add_argument("--force-build", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--state", choices=("b0", "b10", "s0", "s10"))
    parser.add_argument("--lambdas", nargs="+", type=float, default=[0., .05, .1])
    parser.add_argument("--moves", type=int, default=12)
    parser.add_argument("--burnin", type=int, default=4)
    parser.add_argument("--md-steps", type=int, default=100)
    parser.add_argument("--timestep-fs", type=float, default=.5)
    parser.add_argument("--save-every", type=int, default=4)
    parser.add_argument("--platform", choices=("CPU", "CUDA"), default="CUDA")
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--seed", type=int, default=274101)
    args = parser.parse_args()
    for field in ("validation", "parent_status", "force_build", "output"):
        setattr(args, field, getattr(args, field).resolve())
    if not args.lambdas or args.lambdas[0] != 0 or any(not math.isfinite(value) or not 0 <= value <= 1 for value in args.lambdas) or any(a >= b for a, b in zip(args.lambdas, args.lambdas[1:])):
        parser.error("strictly increasing finite bridge lambdas starting at zero required")
    if min(args.moves, args.md_steps, args.save_every) < 1 or args.burnin < 0 or not math.isfinite(args.timestep_fs) or not 0 < args.timestep_fs <= 1:
        parser.error("positive counts, nonnegative burnin and timestep in (0,1] fs required")
    args.output.mkdir(parents=True, exist_ok=True)
    if args.state:
        try:
            run_state(args)
        except Exception as error:
            path = args.output / args.state / "protocol.json"
            failure = json.loads(path.read_text()) if path.exists() else {"state": args.state}
            failure.update(status="ERROR", free_energy_converged=False, error=repr(error))
            write_json(path, failure)
            raise
        return
    if json.loads(args.validation.read_text()).get("status") != "PASSED":
        raise ValueError("full force validation not passed; sampling blocked")
    status = {"status": "RUNNING", "labels_used": False, "free_energy_converged": False,
              "lambdas": args.lambdas, "states": {}, "production_tpie_modified": False}
    write_json(args.output / "status.json", status)

    def worker(device):
        rows = []
        for state in (("b0", "s0") if device == 0 else ("b10", "s10")):
            command = [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--state", state,
                       "--device", str(device), "--seed", str(args.seed + (0 if state == "b0" else 1 if state == "b10" else 2 if state == "s0" else 3))]
            with (args.output / f"{state}.stdout").open("w") as handle:
                result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
            path = args.output / state / "protocol.json"
            rows.append((state, json.loads(path.read_text()) if path.exists() else {"status": "ERROR", "returncode": result.returncode}))
            if result.returncode:
                break
        return rows

    with ThreadPoolExecutor(max_workers=2) as pool:
        for rows in pool.map(worker, (0, 1)):
            for state, record in rows:
                status["states"][state] = record
            write_json(args.output / "status.json", status)
    complete = len(status["states"]) == 4 and all(row["status"] == "PILOT_COMPLETE_NOT_CONVERGED_FREE_ENERGY" for row in status["states"].values())
    status["status"] = "PILOT_COMPLETE_NOT_CONVERGED_FREE_ENERGY" if complete else "ERROR"
    write_json(args.output / "status.json", status)
    if not complete:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
