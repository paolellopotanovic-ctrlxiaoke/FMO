"""Matched ligand-internal MM diagnostics; never an affinity scoring term."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import multiprocessing
from pathlib import Path
import time

import numpy as np
from rdkit import Chem

from fmotpiefinal.core.bound_geometry import load_bound_ligand
from run_cdk8_solution_cycle import digest, write_json
from run_local_sidechain_transfer import directory_key


TARGETS = ("CHK1", "BACE", "CDK8", "HIV1PROT", "Galectin3")


def execute(payload):
    target, name, cohort, source, output = payload
    destination = Path(output) / target / directory_key(name) / "result.json"
    try:
        import openmm
        from openmm import app, unit
        from openmmforcefields.generators import SystemGenerator
        from openff.toolkit import Molecule, Topology
        from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
        matches = [molecule for molecule in Chem.SDMolSupplier(cohort["sdf"], removeHs=False)
            if molecule is not None and molecule.GetProp("_Name").strip() == name]
        if len(matches) != 1:
            raise ValueError("one matching explicit-H molecule required")
        molecule = matches[0]
        Chem.AssignStereochemistryFrom3D(molecule)
        off = Molecule.from_rdkit(molecule, allow_undefined_stereo=False)
        off.assign_partial_charges(partial_charge_method="openff-gnn-am1bcc-1.0.0.pt",
            toolkit_registry=NAGLToolkitWrapper())
        topology = Topology.from_molecules([off]).to_openmm()
        generator = SystemGenerator(forcefields=[], small_molecule_forcefield="openff-2.2.0",
            forcefield_kwargs={"constraints": None, "rigidWater": False},
            nonperiodic_forcefield_kwargs={"nonbondedMethod": app.NoCutoff}, molecules=[off], cache=None)
        system = generator.create_system(topology)
        force_groups = {}
        for index, force in enumerate(system.getForces()):
            force.setForceGroup(index)
            force_groups[index] = f"{index}:{type(force).__name__}"
        integrator = openmm.VerletIntegrator(0.001)
        context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName("Reference"))
        candidate_path = Path(source) / target / directory_key(name) / "minimization/minimized_complex.pdb"
        if not cohort["local_residues"]:
            candidate_path = Path(cohort["complexes"][name])
        geometries = {"input": np.asarray(molecule.GetConformer().GetPositions()),
            "baseline_k10": np.asarray([atom.xyz for atom in load_bound_ligand(cohort["sdf"], name, cohort["complexes"][name]).atoms]),
            "sidechain_B": np.asarray([atom.xyz for atom in load_bound_ligand(cohort["sdf"], name, candidate_path).atoms])}
        energies = {}

        def energy():
            return {"total_kcal_mol": context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole),
                "components_kcal_mol": {label: context.getState(getEnergy=True, groups=1 << index).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
                    for index, label in force_groups.items()}}

        for state, positions in geometries.items():
            if positions.shape != (molecule.GetNumAtoms(), 3) or not np.isfinite(positions).all():
                raise ValueError("matching finite geometry required")
            context.setPositions(positions / 10 * unit.nanometer)
            energies[state] = energy()
        context.setPositions(geometries["input"] / 10 * unit.nanometer)
        openmm.LocalEnergyMinimizer.minimize(context, 1.0 * unit.kilojoule_per_mole / unit.nanometer, 1000)
        energies["isolated_local_minimum"] = energy()
        forces = context.getState(getForces=True).getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
        rms_force = float(np.sqrt(np.mean(np.sum(forces**2, axis=1))))
        serialized = openmm.XmlSerializer.serialize(system)
        del context, integrator
        reference = energies["isolated_local_minimum"]["total_kcal_mol"]
        result = {"status": "SUCCEEDED", "energies": energies,
            "relative_to_input_start_local_minimum_kcal_mol": {state: row["total_kcal_mol"]-reference for state, row in energies.items()},
            "isolated_local_minimum_rms_force_kj_mol_nm": rms_force,
            "isolated_local_minimum_force_tolerance_passed": rms_force <= 1.0,
            "system_sha256": hashlib.sha256(serialized.encode()).hexdigest(),
            "geometry_sha256": {cohort["complexes"][name]: digest(cohort["complexes"][name]), str(candidate_path): digest(candidate_path)},
            "warning": "Within-molecule vacuum MM geometry sensitivity; not a solution free energy or binding score. One input-start local minimum, not a global or ensemble reference."}
    except Exception as error:
        result = {"status": "FAILED", "error": f"{type(error).__name__}: {error}"}
    write_json(destination, result)
    return target, name, result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("workers must be positive")
    source, output = args.source.resolve(), args.output.resolve()
    original = json.loads((source / "frozen_protocol.json").read_text())
    cohorts = {target: original["cohorts"][target] for target in TARGETS}
    hashes = {}
    for target, cohort in cohorts.items():
        for name in cohort["names"]:
            paths = [Path(cohort["sdf"]), Path(cohort["complexes"][name]),
                source / target / directory_key(name) / "minimization/minimized_complex.pdb"]
            if not cohort["local_residues"]:
                paths[-1] = Path(cohort["complexes"][name])
            for path in paths:
                hashes[str(path)] = digest(path)
    plan = {"schema": "weak-cohort-ligand-internal-mm-v1", "targets": list(TARGETS),
        "source_protocol_sha256": digest(source / "frozen_protocol.json"), "source_sha256": hashes,
        "driver_sha256": digest(__file__), "forcefield": "OpenFF 2.2.0", "charges": "NAGL openff-gnn-am1bcc-1.0.0.pt",
        "platform": "Reference", "isolation": "vacuum, NoCutoff, no restraints",
        "local_minimization": {"start": "input", "max_iterations": 1000, "tolerance_kj_mol_nm": 1.0},
        "labels_used": False, "prediction_score_changed": False}
    frozen = output / "frozen_protocol.json"
    if frozen.exists() and json.loads(frozen.read_text()) != plan:
        raise ValueError("refusing changed MM diagnostic protocol")
    write_json(frozen, plan)
    started = time.monotonic()
    rows = {target: {} for target in TARGETS}
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context("spawn")) as pool:
        futures = []
        for target, cohort in cohorts.items():
            for name in cohort["names"]:
                futures.append(pool.submit(execute, (target, name, cohort, str(source), str(output))))
        for future in as_completed(futures):
            target, name, result = future.result()
            rows[target][name] = result
            write_json(output / "progress.json", rows)
            print(target, name, result["status"], result.get("relative_to_input_start_local_minimum_kcal_mol", result.get("error")), flush=True)
    failures = {target: {name: row for name, row in group.items() if row["status"] != "SUCCEEDED"}
        for target, group in rows.items()}
    result = {"status": "COMPLETE" if not any(failures.values()) else "INCOMPLETE_WITH_FAILURES",
        "n_requested": sum(len(cohort["names"]) for cohort in cohorts.values()),
        "n_succeeded": sum(row["status"] == "SUCCEEDED" for group in rows.values() for row in group.values()),
        "failures": failures, "wall_seconds": time.monotonic()-started}
    write_json(output / "final_status.json", result)
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
