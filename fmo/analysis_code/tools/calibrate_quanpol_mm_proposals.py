"""Measure reference-MM HMC discretization errors without labels or new QM calls."""
import argparse
import json
from pathlib import Path
import time

import numpy as np
import openmm
from openmm import unit
from openmmtools.integrators import VelocityVerletIntegrator

from run_matched_openfe import write_json
from sample_quanpol_bridge import GAS_KJ


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--states", nargs="+", default=["b0", "b10"])
    parser.add_argument("--seeds", type=int, default=8)
    parser.add_argument("--device", type=int, default=0)
    args = parser.parse_args()
    report = {"status": "RUNNING", "labels_used": False, "qm_evaluations": 0, "current_pilot_modified": False, "states": {}}
    for name in args.states:
        directory = args.root / name
        protocol = json.loads((directory / "protocol.json").read_text())
        data = np.load(directory / "initial.npz")
        positions = data["positions_nm"]
        masses = data["masses_dalton"]
        temperature = protocol["temperature_kelvin"]
        system = openmm.XmlSerializer.deserialize((directory / "system.xml").read_text())
        rows = []
        for timestep in (.5, .25, .1):
            integrator = VelocityVerletIntegrator(timestep * unit.femtosecond)
            integrator.setConstraintTolerance(1e-8)
            context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName("CUDA"), {"Precision": "double", "DeviceIndex": str(args.device)})
            context.setPeriodicBoxVectors(*data["box_nm"])
            started = time.monotonic()
            for seed in range(args.seeds):
                context.setPositions(positions * unit.nanometer)
                context.setVelocitiesToTemperature(temperature * unit.kelvin, 81403 + seed)
                context.applyVelocityConstraints(1e-8)
                initial = context.getState(getVelocities=True, getEnergy=True)
                velocity = np.asarray(initial.getVelocities(asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond))
                old = initial.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole) + float(np.sum(.5 * masses[:, None] * velocity**2))
                steps = round(50 / timestep)
                integrator.step(steps)
                final = context.getState(getVelocities=True, getEnergy=True)
                velocity = np.asarray(final.getVelocities(asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond))
                new = final.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole) + float(np.sum(.5 * masses[:, None] * velocity**2))
                difference = new - old
                rows.append({"timestep_fs": timestep, "steps": steps, "proposal_fs": 50, "seed": 81403 + seed,
                             "mm_hamiltonian_error_kj_mol": difference,
                             "reference_acceptance_probability": float(np.exp(min(0., -difference / (GAS_KJ * temperature))))})
            elapsed = time.monotonic() - started
            selected = [row for row in rows if row["timestep_fs"] == timestep]
            print(json.dumps({"state": name, "timestep_fs": timestep,
                              "mean_acceptance_probability": float(np.mean([row["reference_acceptance_probability"] for row in selected])),
                              "rms_energy_error_kj_mol": float(np.sqrt(np.mean([row["mm_hamiltonian_error_kj_mol"]**2 for row in selected]))),
                              "wall_seconds": elapsed}), flush=True)
            del context, integrator
        report["states"][name] = rows
        write_json(args.root / "mm_proposal_calibration.json", report)
    report["status"] = "COMPLETE_DIAGNOSTIC_ONLY"
    write_json(args.root / "mm_proposal_calibration.json", report)


if __name__ == "__main__":
    main()
