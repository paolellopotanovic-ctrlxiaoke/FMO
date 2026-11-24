"""Evaluate the original OpenMM endpoint Hamiltonian on the QM/MM pilot frames."""
from __future__ import annotations

import argparse
import bz2
import json
from pathlib import Path

import numpy as np

from run_matched_openfe import write_json


def endpoint_parameters(system, endpoint):
    value = float(endpoint == 10)
    names = set()
    for force in system.getForces():
        if hasattr(force, "getNumGlobalParameters"):
            for index in range(force.getNumGlobalParameters()):
                name = force.getGlobalParameterName(index)
                if name.startswith("lambda_"):
                    names.add(name)
    return {name: (value if name.startswith("lambda_") else None) for name in sorted(names)}


def evaluate_leg(root: Path, endpoint: int, iteration: int):
    import openmm
    from netCDF4 import Dataset

    shared = sorted((root / "shared").glob("*/hybrid_system.xml.bz2"))
    topology = sorted((root / "shared").glob("*/hybrid_system.pdb"))
    trajectories = sorted((root / "shared").glob("*/simulation.nc"))
    if len(shared) != 3 or len(topology) != 3 or len(trajectories) != 3:
        raise ValueError("expected three complete repeats")
    system = openmm.XmlSerializer.deserialize(bz2.decompress(shared[0].read_bytes()).decode())
    integrator = openmm.VerletIntegrator(0.001 * openmm.unit.picoseconds)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName("Reference"))
    parameters = endpoint_parameters(system, endpoint)
    for name, value in parameters.items():
        context.setParameter(name, value)
    rows = []
    for trajectory in trajectories:
        with Dataset(trajectory) as data:
            states = np.asarray(data["states"][iteration])
            replicas = np.flatnonzero(states == endpoint)
            if len(replicas) != 1:
                raise ValueError("endpoint does not have exactly one replica")
            replica = int(replicas[0])
            raw = np.asarray(data["positions"][iteration, replica], dtype=float)
            saved = np.asarray(data["analysis_particle_indices"][:], dtype=int)
            coordinates = np.empty((len(saved), 3), dtype=float)
            coordinates[saved] = raw * 10.0
            box = np.asarray(data["box_vectors"][iteration, replica], dtype=float) * 10.0
        context.setPeriodicBoxVectors(*[openmm.Vec3(*row) * 0.1 for row in box] * openmm.unit.nanometer)
        context.setPositions(coordinates * 0.1 * openmm.unit.nanometer)
        state = context.getState(getEnergy=True)
        rows.append({"trajectory": str(trajectory), "replica": replica,
                     "energy_kcal_mol": state.getPotentialEnergy().value_in_unit(openmm.unit.kilocalorie_per_mole),
                     "parameters": parameters})
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--target", default="CDK8")
    parser.add_argument("--iteration", type=int, default=2000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = {"schema": "same-frame-openmm-endpoint-energy-v1", "target": args.target,
              "iteration": args.iteration, "legs": {}}
    for leg in ("bound", "solution"):
        root = args.run / "production" / args.target / "edge_000" / leg
        result["legs"][leg] = {str(endpoint): evaluate_leg(root, endpoint, args.iteration)
                                for endpoint in (0, 10)}
    means = {}
    for leg in result["legs"]:
        means[leg] = {endpoint: float(np.mean([row["energy_kcal_mol"] for row in result["legs"][leg][endpoint]]))
                      for endpoint in ("0", "10")}
    result["mean_energy_kcal_mol"] = means
    result["endpoint_binding_difference_kcal_mol"] = {
        endpoint: means["bound"][endpoint] - means["solution"][endpoint] for endpoint in ("0", "10")}
    result["endpoint_double_difference_kcal_mol"] = (result["endpoint_binding_difference_kcal_mol"]["10"] -
                                                       result["endpoint_binding_difference_kcal_mol"]["0"])
    result["interpretation"] = [
        "These are the original OpenMM endpoint potential energies on the same saved physical frames.",
        "They are a Hamiltonian audit, not an affinity or free-energy estimate.",
        "The endpoint energy is only comparable to QuanPol after matching the selected environment and force terms.",
    ]
    write_json(args.output, result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
