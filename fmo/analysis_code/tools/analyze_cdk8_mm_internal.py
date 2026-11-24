"""Ligand-internal MM geometry sensitivity, not a binding-energy score."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import openmm
from openmm import app, unit
from openmmforcefields.generators import SystemGenerator
from openff.toolkit import Molecule, Topology
from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
from rdkit import Chem

from run_cdk8_solution_cycle import digest, group, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--bound-run", type=Path, required=True)
    args = parser.parse_args()
    workspace, reference, bound_run = args.workspace.resolve(), args.reference.resolve(), args.bound_run.resolve()
    sdf = workspace / "vendor/public_binding_free_energy_benchmark/fep_benchmark_inputs/structure_inputs/misc/cdk8_koehler_ligands.sdf"
    output = bound_run / "mm_internal"
    output.mkdir(exist_ok=True)
    write_json(output / "plan.json", {"driver_sha256": digest(__file__), "sdf_sha256": digest(sdf),
        "forcefield": "openff-2.2.0", "charges": "openff-gnn-am1bcc-1.0.0.pt",
        "labels_used": False, "warning": "Ligand-only vacuum MM single points; diagnostic, not solution energy or binding score"})
    results = []
    for molecule in Chem.SDMolSupplier(str(sdf), removeHs=False):
        name = molecule.GetProp("_Name")
        off_ligand = Molecule.from_rdkit(molecule, allow_undefined_stereo=False)
        off_ligand.assign_partial_charges(partial_charge_method="openff-gnn-am1bcc-1.0.0.pt", toolkit_registry=NAGLToolkitWrapper())
        topology = Topology.from_molecules([off_ligand]).to_openmm()
        generator = SystemGenerator(forcefields=[], small_molecule_forcefield="openff-2.2.0",
            forcefield_kwargs={"constraints": None, "rigidWater": False},
            nonperiodic_forcefield_kwargs={"nonbondedMethod": app.NoCutoff}, molecules=[off_ligand], cache=None)
        system = generator.create_system(topology)
        groups = {}
        for index, force in enumerate(system.getForces()):
            force.setForceGroup(index)
            groups[index] = type(force).__name__
        integrator = openmm.VerletIntegrator(0.001)
        context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName("Reference"))
        deck = (reference / f"pcm/{name}/ligand/input.inp").read_text()
        rows = group(deck, "FMOXYZ")[1].strip().splitlines()
        if [row.split()[0].upper() for row in rows] != [atom.GetSymbol().upper() for atom in molecule.GetAtoms()]:
            raise ValueError("ligand atom-order mismatch")
        geometries = {"input": molecule.GetConformer().GetPositions(),
            "k10": np.array([row.split()[2:] for row in rows], dtype=float),
            "qm_free": np.array([row.split()[2:] for row in json.loads((reference / f"optimized/{name}/ligand/result.json").read_text())["optimized_coordinates"]], dtype=float)}
        summary = bound_run / name / "summary.json"
        if summary.exists():
            bound = json.loads(summary.read_text())
            if bound["status"] == "SUCCEEDED":
                geometries["qm_bound"] = np.array([row.split()[2:] for row in bound["optimized_coordinates"][-len(rows):]], dtype=float)
        energies = {}
        for state, coordinates in geometries.items():
            context.setPositions(coordinates / 10 * unit.nanometer)
            components = {label: context.getState(getEnergy=True, groups=1 << index).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole) for index,label in groups.items()}
            energies[state] = {"total_kcal_mol": context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole), "components": components}
        del context, integrator
        row = {"ligand": name, "energies": energies, "chemical_state": "same explicit-H input state; same molecule/charges for all geometries"}
        results.append(row)
        write_json(output / "evaluation.json", {"rows": results, "n": len(results), "warning": "Within-molecule energy differences only; never compare absolute ligand totals across different molecules"})
        print(name, energies, flush=True)


if __name__ == "__main__":
    main()
