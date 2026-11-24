"""Rebuild exactly the exported finite FFDATA MM model in OpenMM for parity."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from report_quanpol_endpoint_pilot import parse
from run_matched_openfe import digest, write_json


def read_sections(text):
    lines = text.split(" $FFDATA", 1)[1].strip().splitlines()
    result = {}
    current = None
    for line in lines:
        words = line.split()
        if not words:
            continue
        if words[0] == "$END":
            break
        if words[0] == "STOP":
            current = None
        elif current is None:
            current = words[0]
            result[current] = []
        else:
            result[current].append(words)
    return result


def build_export_system(path):
    import openmm
    from openmm import unit

    sections = read_sections(path.read_text())
    system = openmm.System()
    nonbonded = openmm.NonbondedForce()
    bonds = openmm.HarmonicBondForce()
    angles = openmm.HarmonicAngleForce()
    torsions = openmm.PeriodicTorsionForce()
    for row in sections["PARAMETERS"]:
        _, mass, charge, _, sigma, epsilon, _, _ = row
        system.addParticle(float(mass))
        nonbonded.addParticle(float(charge), float(sigma)*.1, float(epsilon)*4.184)
    pairs = []
    for row in sections.get("BOND", []):
        _, first, second, constant, length = row
        first, second = int(first)-1, int(second)-1
        pairs.append((first, second))
        bonds.addBond(first, second, float(length)*.1, float(constant)*2*418.4)
    for row in sections.get("ANGLE", []):
        _, first, middle, last, constant, angle = row
        angles.addAngle(int(first)-1, int(middle)-1, int(last)-1,
                        np.deg2rad(float(angle)), float(constant)*2*4.184)
    for row in sections.get("DIHROT", []):
        _, first, second, third, fourth, constant, periodicity, phase = row
        torsions.addTorsion(int(first)-1, int(second)-1, int(third)-1, int(fourth)-1,
                            int(periodicity), np.deg2rad(float(phase)), float(constant)*4.184)
    nonbonded.createExceptionsFromBonds(pairs, 1/1.2, .5)
    names = ("bond_stretching_kcal_mol", "bond_angle_bending_kcal_mol",
             "dihedral_rotation_kcal_mol", "nonbonded_kcal_mol")
    for group, force in enumerate((bonds, angles, torsions, nonbonded)):
        force.setForceGroup(group)
        system.addForce(force)
    coordinates = np.array([[float(value) for value in row[2:5]] for row in sections["COORDINATES"]])
    return system, coordinates, names


def evaluate_export(path, include_forces=False):
    import openmm
    from openmm import unit

    system, coordinates, names = build_export_system(path)
    integrator = openmm.VerletIntegrator(.001)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName("CPU"), {"Threads": "2"})
    context.setPositions(coordinates*.1)
    result = {name: context.getState(getEnergy=True, groups=1 << group).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
              for group, name in enumerate(names)}
    result["potential_kcal_mol"] = sum(result.values())
    if include_forces:
        result["forces_kcal_mol_angstrom"] = context.getState(getForces=True).getForces(asNumpy=True).value_in_unit(unit.kilocalorie_per_mole / unit.angstrom).tolist()
    del context, integrator
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    results = []
    for directory in sorted(args.root.glob("*/endpoint_*")):
        deck = next(directory.glob("mm_only_*.inp"))
        reference = evaluate_export(deck)
        quanpol = parse(directory / "mm_matched.stdout")
        difference = quanpol["potential_kcal_mol"] - reference["potential_kcal_mol"]
        result = {"directory": str(directory), "deck_sha256": digest(deck), "openmm_export": reference,
                  "quanpol": quanpol, "difference_kcal_mol": difference,
                  "passed": abs(difference) < .1 and quanpol["status"] == "SUCCEEDED"}
        write_json(directory / "mm_parity.json", result)
        results.append(result)
        print(directory, difference, flush=True)
    write_json(args.root / "mm_parity.json", {"results": results, "passed": all(row["passed"] for row in results),
        "scope": "Reconstructed FFDATA model only; this does not establish equivalence to the periodic OpenFE system."})


if __name__ == "__main__":
    main()
