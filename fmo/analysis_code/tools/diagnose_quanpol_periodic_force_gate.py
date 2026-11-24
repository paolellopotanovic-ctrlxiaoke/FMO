"""Resolve periodic MM finite-difference failures by force group and step size."""
import argparse
import json
import itertools
import math
from pathlib import Path

import numpy as np
import openmm
from openmm import unit

from run_matched_openfe import write_json
from validate_quanpol_embedded_forces import direction_vector


def custom_pair_cutoff_energy(system, source_context, cutoff):
    bundles = []
    cache = {}
    for original in system.getForces():
        if not isinstance(original, openmm.CustomNonbondedForce) or original.getNonbondedMethod() == openmm.CustomNonbondedForce.NoCutoff or original.getUseSwitchingFunction():
            continue
        if abs(original.getCutoffDistance().value_in_unit(unit.nanometer) - cutoff) > 1e-10:
            raise ValueError("different custom/native cutoff radii need separate boundary audits")
        pair_system = openmm.System()
        pair_system.addParticle(1)
        pair_system.addParticle(1)
        pair_force = openmm.CustomNonbondedForce(original.getEnergyFunction())
        for index in range(original.getNumPerParticleParameters()):
            pair_force.addPerParticleParameter(original.getPerParticleParameterName(index))
        for index in range(original.getNumGlobalParameters()):
            name = original.getGlobalParameterName(index)
            pair_force.addGlobalParameter(name, source_context.getParameter(name))
        if hasattr(original, "getNumComputedValues"):
            for index in range(original.getNumComputedValues()):
                pair_force.addComputedValue(*original.getComputedValueParameters(index))
        for index in range(original.getNumTabulatedFunctions()):
            pair_force.addTabulatedFunction(original.getTabulatedFunctionName(index), openmm.XmlSerializer.deserialize(openmm.XmlSerializer.serialize(original.getTabulatedFunction(index))))
        pair_force.addParticle(original.getParticleParameters(0))
        pair_force.addParticle(original.getParticleParameters(0))
        pair_system.addForce(pair_force)
        integrator = openmm.VerletIntegrator(.001)
        context = openmm.Context(pair_system, integrator, openmm.Platform.getPlatformByName("Reference"))
        context.setPositions([[0, 0, 0], [cutoff, 0, 0]])
        groups = [tuple(set(group) for group in original.getInteractionGroupParameters(index)) for index in range(original.getNumInteractionGroups())]
        exclusions = {tuple(sorted(original.getExclusionParticles(index))) for index in range(original.getNumExclusions())}
        bundles.append((original, pair_force, context, integrator, groups, exclusions))

    def evaluate(first, second):
        pair = tuple(sorted((int(first), int(second))))
        if pair in cache:
            return cache[pair]
        value = 0.
        for original, pair_force, context, integrator, groups, exclusions in bundles:
            if pair in exclusions or groups and not any((first in group[0] and second in group[1]) or (second in group[0] and first in group[1]) for group in groups):
                continue
            pair_force.setParticleParameters(0, original.getParticleParameters(int(first)))
            pair_force.setParticleParameters(1, original.getParticleParameters(int(second)))
            pair_force.updateParametersInContext(context)
            value += context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        cache[pair] = value
        return value

    return evaluate


def evaluate_cutoff_jumps(prepared, positions, indices, vector, steps):
    displacement = np.zeros_like(positions)
    displacement[indices] = vector
    first, second, separation, jump, cutoff = (prepared[key] for key in ("first", "second", "separation", "jump", "cutoff"))
    relative = displacement[first] - displacement[second]
    rows = []
    for step in steps:
        if 2 * step * .1 * np.linalg.norm(vector, axis=1).max() > .002:
            raise ValueError("displacement exceeds prepared cutoff audit neighborhood")
        minus = np.linalg.norm(separation - step * .1 * relative, axis=1) < cutoff
        plus = np.linalg.norm(separation + step * .1 * relative, axis=1) < cutoff
        crossing = plus.astype(int) - minus.astype(int)
        nonzero = np.flatnonzero(crossing)
        effective_jump = jump.copy()
        for index in nonzero:
            if "custom_pair_energy" in prepared:
                effective_jump[index] += prepared["custom_pair_energy"](int(first[index]), int(second[index]))
        rows.append({"step_angstrom": step, "crossing_pairs": len(nonzero),
                     "predicted_gradient_discontinuity": float(np.sum(effective_jump * crossing) / (2 * step)),
                     "pairs": [{"first": int(first[index]), "second": int(second[index]), "crossing": int(crossing[index]),
                                "jump_kcal_mol": float(effective_jump[index]), "native_jump_kcal_mol": float(jump[index]),
                                "custom_jump_kcal_mol": float(effective_jump[index] - jump[index]),
                                "distance_nm": float(np.linalg.norm(separation[index]))} for index in nonzero]})
    return rows


def cutoff_jump_diagnostics(system, context, positions, box, indices, vector, steps, cache=None):
    from scipy.spatial import cKDTree

    if cache:
        return evaluate_cutoff_jumps(cache, positions, indices, vector, steps)
    force = next(force for force in system.getForces() if isinstance(force, openmm.NonbondedForce))
    if force.getUseSwitchingFunction():
        raise ValueError("diagnostic currently expects the saved unswitched PME reference")
    cutoff = force.getCutoffDistance().value_in_unit(unit.nanometer)
    alpha, *_ = force.getPMEParametersInContext(context)
    alpha = float(alpha.value_in_unit(unit.nanometer**-1)) if hasattr(alpha, "value_in_unit") else float(alpha)
    parameters = np.array([[charge.value_in_unit(unit.elementary_charge), sigma.value_in_unit(unit.nanometer), epsilon.value_in_unit(unit.kilojoule_per_mole)]
                           for charge, sigma, epsilon in (force.getParticleParameters(index) for index in range(force.getNumParticles()))])
    for offset in range(force.getNumParticleParameterOffsets()):
        name, index, charge, sigma, epsilon = force.getParticleParameterOffset(offset)
        parameters[index] += context.getParameter(name) * np.array([float(charge), float(sigma), float(epsilon)])
    exceptions = {tuple(sorted(force.getExceptionParameters(offset)[:2])) for offset in range(force.getNumExceptions())}
    fractional = positions @ np.linalg.inv(box)
    wrapped = (fractional - np.floor(fractional)) @ box
    lattice = np.array(list(itertools.product((-1, 0, 1), repeat=3))) @ box
    images = (wrapped[None, :, :] + lattice[:, None, :]).reshape(-1, 3)
    tree = cKDTree(images)
    width = .002
    candidates = tree.query_ball_point(wrapped[indices], cutoff + width)
    first = np.repeat(indices, [len(row) for row in candidates])
    image_index = np.concatenate(candidates)
    second = image_index % len(positions)
    separation = wrapped[first] - images[image_index]
    distance = np.linalg.norm(separation, axis=1)
    mask = np.abs(distance - cutoff) < width
    first, second, separation = first[mask], second[mask], separation[mask]
    mobile = set(indices.tolist())
    keep = np.array([not (int(other) in mobile and int(atom) > int(other)) and tuple(sorted((int(atom), int(other)))) not in exceptions
                     for atom, other in zip(first, second)], dtype=bool)
    first, second, separation = first[keep], second[keep], separation[keep]
    qprod = parameters[first, 0] * parameters[second, 0]
    sigma = .5 * (parameters[first, 1] + parameters[second, 1])
    epsilon = np.sqrt(parameters[first, 2] * parameters[second, 2])
    ratio = (sigma / cutoff)**6
    jump = (138.93545764438198 * qprod * math.erfc(alpha * cutoff) / cutoff + 4 * epsilon * ratio * (ratio - 1)) / 4.184
    prepared = {"first": first, "second": second, "separation": separation, "jump": jump, "cutoff": cutoff,
                "custom_pair_energy": custom_pair_cutoff_energy(system, context, cutoff)}
    if cache is not None:
        cache.update(prepared)
    return evaluate_cutoff_jumps(prepared, positions, indices, vector, steps)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chain", type=Path, required=True)
    parser.add_argument("--validation", type=Path, required=True)
    args = parser.parse_args()
    data = np.load(args.chain / "initial.npz")
    positions = data["full_reference_nm"]
    indices = data["shell_to_system"]
    record = next(row for row in json.loads(args.validation.read_text())["states"] if row["state"] == args.chain.name)
    system = openmm.XmlSerializer.deserialize((args.chain / "system.xml").read_text())
    settings = []
    for index, force in enumerate(system.getForces()):
        force.setForceGroup(index)
        row = {"index": index, "name": force.getName(), "type": type(force).__name__}
        if hasattr(force, "getCutoffDistance"):
            row["cutoff_nm"] = force.getCutoffDistance().value_in_unit(unit.nanometer)
        if hasattr(force, "getUseSwitchingFunction"):
            row["switching"] = force.getUseSwitchingFunction()
        if hasattr(force, "getSwitchingDistance"):
            row["switching_nm"] = force.getSwitchingDistance().value_in_unit(unit.nanometer)
        settings.append(row)
    reciprocal_group = system.getNumForces()
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            force.setReciprocalSpaceForceGroup(reciprocal_group)
    settings.append({"index": reciprocal_group, "name": "PME_reciprocal", "type": "reciprocal"})
    integrator = openmm.VerletIntegrator(.001)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName("CUDA"), {"Precision": "double", "DeviceIndex": "0"})
    context.setPeriodicBoxVectors(*data["box_nm"])
    context.setPositions(positions)
    report = {"settings": settings, "directions": []}
    for ordinal in (4, 5):
        direction = record["finite_differences"][::2][ordinal]
        vector = direction_vector(direction, data["shell_reference_angstrom"].shape, record["n_qm"])
        analytic = []
        for index in range(len(settings)):
            forces = context.getState(getForces=True, groups=1 << index).getForces(asNumpy=True).value_in_unit(unit.kilocalorie_per_mole / unit.angstrom)
            analytic.append(float(-np.sum(np.asarray(forces)[indices] * vector)))
        rows = []
        cutoff = cutoff_jump_diagnostics(system, context, positions, data["box_nm"], indices, vector, [.01, .001, .0002, .00005, .00001, .000002, .0000004])
        print(json.dumps({"cutoff_diagnosis": direction["kind"], "jumps": [{key: value for key, value in row.items() if key != "pairs"} for row in cutoff]}), flush=True)
        for step in (.01, .001, .0002, .00005, .00001, .000002, .0000004):
            energies = []
            for sign in (-1, 1):
                proposed = positions.copy()
                proposed[indices] += sign * step * .1 * vector
                context.setPositions(proposed)
                energies.append([context.getState(getEnergy=True, groups=1 << index).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole) for index in range(len(settings))])
            numeric = (np.array(energies[1]) - np.array(energies[0])) / (2 * step)
            row = {"step_angstrom": step, "analytic": analytic, "numerical": numeric.tolist(), "errors": (numeric - analytic).tolist()}
            rows.append(row)
            print(json.dumps({"direction": direction["kind"], "step": step, "total_error": float(sum(row["errors"])),
                              "groups": {settings[index]["name"]: value for index, value in enumerate(row["errors"]) if abs(value) > .01}}), flush=True)
        context.setPositions(positions)
        report["directions"].append({"kind": direction["kind"], "rows": rows, "cutoff_discontinuities": cutoff})
    write_json(args.chain / "periodic_force_gate_diagnosis.json", report)


if __name__ == "__main__":
    main()
