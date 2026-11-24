"""Trace isolated-ligand Coulomb changes to atom pairs on frozen frames."""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np

from analyze_endpoint_internal_energy import parameterized_context
from run_matched_openfe import digest, write_json


def pair_parameters(force):
    from openmm import unit

    charges = np.array([force.getParticleParameters(index)[0].value_in_unit(unit.elementary_charge)
                        for index in range(force.getNumParticles())])
    first, second = np.triu_indices(len(charges), 1)
    products = charges[first] * charges[second]
    locations = {(int(left), int(right)): index for index, (left, right) in enumerate(zip(first, second, strict=True))}
    for index in range(force.getNumExceptions()):
        left, right, product, _, _ = force.getExceptionParameters(index)
        products[locations[tuple(sorted((left, right)))]] = product.value_in_unit(unit.elementary_charge ** 2)
    return first, second, products


def main():
    from gufe import LigandAtomMapping
    import openmm
    from openmmtools.constants import ONE_4PI_EPS0

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--frozen", type=Path, required=True)
    parser.add_argument("--mm", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.frozen.read_text())
    mm = json.loads(args.mm.read_text())
    if plan["diagnostics_sha256"] != mm["diagnostics_sha256"] or len(plan["fixed_iterations"]) != 41:
        raise ValueError("matching complete fixed 41-frame plan required")
    output = {"schema": "endpoint-isolated-coulomb-pair-diagnostic-v1", "targets": {},
              "sources": {str(path): digest(path) for path in (args.frozen, args.mm)},
              "driver_sha256": digest(__file__), "labels_used": False,
              "warning": "Intramolecular isolated-vacuum Coulomb energy, not protein-ligand tPIE or binding free energy; no score changed."}
    for target in ("BACE", "CDK8"):
        root = args.run / "production" / target / "edge_000/solution"
        mapping = LigandAtomMapping.from_json(root / "mapping.json")
        output["targets"][target] = {}
        for endpoint, component in ((0, mapping.componentA), (10, mapping.componentB)):
            state = "A" if endpoint == 0 else "B"
            database = sorted((root / "shared").glob(f"*/{state}_db.json"))[0]
            molecule = component.to_rdkit()
            context, integrator, _, _ = parameterized_context(database, molecule)
            try:
                force = next(force for force in context.getSystem().getForces() if isinstance(force, openmm.NonbondedForce))
                first, second, products = pair_parameters(force)
            finally:
                del context, integrator
            energies = defaultdict(list)
            distances = defaultdict(list)
            closure_errors = []
            for task in plan["tasks"]:
                if task["target"] != target or task["endpoint"] != endpoint:
                    continue
                positions_nm = (np.asarray(task["xyz_angstrom"], dtype=np.float32) * .1).astype(float)
                distance_nm = np.linalg.norm(positions_nm[first] - positions_nm[second], axis=1)
                distance = distance_nm * 10
                energy = ONE_4PI_EPS0 * products / distance_nm / 4.184
                repeats = mm["targets"][target]["legs"][task["leg"]]["repeats"]
                repeat = next(row for row in repeats if row["endpoint"] == endpoint and row["trajectory"] == task["trajectory"])
                frame = next(row for row in repeat["frames"] if row["iteration"] == task["iteration"])
                error = abs(float(energy.sum()) - frame["energies_kcal_mol"]["Nonbonded_Coulomb"])
                closure_errors.append(error)
                if error > 1e-5:
                    raise ValueError(f'pairwise Coulomb sum mismatch: {target} {endpoint} {task["leg"]} {task["iteration"]}; pair={energy.sum()}, original={frame["energies_kcal_mol"]["Nonbonded_Coulomb"]}, constant={ONE_4PI_EPS0}')
                energies[task["leg"]].append(energy)
                distances[task["leg"]].append(distance)
            averages = {leg: np.mean(values, axis=0) for leg, values in energies.items()}
            change = averages["bound"] - averages["solution"]
            order = np.argsort(-np.abs(change))[:15]
            pairs = []
            for index in order:
                atoms = [molecule.GetAtomWithIdx(int(atom)) for atom in (first[index], second[index])]
                pairs.append({"atom_indices_zero_based": [int(first[index]), int(second[index])],
                    "elements": [atom.GetSymbol() for atom in atoms],
                    "charges_e": [atom.GetDoubleProp("PartialCharge") for atom in atoms],
                    "bound_minus_solution_coulomb_kcal_mol": float(change[index]),
                    "mean_distance_angstrom": {leg: float(np.mean(np.asarray(values)[:, index])) for leg, values in distances.items()}})
            output["targets"][target][str(endpoint)] = {"name": component.name,
                "total_bound_minus_solution_coulomb_kcal_mol": float(change.sum()),
                "maximum_pair_sum_closure_error_kcal_mol": max(closure_errors),
                "top_15_absolute_internal_pair_changes": pairs,
                "cancellation_warning": "Large individual pair changes can cancel; do not select only favorable pairs."}
    write_json(args.output, output)


if __name__ == "__main__":
    main()
