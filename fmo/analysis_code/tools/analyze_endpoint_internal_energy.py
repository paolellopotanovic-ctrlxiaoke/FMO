"""Evaluate cached ligand MM Hamiltonians on fixed, physical endpoint frames."""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import time
import xml.etree.ElementTree as ET

import numpy as np

from analyze_openfe_endpoints import ligand_bonds, read_topology, summarize, unwrap_ligand
from run_matched_openfe import digest, write_json


def physical_atom_maps(atoms, mapping):
    first = [index for index, atom in enumerate(atoms) if atom["class"] in (.25, .5)]
    added = [index for index, atom in enumerate(atoms) if atom["class"] == .75]
    molecules = {0: mapping.componentA.to_rdkit(), 10: mapping.componentB.to_rdkit()}
    if len(first) != molecules[0].GetNumAtoms():
        raise ValueError("state A atom inventory mismatch")
    second = {new: first[old] for old, new in mapping.componentA_to_componentB.items()}
    unmatched = sorted(set(range(molecules[10].GetNumAtoms())) - set(second))
    if len(unmatched) != len(added):
        raise ValueError("state B unique-atom inventory mismatch")
    second.update(zip(unmatched, added, strict=True))
    indices = {0: np.array(first), 10: np.array([second[index] for index in range(len(second))])}
    for endpoint, molecule in molecules.items():
        for local, saved in enumerate(indices[endpoint]):
            if molecule.GetAtomWithIdx(local).GetSymbol().upper() != atoms[saved]["element"].upper():
                raise ValueError("mapped element differs from registered atom")
    return indices, molecules


def verify_bonds(topology_path, indices, molecule):
    actual = set(ligand_bonds(topology_path, indices))
    expected = {tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))) for bond in molecule.GetBonds()}
    if actual != expected:
        raise ValueError("endpoint bond graph differs from registered ligand")
    return sorted(actual)


def parameterized_context(database_path, molecule):
    import openmm
    from openmm import app, unit

    database = json.loads(database_path.read_text())
    templates = [entry["ffxml"] for entry in database["openff-2.2.0"].values()]
    selected = []
    for text in templates:
        root = ET.fromstring(text)
        residue = root.find("./Residues/Residue")
        if len(residue.findall("Atom")) == molecule.GetNumAtoms():
            selected.append((text, root, residue))
    if len(selected) != 1:
        raise ValueError("unique cached ligand force-field template required")
    text, root, residue = selected[0]
    type_elements = {entry.attrib["name"]: entry.attrib["element"] for entry in root.findall("./AtomTypes/Type")}
    topology = app.Topology()
    chain = topology.addChain()
    top_residue = topology.addResidue(residue.attrib["name"], chain)
    top_atoms = []
    charges = []
    for index, atom in enumerate(residue.findall("Atom")):
        symbol = type_elements[atom.attrib["type"]]
        if symbol != molecule.GetAtomWithIdx(index).GetSymbol():
            raise ValueError("cached template atom order differs from registered molecule")
        top_atoms.append(topology.addAtom(atom.attrib["name"], app.Element.getBySymbol(symbol), top_residue))
        charges.append(float(atom.attrib["charge"]))
    expected_charges = [atom.GetDoubleProp("PartialCharge") for atom in molecule.GetAtoms()]
    if not np.allclose(charges, expected_charges, atol=1e-10):
        raise ValueError("cached charges differ from frozen component")
    for bond in molecule.GetBonds():
        topology.addBond(top_atoms[bond.GetBeginAtomIdx()], top_atoms[bond.GetEndAtomIdx()])
    forcefield = app.ForceField(io.StringIO(text))
    system = forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff, constraints=None,
                                    rigidWater=False, removeCMMotion=False)
    if system.getNumConstraints() != 0 or system.getNumParticles() != molecule.GetNumAtoms():
        raise ValueError("unexpected isolated ligand system")
    groups = {}
    for index, force in enumerate(system.getForces()):
        label = type(force).__name__
        if label not in ("HarmonicBondForce", "HarmonicAngleForce", "PeriodicTorsionForce", "NonbondedForce"):
            raise ValueError(f"unsupported force group: {label}")
        force.setForceGroup(index)
        groups[index] = label
    nonbonded = next(force for force in system.getForces() if isinstance(force, openmm.NonbondedForce))
    lennard_jones = openmm.XmlSerializer.deserialize(openmm.XmlSerializer.serialize(nonbonded))
    for index in range(lennard_jones.getNumParticles()):
        _, sigma, epsilon = lennard_jones.getParticleParameters(index)
        lennard_jones.setParticleParameters(index, 0, sigma, epsilon)
    for index in range(lennard_jones.getNumExceptions()):
        first, second, _, sigma, epsilon = lennard_jones.getExceptionParameters(index)
        lennard_jones.setExceptionParameters(index, first, second, 0, sigma, epsilon)
    lennard_jones.setForceGroup(31)
    system.addForce(lennard_jones)
    integrator = openmm.VerletIntegrator(.001 * unit.picosecond)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName("Reference"))
    metadata = {"cached_database_sha256": digest(database_path), "constraints": 0,
                "system_xml_sha256": hashlib.sha256(openmm.XmlSerializer.serialize(system).encode()).hexdigest(),
                "charges_from_frozen_cache": True, "platform": "Reference",
                "diagnostic_force_group_31": "LJ-only copy excluded from original potential energy"}
    return context, integrator, groups, metadata


def evaluate(context, groups, coordinates):
    from openmm import unit

    context.setPositions(coordinates * .1 * unit.nanometer)
    components = {name: float(context.getState(getEnergy=True, groups=1 << group).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole))
                  for group, name in groups.items()}
    original_mask = sum(1 << group for group in groups)
    total = float(context.getState(getEnergy=True, groups=original_mask).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole))
    if not np.isfinite(total) or abs(sum(components.values()) - total) > 1e-6:
        raise ValueError("nonfinite energy or force-group closure failure")
    lennard_jones = float(context.getState(getEnergy=True, groups=1 << 31).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole))
    return {"total": total, **components, "Nonbonded_LennardJones": lennard_jones,
            "Nonbonded_Coulomb": components["NonbondedForce"] - lennard_jones}


def mean_difference(bound, solution):
    keys = bound[0]["mean_kcal_mol"].keys()
    return {key: float(np.mean([row["mean_kcal_mol"][key] for row in bound]) -
                       np.mean([row["mean_kcal_mol"][key] for row in solution])) for key in keys}


def frame_coordinates(data, frozen_frame, endpoint, indices, bonds):
    iteration, replica = frozen_frame["iteration"], frozen_frame["replica"]
    if int(data["states"][iteration, replica]) != endpoint:
        raise ValueError("frozen frame endpoint mismatch")
    raw = data["positions"][iteration, replica]
    if np.ma.getmaskarray(raw).any():
        raise ValueError("missing saved coordinates")
    box = np.asarray(data["box_vectors"][iteration, replica]) * 10
    return unwrap_ligand(np.asarray(raw)[indices] * 10, bonds, box)


def analyze_target(run, target, diagnostics):
    from gufe import LigandAtomMapping
    from netCDF4 import Dataset

    result = {"endpoints": {}, "legs": {}}
    for leg in ("solution", "bound"):
        root = run / "production" / target / "edge_000" / leg
        mapping = LigandAtomMapping.from_json(root / "mapping.json")
        topology_path = sorted((root / "shared").glob("*/hybrid_system.pdb"))[0]
        atoms = read_topology(topology_path)
        indices, molecules = physical_atom_maps(atoms, mapping)
        bonds = {endpoint: verify_bonds(topology_path, selected, molecules[endpoint]) for endpoint, selected in indices.items()}
        result["legs"][leg] = {"repeats": [], "mapping_sha256": digest(root / "mapping.json"),
                                "topology_sha256": digest(topology_path), "systems": {}}
        for endpoint, component in ((0, mapping.componentA), (10, mapping.componentB)):
            state = "A" if endpoint == 0 else "B"
            database_path = topology_path.with_name(f"{state}_db.json")
            context, integrator, groups, metadata = parameterized_context(database_path, molecules[endpoint])
            result["legs"][leg]["systems"][str(endpoint)] = metadata
            try:
                for repeat in diagnostics[leg]["complete_repeats"]:
                    if repeat["endpoint"] != endpoint:
                        continue
                    trajectory = Path(repeat["trajectory"])
                    frames = []
                    with Dataset(trajectory) as data:
                        for frozen_frame in repeat["frames"]:
                            iteration, replica = frozen_frame["iteration"], frozen_frame["replica"]
                            coordinates = frame_coordinates(data, frozen_frame, endpoint, indices[endpoint], bonds[endpoint])
                            frames.append({"iteration": iteration, "replica": replica,
                                           "energies_kcal_mol": evaluate(context, groups, coordinates)})
                    summary = {key: summarize([frame["energies_kcal_mol"][key] for frame in frames])
                               for key in frames[0]["energies_kcal_mol"]}
                    result["legs"][leg]["repeats"].append({"endpoint": endpoint, "name": component.name,
                        "trajectory": str(trajectory), "n_frames": len(frames), "summary_kcal_mol": summary,
                        "mean_kcal_mol": {key: value["mean"] for key, value in summary.items()}, "frames": frames})
            finally:
                del context, integrator
    for endpoint in (0, 10):
        if result["legs"]["bound"]["systems"][str(endpoint)]["system_xml_sha256"] != result["legs"]["solution"]["systems"][str(endpoint)]["system_xml_sha256"]:
            raise ValueError("bound/solution isolated ligand Hamiltonians differ")
        legs = {leg: [row for row in result["legs"][leg]["repeats"] if row["endpoint"] == endpoint]
                for leg in ("solution", "bound")}
        if any(not rows for rows in legs.values()):
            raise ValueError("completed bound and solution repeats required")
        result["endpoints"][str(endpoint)] = {"name": diagnostics["endpoint_names"][str(endpoint)],
            "bound_minus_solution_mean_internal_energy_kcal_mol": mean_difference(legs["bound"], legs["solution"]),
            "n_repeats": {leg: len(rows) for leg, rows in legs.items()}}
    first = result["endpoints"]["0"]["bound_minus_solution_mean_internal_energy_kcal_mol"]
    second = result["endpoints"]["10"]["bound_minus_solution_mean_internal_energy_kcal_mol"]
    result["B_minus_A_internal_energy_double_difference_kcal_mol"] = {key: second[key] - first[key] for key in first}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--diagnostics", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    diagnostics = json.loads(args.diagnostics.read_text())
    if diagnostics["affinity_used_for_frame_selection"] or diagnostics["frozen_protocol_sha256"] != digest(args.run / "frozen_protocol.json"):
        raise ValueError("matching frozen, label-independent trajectory diagnostic required")
    report = {"schema": "openfe-endpoint-internal-energy-v1", "targets": {},
              "diagnostics_sha256": digest(args.diagnostics), "analysis_driver_sha256": digest(__file__),
              "selection": diagnostics["selection"], "labels_used": False, "main_pipeline_changed": False,
              "definition": "Mean isolated-ligand vacuum MM potential on bound frames minus mean on solution frames; double difference is B minus A.",
              "limitations": ["Not a binding free energy, entropy or a QM/MM reweighting correction.",
                  "Nonbonded_Coulomb and Nonbonded_LennardJones partition NonbondedForce; do not add all three.",
                  "No solvent energy or ligand-protein interaction in this isolated potential.",
                  "Bonds are re-evaluated without constraints; constrained-trajectory geometry may bias bond terms.",
                  "Frames are correlated; raw frame variance is not independent-repeat uncertainty.",
                  "Within-molecule differences are meaningful diagnostics; absolute energies across different ligands are not affinity scores."]}
    for target in ("BACE", "CDK8"):
        report["targets"][target] = analyze_target(args.run, target, diagnostics["targets"][target])
        write_json(args.output, report)
        print(target, json.dumps(report["targets"][target]["B_minus_A_internal_energy_double_difference_kcal_mol"]), flush=True)
    report["wall_seconds"] = time.monotonic() - started
    report["status"] = "COMPLETE"
    write_json(args.output, report)


if __name__ == "__main__":
    main()
