"""Topology-validated restrained local amide-motif repair; remote atoms fixed."""
import hashlib
import json
from pathlib import Path

import numpy as np


def identify_amide_nitrogens(tasks):
    selected = set()
    for task in tasks:
        if task["code"] not in {"INVALID_COVALENT_LINK_GEOMETRY", "PROTEIN_BOND_LENGTH_SCREEN"}:
            continue
        edge = task["evidence"]
        first, second = edge["source"], edge["target"]
        first_residue, first_name = first.rsplit("/", 1)
        second_residue, second_name = second.rsplit("/", 1)
        if first_residue == second_residue and {first_name, second_name} == {"N", "CA"}:
            selected.add(first if first_name == "N" else second)
        elif first_residue != second_residue and {first_name, second_name} == {"C", "N"}:
            selected.add(first if first_name == "N" else second)
    return sorted(selected)


def repair_amide_nitrogens(source, tasks, output):
    import openmm
    from openmm import app, unit
    source, output = Path(source), Path(output)
    selected = identify_amide_nitrogens(tasks)
    if not selected:
        raise ValueError("no graph-selected anomalous amide N")
    pdb = app.PDBFile(str(source))
    atoms = list(pdb.topology.atoms())
    identifiers = [f"{atom.residue.chain.id.strip() or '_'}:{atom.residue.name}:{atom.residue.id}{atom.residue.insertionCode.strip()}/{atom.name}" for atom in atoms]
    if len(identifiers) != len(set(identifiers)) or not set(selected) <= set(identifiers):
        raise ValueError("ambiguous or missing graph-to-topology atom identity")
    mobile = {identifiers.index(name) for name in selected}
    neighbors = {index: [] for index in mobile}
    for first, second in pdb.topology.bonds():
        for nitrogen, partner in ((first, second), (second, first)):
            if nitrogen.index in neighbors:
                neighbors[nitrogen.index].append(partner)
    mobile_residues = set()
    for index, bonded in neighbors.items():
        heavy = [atom for atom in bonded if atom.element.symbol != "H"]
        if sorted(atom.name for atom in heavy) != ["C", "CA"]:
            raise ValueError("amide N repair requires validated internal peptide C-N-CA connectivity")
        mobile_residues.add(atoms[index].residue)
        mobile_residues.update(atom.residue for atom in heavy)
    mobile = {atom.index for residue in mobile_residues for atom in residue.atoms()}
    forcefield = app.ForceField("amber14/protein.ff14SB.xml", "amber14/tip3p.xml")
    forcefield.getMatchingTemplates(pdb.topology)
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False)
    fixed = sorted(set(range(len(atoms))) - mobile)
    for index in fixed:
        system.setParticleMass(index, 0.)
    restraint = openmm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint.addGlobalParameter("k", 1000.)
    for label in ("x0", "y0", "z0"):
        restraint.addPerParticleParameter(label)
    original_positions = np.asarray(pdb.positions.value_in_unit(unit.nanometer))
    for index in sorted(mobile):
        if atoms[index].element.symbol != "H":
            restraint.addParticle(index, original_positions[index])
    system.addForce(restraint)
    integrator = openmm.VerletIntegrator(.001 * unit.picoseconds)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName("CPU"))
    context.setPositions(pdb.positions)
    before = float(context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    openmm.LocalEnergyMinimizer.minimize(context, 1. * unit.kilojoule_per_mole / unit.nanometer, maxIterations=0)
    state = context.getState(getPositions=True, getEnergy=True)
    after = float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    initial = np.asarray(pdb.positions.value_in_unit(unit.angstrom))
    final = np.asarray(state.getPositions().value_in_unit(unit.angstrom))
    if not np.isfinite(final).all() or not np.isfinite([before, after]).all() or after > before + 1. or np.linalg.norm(final[fixed] - initial[fixed], axis=1).max() > 1e-5:
        raise ValueError("local amide repair failed energy/fixed-atom contract")
    lengths = []
    for index, bonded in neighbors.items():
        for partner in bonded:
            if partner.element.symbol == "H":
                continue
            distance = float(np.linalg.norm(final[index] - final[partner.index]))
            lower, upper = (1.15, 1.5) if partner.name == "C" else (1.3, 1.65)
            if not lower <= distance <= upper:
                raise ValueError(f"local restrained-motif repair could not restore amide bond geometry: {distance}")
            lengths.append({"source": identifiers[index], "target": identifiers[partner.index],
                "before_angstrom": float(np.linalg.norm(initial[index] - initial[partner.index])), "after_angstrom": distance})
    lines = source.read_text().splitlines()
    atom_indices = [index for index, line in enumerate(lines) if line.startswith(("ATOM  ", "HETATM"))]
    if len(atom_indices) != len(atoms):
        raise ValueError("source/topology atom inventory changed")
    for index in mobile:
        line_index = atom_indices[index]
        lines[line_index] = lines[line_index][:30] + "".join(f"{value:8.3f}" for value in final[index]) + lines[line_index][54:]
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError("local repair must use a fresh receptor path")
    output.write_text("\n".join(lines) + "\n")
    report = {"schema": "graph-selected-restrained-amide-motif-repair-v1", "experimental_labels_used": False,
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(), "output_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "n_atoms": len(atoms), "selected_nitrogens": selected, "mobile_atoms": [identifiers[index] for index in sorted(mobile)],
        "potential_kj_mol": [before, after], "bonds": lengths, "all_other_atoms_fixed": True,
        "mobile_heavy_restraint_k_kj_mol_angstrom2": 10.,
        "limitation": "local repair only; not missing-loop reconstruction, covalent chemistry reassignment or general receptor relaxation"}
    output.with_suffix(".repair.json").write_text(json.dumps(report, indent=2) + "\n")
    del context, integrator
    return report
