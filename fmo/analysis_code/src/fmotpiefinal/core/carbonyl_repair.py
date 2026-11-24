"""Repair graph-selected distorted backbone carbonyl O atoms with other atoms fixed."""
import json
from pathlib import Path

import numpy as np

from .component_inventory import sha


def selected_carbonyl_oxygens(tasks):
    selected = set()
    for task in tasks:
        if task['code'] != 'PROTEIN_BOND_LENGTH_SCREEN':
            continue
        edge = task['evidence']
        first_residue, first_name = edge['source'].rsplit('/', 1)
        second_residue, second_name = edge['target'].rsplit('/', 1)
        if first_residue == second_residue and {first_name, second_name} == {'C', 'O'}:
            selected.add(edge['source'] if first_name == 'O' else edge['target'])
    return sorted(selected)


def repair_carbonyl_oxygens(source, tasks, output):
    import openmm
    from openmm import app, unit
    source, output = Path(source), Path(output)
    selected = selected_carbonyl_oxygens(tasks)
    if not selected or output.exists():
        raise ValueError('graph-selected carbonyl oxygens and fresh output required')
    if any(source.with_suffix(suffix).exists() for suffix in ('.model.json', '.components.json')):
        raise ValueError('combined carbonyl repair and special-component forcefield requires joint validation')
    pdb = app.PDBFile(str(source))
    atoms = list(pdb.topology.atoms())
    identifiers = [f"{atom.residue.chain.id.strip() or '_'}:{atom.residue.name}:{atom.residue.id}{atom.residue.insertionCode.strip()}/{atom.name}" for atom in atoms]
    if len(set(identifiers)) != len(identifiers) or not set(selected) <= set(identifiers):
        raise ValueError('ambiguous graph-to-topology carbonyl identity')
    mobile = {identifiers.index(name) for name in selected}
    neighbors = {index: [] for index in mobile}
    for first, second in pdb.topology.bonds():
        for oxygen, partner in ((first, second), (second, first)):
            if oxygen.index in mobile:
                neighbors[oxygen.index].append(partner)
    for index, bonded in neighbors.items():
        if (atoms[index].name != 'O' or atoms[index].element.symbol != 'O' or len(bonded) != 1 or
                bonded[0].name != 'C' or bonded[0].element.symbol != 'C' or bonded[0].residue != atoms[index].residue):
            raise ValueError('carbonyl repair requires unique same-residue backbone C=O topology')
    forcefield = app.ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml')
    forcefield.getMatchingTemplates(pdb.topology)
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False)
    fixed = sorted(set(range(len(atoms))) - mobile)
    if not fixed:
        raise ValueError('carbonyl repair requires fixed surrounding topology')
    for index in fixed:
        system.setParticleMass(index, 0.)
    restraint = openmm.CustomExternalForce('0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
    restraint.addGlobalParameter('k', 1000.)
    for coordinate in ('x0', 'y0', 'z0'):
        restraint.addPerParticleParameter(coordinate)
    initial_nm = np.asarray(pdb.positions.value_in_unit(unit.nanometer))
    for index in sorted(mobile):
        restraint.addParticle(index, initial_nm[index])
    system.addForce(restraint)
    integrator = openmm.VerletIntegrator(.001 * unit.picoseconds)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('CPU'))
    context.setPositions(pdb.positions)
    before = float(context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    openmm.LocalEnergyMinimizer.minimize(context, 1. * unit.kilojoule_per_mole / unit.nanometer, maxIterations=1000)
    state = context.getState(getPositions=True, getEnergy=True)
    after = float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    initial = initial_nm * 10.
    final = np.asarray(state.getPositions().value_in_unit(unit.angstrom))
    if (not np.isfinite(final).all() or not np.isfinite([before, after]).all() or after > before + 1. or
            np.linalg.norm(final[fixed] - initial[fixed], axis=1).max() > 1e-5):
        raise ValueError('carbonyl repair failed energy or fixed-atom contract')
    bonds = []
    for index, bonded in neighbors.items():
        partner = bonded[0].index
        distance = float(np.linalg.norm(final[index] - final[partner]))
        movement = float(np.linalg.norm(final[index] - initial[index]))
        if not 1.15 <= distance <= 1.35 or movement > .75:
            raise ValueError('bounded carbonyl repair could not restore validated geometry')
        bonds.append({'oxygen': identifiers[index], 'carbon': identifiers[partner],
            'before_angstrom': float(np.linalg.norm(initial[index] - initial[partner])),
            'after_angstrom': distance, 'oxygen_displacement_angstrom': movement})
    lines = source.read_text().splitlines()
    indices = [index for index, line in enumerate(lines) if line.startswith(('ATOM  ', 'HETATM'))]
    if len(indices) != len(atoms):
        raise ValueError('carbonyl PDB atom-order coverage mismatch')
    for index in sorted(mobile):
        position = indices[index]
        lines[position] = lines[position][:30] + ''.join(f'{value:8.3f}' for value in final[index]) + lines[position][54:]
    output.write_text('\n'.join(lines) + '\n')
    report = {'schema': 'graph-selected-backbone-carbonyl-O-repair-v1', 'experimental_labels_used': False,
        'source_sha256': sha(source), 'output_sha256': sha(output), 'selected_oxygens': selected,
        'potential_kj_mol': [before, after], 'bonds': bonds, 'all_other_atoms_byte_identical': True,
        'n_atoms_deleted': 0, 'restraint_k_kj_mol_angstrom2': 10., 'maximum_allowed_displacement_angstrom': .75,
        'limitation': 'complete standard topology and anomalous backbone C=O only; not global receptor relaxation or affinity optimization'}
    output.with_suffix('.carbonyl_repair.json').write_text(json.dumps(report, indent=2) + '\n')
    del context, integrator
    return report
