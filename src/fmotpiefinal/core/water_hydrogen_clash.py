"""Repair only graph-diagnosed water-H collisions; keep all heavy atoms fixed."""
import json
import math
from pathlib import Path
import shutil

import numpy as np
from scipy.spatial import cKDTree

from .component_inventory import sha
from .water_hydrogen_transfer import waters, water_key, valid_water, xyz


def hydrogen_clashes(path, threshold=1.2):
    rows = [line for line in Path(path).read_text().splitlines() if line.startswith(('ATOM  ', 'HETATM'))
            and line[76:78].strip().upper() == 'H']
    if len(rows) < 2:
        return []
    positions = np.array([xyz(line) for line in rows])
    found = []
    for first, second in sorted(cKDTree(positions).query_pairs(threshold)):
        left, right = rows[first], rows[second]
        if water_key(left) == water_key(right) or not any(line[17:20] == 'HOH' for line in (left, right)):
            continue
        selected = sorted({water_key(line) for line in (left, right) if line[17:20] == 'HOH'})
        found.append({'atoms': [f'{line[21].strip() or "_"}:{line[17:20].strip()}:{line[22:27].strip()}/{line[12:16].strip()}'
                               for line in (left, right)],
                      'water_keys': [list(key) for key in selected],
                      'distance_angstrom': float(np.linalg.norm(positions[first] - positions[second])),
                      'threshold_angstrom': threshold})
    return found


def water_geometry_valid(atoms):
    if not valid_water(atoms):
        return False
    oxygen = np.array(xyz(atoms['O']))
    vectors = [np.array(xyz(atoms[name])) - oxygen for name in ('H1', 'H2')]
    cosine = np.dot(*vectors) / np.prod([np.linalg.norm(vector) for vector in vectors])
    return 80 <= math.degrees(math.acos(float(np.clip(cosine, -1, 1)))) <= 130


def repair_water_hydrogen_clashes(prepared, source_water, output):
    prepared, output = Path(prepared), Path(output)
    if output.exists():
        raise FileExistsError('water-H repair requires a fresh output')
    before = hydrogen_clashes(prepared)
    selected = sorted({tuple(key) for event in before for key in event['water_keys']})
    source = waters(Path(source_water)) if source_water is not None else {}
    target = waters(prepared)
    lines = prepared.read_text().splitlines(keepends=True)
    original = list(lines)
    transferred = []
    for key in selected:
        incoming = source.get(key)
        if incoming is None or not water_geometry_valid(incoming):
            continue
        if not water_geometry_valid(target[key]) or math.dist(xyz(incoming['O']), xyz(target[key]['O'])) > 0.01:
            raise ValueError('water-H repair source oxygen mapping or geometry changed')
        for index, line in enumerate(lines):
            if line.startswith('HETATM') and water_key(line) == key and line[12:16].strip() in {'H1', 'H2'}:
                lines[index] = line[:30] + incoming[line[12:16].strip()][30:54] + line[54:]
        transferred.append(list(key))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(''.join(lines))
    remaining = hydrogen_clashes(output)
    mm_used = False
    if remaining:
        import openmm
        from openmm import app, unit
        from .special_chemistry import parameter_paths

        pdb = app.PDBFile(str(output))
        atoms = list(pdb.topology.atoms())
        atom_lines = [index for index, line in enumerate(lines) if line.startswith(('ATOM  ', 'HETATM'))]
        if len(atoms) != len(atom_lines):
            raise ValueError('water-H repair topology/source atom inventory differs')
        mobile = {index for index, line_index in enumerate(atom_lines) if
                  lines[line_index][76:78].strip() == 'H' and water_key(lines[line_index]) in selected}
        forcefield = app.ForceField(*map(str, parameter_paths()))
        forcefield.getMatchingTemplates(pdb.topology)
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff,
                                        constraints=None, rigidWater=False)
        for index in range(len(atoms)):
            if index not in mobile:
                system.setParticleMass(index, 0)
        integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
        context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('CPU'), {'Threads': '1'})
        context.setPositions(pdb.positions)
        energy_before = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        openmm.LocalEnergyMinimizer.minimize(context, 1 * unit.kilojoule_per_mole / unit.nanometer, 500)
        state = context.getState(getPositions=True, getEnergy=True)
        energy_after = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        final = np.asarray(state.getPositions().value_in_unit(unit.angstrom))
        initial = np.asarray(pdb.positions.value_in_unit(unit.angstrom))
        fixed = sorted(set(range(len(atoms))) - mobile)
        if (not np.isfinite(final).all() or not np.isfinite([energy_before, energy_after]).all()
                or energy_after > energy_before + 1e-5 or np.max(np.abs(final[fixed] - initial[fixed])) > 1e-5):
            raise ValueError('water-H MM repair failed finite-energy/fixed-atom contract')
        for index in mobile:
            line_index = atom_lines[index]
            lines[line_index] = lines[line_index][:30] + ''.join(f'{value:8.3f}' for value in final[index]) + lines[line_index][54:]
        output.write_text(''.join(lines))
        mm_used = True
        del context, integrator
    after = hydrogen_clashes(output)
    updated = waters(output)
    if after or any(not water_geometry_valid(updated[key]) for key in selected):
        raise ValueError('unresolved water-H clash or invalid repaired water geometry')
    changed = []
    for index, (first, second) in enumerate(zip(original, lines, strict=True)):
        if first == second:
            continue
        if (first[:30] != second[:30] or first[54:] != second[54:]
                or first[76:78].strip() != 'H' or water_key(first) not in selected):
            raise ValueError('water-H repair changed atom inventory or a fixed atom')
        changed.append(index)
    for suffix in ('.model.json', '.components.json', '.components.sdf'):
        parent = prepared.with_suffix(suffix)
        if not parent.exists():
            continue
        destination = output.with_suffix(suffix)
        if suffix.endswith('.json'):
            model = json.loads(parent.read_text())
            if model.get('prepared_sha256') != sha(prepared):
                raise ValueError('registered component parent receptor changed')
            model.update(prepared_sha256=sha(output), coordinate_repair_parent_sha256=sha(prepared),
                         coordinate_repair='validated-selected-water-H-only')
            destination.write_text(json.dumps(model, indent=2) + '\n')
        else:
            shutil.copyfile(parent, destination)
    report = {'schema': 'graph-selected-water-H-clash-repair-v1', 'source_sha256': sha(prepared),
              'source_water_sha256': sha(Path(source_water)) if source_water is not None else None,
              'output_sha256': sha(output), 'before': before, 'after': after,
              'selected_water_keys': [list(key) for key in selected], 'transferred_source_water_keys': transferred,
              'mm_fallback_used': mm_used, 'n_changed_water_hydrogens': len(changed),
              'all_other_atoms_byte_identical': True, 'experimental_labels_used': False,
              'water_or_protein_atoms_deleted': 0}
    output.with_suffix('.water_H_clash.json').write_text(json.dumps(report, indent=2) + '\n')
    return report
