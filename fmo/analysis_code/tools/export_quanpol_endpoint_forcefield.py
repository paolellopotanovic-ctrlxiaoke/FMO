"""Export fixed physical OpenFE environment and cached ligand parameters to QuanPol."""
from __future__ import annotations

import argparse
import bz2
import json
from pathlib import Path

import numpy as np

from analyze_endpoint_internal_energy import parameterized_context
from quanpol_constants import DEFAULT_SWRB_ANGSTROM
from run_matched_openfe import digest, write_json


def section(name, rows):
    return [name, *rows, "STOP"] if rows else []


def export(task):
    import openmm
    from gufe import LigandAtomMapping
    from openmm import unit

    directory = Path(task["deck"]).parent
    data = np.load(directory / "physical_coordinates.npz")
    topology = Path(task["topology"])
    source = topology.with_name("hybrid_system.xml.bz2")
    original = openmm.XmlSerializer.deserialize(bz2.decompress(source.read_bytes()).decode())
    mapping = LigandAtomMapping.from_json(topology.parents[2] / "mapping.json")
    component = mapping.componentA if task["endpoint"] == 0 else mapping.componentB
    molecule = component.to_rdkit()
    state = "A" if task["endpoint"] == 0 else "B"
    context, integrator, groups, metadata = parameterized_context(topology.with_name(f"{state}_db.json"), molecule)
    ligand_system = openmm.XmlSerializer.deserialize(openmm.XmlSerializer.serialize(context.getSystem()))
    del context, integrator
    ligand_count = len(data["ligand_xyz"])
    environment_indices = data["environment_original_indices"].astype(int).tolist()
    environment_map = {index: ligand_count + offset + 1 for offset, index in enumerate(environment_indices)}
    ligand_map = {index: index + 1 for index in range(ligand_count)}
    coordinates = np.vstack((data["ligand_xyz"], data["environment_xyz"]))
    parameters = []
    bond_rows, angle_rows, torsion_rows = [], [], []
    existing_bonds = set()
    export_counts = []
    for system, selected, label in ((ligand_system, ligand_map, "cached_ligand"),
                                    (original, environment_map, "physical_environment")):
        nonbonded = next(force for force in system.getForces()
                         if isinstance(force, openmm.NonbondedForce) and force.getForceGroup() != 31)
        if label == "physical_environment":
            for offset in range(nonbonded.getNumParticleParameterOffsets()):
                _, index, charge, sigma, epsilon = nonbonded.getParticleParameterOffset(offset)
                if index in selected and any(abs(value) > 1e-12 for value in (charge, sigma, epsilon)):
                    raise ValueError("alchemical parameter offset in environment")
        for index in selected:
            charge, sigma, epsilon = nonbonded.getParticleParameters(index)
            parameters.append([system.getParticleMass(index).value_in_unit(unit.dalton),
                               charge.value_in_unit(unit.elementary_charge),
                               sigma.value_in_unit(unit.angstrom),
                               epsilon.value_in_unit(unit.kilocalorie_per_mole)])
        counts = {"system": label, "custom_terms_with_retained_environment_atoms": {}}
        for force in system.getForces():
            if force.getForceGroup() == 31:
                continue
            if isinstance(force, openmm.HarmonicBondForce):
                for offset in range(force.getNumBonds()):
                    first, second, length, constant = force.getBondParameters(offset)
                    if first not in selected and second not in selected:
                        continue
                    if first not in selected or second not in selected:
                        raise ValueError("export would cut a covalent bond")
                    pair = [selected[first], selected[second]]
                    existing_bonds.add(tuple(sorted(pair)))
                    coefficient = .5 * constant.value_in_unit(unit.kilocalorie_per_mole / unit.angstrom**2)
                    bond_rows.append(f"{len(bond_rows)+1} {' '.join(map(str, pair))} {coefficient:.12g} {length.value_in_unit(unit.angstrom):.12g}")
            elif isinstance(force, openmm.HarmonicAngleForce):
                for offset in range(force.getNumAngles()):
                    first, middle, last, angle, constant = force.getAngleParameters(offset)
                    indices = (first, middle, last)
                    if not any(index in selected for index in indices):
                        continue
                    if not all(index in selected for index in indices):
                        raise ValueError("export would cut an angle")
                    coefficient = .5 * constant.value_in_unit(unit.kilocalorie_per_mole / unit.radian**2)
                    angle_rows.append(f"{len(angle_rows)+1} {' '.join(str(selected[index]) for index in indices)} {coefficient:.12g} {angle.value_in_unit(unit.degree):.12g}")
            elif isinstance(force, openmm.PeriodicTorsionForce):
                for offset in range(force.getNumTorsions()):
                    first, second, third, fourth, periodicity, phase, constant = force.getTorsionParameters(offset)
                    indices = (first, second, third, fourth)
                    if not any(index in selected for index in indices):
                        continue
                    if not all(index in selected for index in indices):
                        raise ValueError("export would cut a torsion")
                    torsion_rows.append(f"{len(torsion_rows)+1} {' '.join(str(selected[index]) for index in indices)} {constant.value_in_unit(unit.kilocalorie_per_mole):.12g} {periodicity} {phase.value_in_unit(unit.degree):.12g}")
            elif label == "physical_environment" and isinstance(force, (openmm.CustomBondForce, openmm.CustomAngleForce, openmm.CustomTorsionForce)):
                kind = type(force).__name__
                if kind == "CustomBondForce":
                    count_method, read_method, width = force.getNumBonds, force.getBondParameters, 2
                elif kind == "CustomAngleForce":
                    count_method, read_method, width = force.getNumAngles, force.getAngleParameters, 3
                elif kind == "CustomTorsionForce":
                    count_method, read_method, width = force.getNumTorsions, force.getTorsionParameters, 4
                retained = sum(any(index in selected for index in read_method(offset)[:width])
                               for offset in range(count_method()))
                counts["custom_terms_with_retained_environment_atoms"][force.getName()] = retained
                if retained:
                    raise ValueError(f"unsupported retained custom bonded terms: {force.getName()} ({retained})")
        for offset in range(system.getNumConstraints()):
            first, second, length = system.getConstraintParameters(offset)
            if first not in selected and second not in selected:
                continue
            if first not in selected or second not in selected:
                raise ValueError("export would cut a constraint")
            pair = [selected[first], selected[second]]
            if tuple(sorted(pair)) not in existing_bonds:
                existing_bonds.add(tuple(sorted(pair)))
                bond_rows.append(f"{len(bond_rows)+1} {' '.join(map(str, pair))} 0.0 {length.value_in_unit(unit.angstrom):.12g}")
        export_counts.append(counts)
    from openmm import app
    pdb_lines = (directory / "environment.pdb").read_text().splitlines()
    elements = [atom.GetSymbol() for atom in molecule.GetAtoms()]
    elements += [line[76:78].strip().title() for line in pdb_lines if line.startswith("ATOM  ")]
    atomic_numbers = [app.Element.getBySymbol(symbol).atomic_number for symbol in elements]
    coord_rows = [f"{symbol} {number}.0 " + " ".join(f"{value:.9f}" for value in position)
                  for symbol, number, position in zip(elements, atomic_numbers, coordinates, strict=True)]
    parameter_rows = [f"{symbol} {mass:.9g} {charge:.12g} 0 {sigma:.12g} {epsilon:.12g} 0 0"
                      for symbol, (mass, charge, sigma, epsilon) in zip(elements, parameters, strict=True)]
    forcefield = [" $FFDATA", *section("COORDINATES", coord_rows), *section("PARAMETERS", parameter_rows),
                  *section("BOND", bond_rows), *section("ANGLE", angle_rows), *section("DIHROT", torsion_rows), " $END", ""]
    original_deck = Path(task["deck"]).read_text()
    before, remaining = original_deck.split(" $QUANPO", 1)
    _, after = remaining.split(" $END", 1)
    data_group = after.split(" $FFPDB", 1)[0]
    capacity = max(3400, len(coordinates) + 200)
    controls = (f" $QUANPO NSTEP=0 NFFTYP=0 IDOPOL=0 DFTBMM=0.70 LJQMMM=1\n"
                f" LJSIGMA=1 WT14CH=0.8333333333333333 WT14LJ=0.5\n"
                f" MXLIST1={capacity} MXLIST2={capacity} BUFWID1=1.0 BUFWID2=4.0\n"
                f" ISHIFT=0 ISWITCH=0 SWRA=1000 SWRB={DEFAULT_SWRB_ANGSTROM:g} $END\n")
    deck = directory / f"matched_{Path(task['deck']).name}"
    deck.write_text(before + controls + data_group + "\n" + "\n".join(forcefield))
    write_json(directory / "parameter_export.json", {"status": "EXPORTED_NOT_VALIDATED", "deck": str(deck),
        "deck_sha256": digest(deck), "source_system_sha256": digest(source), "ligand": metadata,
        "atom_count": len(coordinates), "neighbor_list_capacity": capacity,
        "bond_count": len(bond_rows), "angle_count": len(angle_rows),
        "torsion_count": len(torsion_rows), "export_counts": export_counts,
        "polarization_model": "fixed-charge, nonpolarizable OpenMM/AMBER environment; FFDATA POL=0 for every exported site",
        "nonzero_polarizable_sites": 0,
        "limitations": ["Finite shell without periodic electrostatics, not the OpenFE PME Hamiltonian.",
                        "This OpenFE AMBER environment has no atomic polarizability parameters; IDOPOL=100 would silently reset to 0 in GAMESS unless nonzero POL values are supplied.",
                        "No-cutoff diagnostic: FFDATA parameter and MM energy parity require validation.",
                        "QM/MM replaces ligand MM terms; damping is part of a different Hamiltonian."]})
    print(task["leg"], task["endpoint"], "EXPORTED", len(coordinates), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    args = parser.parse_args()
    for task in json.loads(args.plan.read_text())["tasks"]:
        export(task)


if __name__ == "__main__":
    main()
