"""Open-source restrained minimization replacing unavailable MacroModel/OPLS3e."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import hashlib
import json
import math
import tempfile

import numpy as np
import openmm
import openmm.app as app
import openmm.unit as unit
from openmmforcefields.generators import SystemGenerator
from openff.toolkit import Molecule, Topology
from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
from rdkit import Chem

from .bound_geometry import load_bound_ligand


@dataclass(frozen=True)
class MinimizationResult:
    complex_pdb: Path
    provenance_json: Path
    initial_energy_kj_mol: float
    final_energy_kj_mol: float
    ligand_rmsd_angstrom: float
    ligand_max_displacement_angstrom: float


def _selected_rdkit_ligand(ligand_sdf: str | Path, ligand_name: str):
    supplier = Chem.SDMolSupplier(str(ligand_sdf), removeHs=False, sanitize=True)
    matches = [molecule for molecule in supplier
               if molecule is not None and molecule.GetProp("_Name").strip() == ligand_name]
    if len(matches) != 1:
        raise ValueError(f"expected exactly one ligand {ligand_name!r}, found {len(matches)}")
    molecule = matches[0]
    Chem.AssignStereochemistryFrom3D(molecule)
    return molecule


def _water_compatible_pdb(pdb_path: str | Path) -> Path:
    source = Path(pdb_path)
    text = source.read_text(encoding="utf-8")
    lines = []
    water_hydrogen_index: dict[tuple[str, str], int] = {}
    for line in text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and line[17:20] == "WAT":
            line = line[:17] + "HOH" + line[20:]
            water_key = (line[21], line[22:26])
            atom_name = line[12:16].strip().upper()
            if atom_name in {"OW", "O", "OX"}:
                replacement = " O  "
            elif atom_name.startswith("H"):
                hydrogen_index = water_hydrogen_index.get(water_key, 0)
                replacement = " H1 " if hydrogen_index == 0 else " H2 "
                water_hydrogen_index[water_key] = hydrogen_index + 1
            else:
                raise ValueError(f"unsupported water atom name {atom_name!r}")
            line = line[:12] + replacement + line[16:]
        lines.append(line)
    identity = hashlib.sha256(str(source.resolve()).encode("utf-8")).hexdigest()[:16]
    output = Path(tempfile.gettempdir()) / f"fmotpiefinal-{source.stem}-{identity}-omm.pdb"
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=output.parent,
                                     prefix=output.name + ".", delete=False) as handle:
        temporary = Path(handle.name)
        handle.write("\n".join(lines) + "\n")
    try:
        temporary.replace(output)
    finally:
        temporary.unlink(missing_ok=True)
    return output


def _restraint_atom_indices(topology_atoms, ligand_atom_indices, restrain_hydrogens: bool):
    return [atom_index for atom_index in ligand_atom_indices
            if restrain_hydrogens or topology_atoms[atom_index].element.symbol.upper() != "H"]


def _local_sidechain_indices(topology_atoms, protein_atom_count, residues):
    backbone = {"N", "CA", "C", "O", "OXT", "H", "H1", "H2", "H3", "HA", "HA2", "HA3"}
    selected = set(residues)
    return [index for index, atom in enumerate(topology_atoms[:protein_atom_count])
            if (atom.residue.chain.id, atom.residue.id, atom.residue.name) in selected
            and atom.residue.name not in {"HOH", "WAT"} and atom.name not in backbone]


def minimize_complex(protein_pdb: str | Path, ligand_sdf: str | Path, ligand_name: str,
                     output_dir: str | Path, restraint_k_kj_mol_a2: float = 100.0,
                     tolerance_kj_mol_nm: float = 1.0, *, move_protein: bool = True,
                     fix_ligand_heavy: bool = False,
                     relax_water_hydrogens: bool = False,
                     fix_ligand_all: bool = False,
                     restrain_ligand_hydrogens: bool = True,
                     starting_ligand_complex: str | Path | None = None,
                     local_flexible_residues: tuple[tuple[str, str, str], ...] = (),
                     receptor_restraint_k_kj_mol_a2: float = 2.5) -> MinimizationResult:
    """Minimize a prepared complex while restraining every ligand atom.

    This low-level helper defaults to movable receptor and waters for the
    whole-complex SOP. The unified blind workflow explicitly passes
    ``move_protein=False`` for its standard fixed-receptor k10 arm and
    ``move_protein=True`` only for a registered clash/SCF rescue arm.
    """
    if restraint_k_kj_mol_a2 <= 0:
        raise ValueError("restraint force constant must be positive")
    if local_flexible_residues and (move_protein or relax_water_hydrogens or fix_ligand_heavy or fix_ligand_all):
        raise ValueError("local sidechain flexibility requires fixed backbone/waters and a mobile ligand")
    if local_flexible_residues and (not math.isfinite(receptor_restraint_k_kj_mol_a2) or receptor_restraint_k_kj_mol_a2 <= 0):
        raise ValueError("local receptor restraint must be finite and positive")
    if move_protein and fix_ligand_heavy:
        raise ValueError("fixed ligand heavy atoms require a fixed prepared receptor")
    if relax_water_hydrogens and move_protein:
        raise ValueError("water-hydrogen-only receptor relaxation requires fixed protein and water oxygens")
    if fix_ligand_all and (not relax_water_hydrogens or fix_ligand_heavy):
        raise ValueError("fully fixed ligand requires water-hydrogen-only relaxation, not fixed-heavy mode")
    rdkit_ligand = _selected_rdkit_ligand(ligand_sdf, ligand_name)
    if starting_ligand_complex is not None:
        ligand_pose = load_bound_ligand(ligand_sdf, ligand_name, starting_ligand_complex)
        for atom_index, atom in enumerate(ligand_pose.atoms):
            rdkit_ligand.GetConformer().SetAtomPosition(atom_index, atom.xyz)
    off_ligand = Molecule.from_rdkit(rdkit_ligand, allow_undefined_stereo=False)
    off_ligand.assign_partial_charges(
        partial_charge_method="openff-gnn-am1bcc-1.0.0.pt",
        toolkit_registry=NAGLToolkitWrapper(),
    )
    from .special_chemistry import registered_phosphoresidue_model, load_phosphoresidue_bonds, parameter_paths
    special_model = registered_phosphoresidue_model(protein_pdb)
    if special_model is not None:
        load_phosphoresidue_bonds()
    protein_pdb_path = _water_compatible_pdb(protein_pdb)
    protein = app.PDBFile(str(protein_pdb_path))
    ligand_topology = Topology.from_molecules([off_ligand]).to_openmm()
    ligand_coordinates = np.asarray(rdkit_ligand.GetConformer().GetPositions(), dtype=float) / 10.0
    modeller = app.Modeller(protein.topology, protein.positions)
    modeller.add(ligand_topology, ligand_coordinates * unit.nanometer)

    forcefield_files = ["amber14/protein.ff14SB.xml", "amber14/tip3p.xml"]
    if special_model is not None:
        forcefield_files = list(map(str, parameter_paths()))
    from .tartrate import component_molecules
    receptor_molecules, receptor_models = component_molecules(protein_pdb)
    generator = SystemGenerator(
        forcefields=forcefield_files,
        small_molecule_forcefield="openff-2.2.0",
        forcefield_kwargs={"constraints": None, "rigidWater": False},
        nonperiodic_forcefield_kwargs={"nonbondedMethod": app.NoCutoff},
        molecules=[off_ligand, *receptor_molecules],
        cache=None,
    )
    system = generator.create_system(modeller.topology)
    topology_atoms = list(modeller.topology.atoms())
    protein_atom_count = len(list(protein.topology.atoms()))
    ligand_atom_indices = list(range(protein_atom_count, len(topology_atoms)))
    local_mobile = _local_sidechain_indices(topology_atoms, protein_atom_count, local_flexible_residues)
    if local_flexible_residues:
        represented = {(atom.residue.chain.id, atom.residue.id, atom.residue.name)
                       for atom in topology_atoms[:protein_atom_count]}
        if set(local_flexible_residues) - represented or not local_mobile:
            raise ValueError("local flexible residue selection does not match the receptor")
    local_mobile_set = set(local_mobile)
    if not move_protein:
        for atom_index in range(protein_atom_count):
            atom = topology_atoms[atom_index]
            water_hydrogen = (atom.residue.name in {"HOH", "WAT"} and
                              atom.element.symbol.upper() == "H")
            if atom_index not in local_mobile_set and not (relax_water_hydrogens and water_hydrogen):
                system.setParticleMass(atom_index, 0.0)
    if fix_ligand_all:
        for atom_index in ligand_atom_indices:
            system.setParticleMass(atom_index, 0.0)
    if fix_ligand_heavy:
        for atom_index in ligand_atom_indices:
            if topology_atoms[atom_index].element.symbol.upper() != "H":
                system.setParticleMass(atom_index, 0.0)

    positions = np.asarray(modeller.positions.value_in_unit(unit.nanometer), dtype=float)
    initial_positions_nm = positions.copy()
    restraint = openmm.CustomExternalForce(
        "0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
    )
    restraint.addGlobalParameter(
        "k", restraint_k_kj_mol_a2 * 100.0 * unit.kilojoule_per_mole / unit.nanometer**2
    )
    for label in ("x0", "y0", "z0"):
        restraint.addPerParticleParameter(label)
    for atom_index in _restraint_atom_indices(topology_atoms, ligand_atom_indices, restrain_ligand_hydrogens):
        restraint.addParticle(atom_index, initial_positions_nm[atom_index])
    system.addForce(restraint)
    if local_mobile:
        receptor_restraint = openmm.CustomExternalForce(
            "0.5*k_receptor*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        receptor_restraint.addGlobalParameter("k_receptor", receptor_restraint_k_kj_mol_a2 * 100.0)
        for label in ("x0", "y0", "z0"):
            receptor_restraint.addPerParticleParameter(label)
        for atom_index in local_mobile:
            if topology_atoms[atom_index].element.symbol.upper() != "H":
                receptor_restraint.addParticle(atom_index, initial_positions_nm[atom_index])
        system.addForce(receptor_restraint)

    integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
    platform = openmm.Platform.getPlatformByName("CPU")
    context = openmm.Context(system, integrator, platform)
    context.setPositions(initial_positions_nm)
    initial_energy = context.getState(getEnergy=True).getPotentialEnergy()
    openmm.LocalEnergyMinimizer.minimize(
        context,
        tolerance=tolerance_kj_mol_nm * unit.kilojoule_per_mole / unit.nanometer,
        maxIterations=0,
    )
    state = context.getState(getPositions=True, getEnergy=True)
    final_energy = state.getPotentialEnergy()
    minimized_nm = np.asarray(state.getPositions().value_in_unit(unit.nanometer), dtype=float)
    final_energy_kj = float(final_energy.value_in_unit(unit.kilojoule_per_mole))
    initial_energy_kj = float(initial_energy.value_in_unit(unit.kilojoule_per_mole))
    if final_energy_kj > initial_energy_kj + max(1.0, abs(initial_energy_kj) * 1e-7):
        raise RuntimeError(
            "restrained minimization did not lower the complex potential energy: "
            f"{initial_energy_kj:.6f} -> {final_energy_kj:.6f} kJ/mol"
        )
    displacement_a = np.linalg.norm(
        (minimized_nm - initial_positions_nm)[ligand_atom_indices] * 10.0, axis=1
    )
    ligand_initial = initial_positions_nm[ligand_atom_indices] * 10.0
    ligand_final = minimized_nm[ligand_atom_indices] * 10.0
    ligand_initial -= ligand_initial.mean(axis=0)
    ligand_final -= ligand_final.mean(axis=0)
    rmsd_a = math.sqrt(np.mean(np.sum((ligand_final - ligand_initial) ** 2, axis=1)))

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    complex_pdb = output / "minimized_complex.pdb"
    with complex_pdb.open("w", encoding="utf-8") as handle:
        app.PDBFile.writeFile(modeller.topology, state.getPositions(), handle, keepIds=True)
    provenance = {
        "protocol": "SophosQM open-source restrained minimization substitute",
        "protein_force_field": "AMBER14 protein.ff14SB",
        "ligand_force_field": "OpenFF 2.2.0",
        "ligand_charge_model": "openff-gnn-am1bcc-1.0.0.pt",
        "water_model": "TIP3P",
        "protein_input": str(protein_pdb),
        "ligand_sdf": str(ligand_sdf),
        "ligand_name": ligand_name,
        "restraint": "harmonic positional restraint on all ligand atoms",
        "restraint_k_kj_mol_angstrom2": restraint_k_kj_mol_a2,
        "openmm_version": openmm.__version__,
        "initial_potential_energy_kj_mol": initial_energy_kj,
        "final_potential_energy_kj_mol": final_energy_kj,
        "ligand_rmsd_angstrom": rmsd_a,
        "ligand_max_displacement_angstrom": float(displacement_a.max()),
        "substitution": "OpenMM/OpenFF replaces unavailable Schrodinger MacroModel/OPLS3e",
        "protein_coordinate_policy": "free" if move_protein else "fixed_prepared_receptor_and_waters",
        "retained_receptor_component_models": receptor_models,
    }
    if starting_ligand_complex is not None:
        provenance["starting_ligand_complex"] = str(Path(starting_ligand_complex).resolve())
        provenance["starting_ligand_complex_sha256"] = hashlib.sha256(
            Path(starting_ligand_complex).read_bytes()).hexdigest()
    if local_mobile:
        fixed_indices = [index for index in range(protein_atom_count) if index not in local_mobile_set]
        fixed_displacement = float(np.max(np.linalg.norm(
            (minimized_nm-initial_positions_nm)[fixed_indices] * 10.0, axis=1)))
        if fixed_displacement > 0.0001:
            raise ValueError("frozen backbone, waters or remote receptor atoms moved")
        provenance.update(protein_coordinate_policy="restrained_local_sidechains_fixed_backbone_waters_remote",
            local_flexible_residues=[list(residue) for residue in local_flexible_residues],
            n_mobile_receptor_atoms=len(local_mobile),
            receptor_restraint_k_kj_mol_angstrom2=receptor_restraint_k_kj_mol_a2,
            max_fixed_receptor_displacement_angstrom=fixed_displacement,
            max_mobile_receptor_displacement_angstrom=float(np.max(np.linalg.norm(
                (minimized_nm-initial_positions_nm)[local_mobile] * 10.0, axis=1))))
    if not restrain_ligand_hydrogens:
        provenance["restraint"] = "harmonic positional restraint on ligand heavy atoms only"
        provenance["ligand_hydrogen_positional_restraints"] = False
        provenance["n_position_restrained_ligand_atoms"] = restraint.getNumParticles()
    if relax_water_hydrogens:
        fixed_indices = [index for index, atom in enumerate(topology_atoms[:protein_atom_count])
                         if not (atom.residue.name in {"HOH", "WAT"} and
                                 atom.element.symbol.upper() == "H")]
        fixed_displacement = max(np.linalg.norm(minimized_nm[index] - initial_positions_nm[index]) * 10.0
                                 for index in fixed_indices)
        if fixed_displacement > 0.0001:
            raise ValueError("non-water-H receptor coordinates changed")
        provenance["protein_coordinate_policy"] = "fixed_protein_and_water_oxygens_mobile_water_hydrogens"
        provenance["n_mobile_water_hydrogens"] = protein_atom_count - len(fixed_indices)
        provenance["max_fixed_receptor_displacement_angstrom"] = float(fixed_displacement)
    if fix_ligand_all:
        if displacement_a.max() > 0.0001:
            raise ValueError("fully fixed ligand moved during water hydrogen relaxation")
        provenance["ligand_coordinate_policy"] = "all_atoms_fixed_to_starting_pose"
    if fix_ligand_heavy:
        heavy_displacement = max(displacement_a[atom.GetIdx()]
                                 for atom in rdkit_ligand.GetAtoms()
                                 if atom.GetSymbol().upper() != "H")
        if heavy_displacement > 0.0001:
            raise ValueError("fixed ligand heavy atom moved during hydrogen minimization")
        provenance["ligand_coordinate_policy"] = "fixed_heavy_atoms_hydrogens_only"
        provenance["max_ligand_heavy_displacement_angstrom"] = float(heavy_displacement)
    provenance_json = output / "minimization.json"
    provenance_json.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    return MinimizationResult(
        complex_pdb, provenance_json,
        provenance["initial_potential_energy_kj_mol"],
        provenance["final_potential_energy_kj_mol"], rmsd_a,
        provenance["ligand_max_displacement_angstrom"],
    )
