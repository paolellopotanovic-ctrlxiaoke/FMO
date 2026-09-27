"""Open-source restrained minimization replacing unavailable MacroModel/OPLS3e."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json
import math
import os
import tempfile

import numpy as np
import openmm
import openmm.app as app
import openmm.unit as unit
from openmmforcefields.generators import SystemGenerator
from openff.toolkit import Molecule, Topology
from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
from rdkit import Chem


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
    output = (
        Path(tempfile.gettempdir())
        / f"sophosqm-{source.stem}-{os.getpid()}-omm.pdb"
    )
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output


def _explicit_water_residues(topology: app.Topology) -> int:
    return sum(residue.name in {"HOH", "WAT", "TIP"} for residue in topology.residues())


def minimize_complex(protein_pdb: str | Path, ligand_sdf: str | Path, ligand_name: str,
                     output_dir: str | Path, restraint_k_kj_mol_a2: float = 100.0,
                     tolerance_kj_mol_nm: float = 1.0) -> MinimizationResult:
    """Minimize a prepared complex while restraining every ligand atom.

    The protocol keeps the PDB2PQR/PROPKA-prepared receptor, waters, and ligand
    chemically consistent while restraining every ligand atom to its input pose.
    """
    if restraint_k_kj_mol_a2 <= 0:
        raise ValueError("restraint force constant must be positive")
    rdkit_ligand = _selected_rdkit_ligand(ligand_sdf, ligand_name)
    off_ligand = Molecule.from_rdkit(rdkit_ligand, allow_undefined_stereo=False)
    off_ligand.assign_partial_charges(
        partial_charge_method="openff-gnn-am1bcc-1.0.0.pt",
        toolkit_registry=NAGLToolkitWrapper(),
    )
    protein_pdb_path = _water_compatible_pdb(protein_pdb)
    protein = app.PDBFile(str(protein_pdb_path))
    explicit_water_residues = _explicit_water_residues(protein.topology)
    ligand_topology = Topology.from_molecules([off_ligand]).to_openmm()
    ligand_coordinates = np.asarray(rdkit_ligand.GetConformer().GetPositions(), dtype=float) / 10.0
    modeller = app.Modeller(protein.topology, protein.positions)
    modeller.add(ligand_topology, ligand_coordinates * unit.nanometer)

    forcefield_files = ["amber14/protein.ff14SB.xml", "amber14/tip3p.xml"]
    generator = SystemGenerator(
        forcefields=forcefield_files,
        small_molecule_forcefield="openff-2.2.0",
        forcefield_kwargs={"constraints": None, "rigidWater": False},
        nonperiodic_forcefield_kwargs={"nonbondedMethod": app.NoCutoff},
        molecules=[off_ligand],
        cache=None,
    )
    system = generator.create_system(modeller.topology)
    topology_atoms = list(modeller.topology.atoms())
    protein_atom_count = len(list(protein.topology.atoms()))
    ligand_atom_indices = list(range(protein_atom_count, len(topology_atoms)))
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
    for atom_index in ligand_atom_indices:
        restraint.addParticle(atom_index, initial_positions_nm[atom_index])
    system.addForce(restraint)

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
        "explicit_water_residues": explicit_water_residues,
        "explicit_water_policy": "retained" if explicit_water_residues else "none",
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
        "protein_coordinate_policy": "free",
    }
    provenance_json = output / "minimization.json"
    provenance_json.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    return MinimizationResult(
        complex_pdb, provenance_json,
        provenance["initial_potential_energy_kj_mol"],
        provenance["final_potential_energy_kj_mol"], rmsd_a,
        provenance["ligand_max_displacement_angstrom"],
    )
