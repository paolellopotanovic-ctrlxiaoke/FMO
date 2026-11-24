"""Build label-independent finite-shell QuanPol endpoint diagnostics, not FEP."""
from __future__ import annotations

import argparse
from collections import defaultdict
import io
import json
from pathlib import Path
import tempfile

import numpy as np

from analyze_endpoint_internal_energy import physical_atom_maps, verify_bonds
from analyze_openfe_endpoints import minimum_image, read_topology, unwrap_ligand
from run_matched_openfe import digest, write_json


def whole_components(coordinates, bonds, box):
    adjacency = defaultdict(list)
    for first, second in bonds:
        adjacency[first].append(second)
        adjacency[second].append(first)
    whole = np.array(coordinates, copy=True)
    visited = set()
    for anchor in range(len(coordinates)):
        if anchor in visited:
            continue
        visited.add(anchor)
        pending = [anchor]
        while pending:
            first = pending.pop()
            for second in adjacency[first]:
                if second not in visited:
                    whole[second] = whole[first] + minimum_image(coordinates[second] - coordinates[first], box)
                    visited.add(second)
                    pending.append(second)
    return whole


def build(args):
    from gufe import LigandAtomMapping
    from netCDF4 import Dataset
    from openmm import app
    from rdkit import Chem

    args.output.mkdir(parents=True, exist_ok=True)
    alias = Path(tempfile.mkdtemp(prefix="qp-endpoint-"))
    for name, source in (("qp", args.gamess / "auxdata/QUANPOL"),
                         ("sk", args.gamess / "auxdata/DFTB/3OB-3-1")):
        (alias / name).symlink_to(source.resolve(), target_is_directory=True)
    diagnostics_path = args.run / "endpoint_structure_diagnostics_complete_20261026.json"
    diagnostics = json.loads(diagnostics_path.read_text())
    if diagnostics["affinity_used_for_frame_selection"]:
        raise ValueError("label-independent input required")
    tasks = []
    for leg in ("bound", "solution"):
        root = args.run / "production" / args.target / "edge_000" / leg
        topology_path = sorted((root / "shared").glob("*/hybrid_system.pdb"))[0]
        atoms = read_topology(topology_path)
        topology_lines = topology_path.read_text().splitlines()
        ligand_count = 0
        for offset, line in enumerate(topology_lines):
            if line.startswith(("ATOM  ", "HETATM")) and line[17:20] == "UNK":
                topology_lines[offset] = line[:12] + f"X{ligand_count:03d}" + line[16:]
                ligand_count += 1
        pdb = app.PDBFile(io.StringIO("\n".join(topology_lines)))
        if pdb.topology.getNumAtoms() != len(atoms):
            raise ValueError("PDB reader changed hybrid atom inventory")
        mapping = LigandAtomMapping.from_json(root / "mapping.json")
        indices, molecules = physical_atom_maps(atoms, mapping)
        residues = list(pdb.topology.residues())
        pdb_atoms = list(pdb.topology.atoms())
        for endpoint in (0, 10):
            repeats = [row for row in diagnostics["targets"][args.target][leg]["complete_repeats"]
                       if row["endpoint"] == endpoint]
            repeats = sorted(repeats, key=lambda row: row["trajectory"])
            if args.repeat < 0 or args.repeat >= len(repeats):
                raise ValueError(f"repeat index {args.repeat} outside 0..{len(repeats)-1} for {leg} endpoint {endpoint}")
            repeat = repeats[args.repeat]
            frame = next(row for row in repeat["frames"] if row["iteration"] == args.iteration)
            with Dataset(repeat["trajectory"]) as data:
                if int(data["states"][args.iteration, frame["replica"]]) != endpoint:
                    raise ValueError("endpoint replica mismatch")
                saved = np.asarray(data["analysis_particle_indices"][:], dtype=int)
                reference = np.array([atom["xyz"] for atom in atoms])
                setup = np.load(topology_path.with_name("hybrid_positions.npy")) * 10
                if not np.allclose(setup[saved], reference, atol=.002):
                    raise ValueError("saved coordinate/PDB mapping mismatch")
                raw = data["positions"][args.iteration, frame["replica"]]
                if np.ma.getmaskarray(raw).any():
                    raise ValueError("frame not saved")
                coordinates = np.asarray(raw, dtype=float) * 10
                box = np.asarray(data["box_vectors"][args.iteration, frame["replica"]]) * 10
            ligand = indices[endpoint]
            molecule = molecules[endpoint]
            bonds = verify_bonds(topology_path, ligand, molecule)
            ligand_xyz = unwrap_ligand(coordinates[ligand], bonds, box)
            origin = ligand_xyz.mean(axis=0)
            ligand_xyz -= origin
            centered = minimum_image(coordinates - origin, box)
            topology_bonds = [(first.index, second.index) for first, second in pdb.topology.bonds()]
            centered = whole_components(centered, topology_bonds, box)
            selected_residues = []
            for residue in residues:
                selected = [atom.index for atom in residue.atoms()]
                if residue.name == "UNK":
                    continue
                heavy = [index for index in selected if atoms[index]["element"].upper() != "H"]
                is_solvent = residue.name.upper() in ("HOH", "WAT", "NA", "CL")
                if not is_solvent or np.linalg.norm(centered[heavy, None] - ligand_xyz[None, :], axis=-1).min() <= args.shell:
                    selected_residues.append((residue, selected))
            selected_residues.sort(key=lambda item: item[0].name in ("HOH", "WAT"))
            lines = []
            inventory = []
            serial = 0
            for number, (residue, selected) in enumerate(selected_residues, 1):
                selected.sort(key=lambda index: atoms[index]["element"].upper() == "H")
                inventory.append({"chain": residue.chain.id, "residue": residue.id,
                                  "name": residue.name, "indices": selected})
                for index in selected:
                    serial += 1
                    atom = pdb_atoms[index]
                    name = atom.name if len(atom.name) == 4 else f" {atom.name:<3}"
                    chain = residue.chain.id[:1] or "A"
                    element = atom.element.symbol
                    position = centered[index]
                    lines.append(f"ATOM  {serial:5d} {name:4s} {residue.name:3s} {chain}{number:4d}    "
                                 f"{position[0]:8.3f}{position[1]:8.3f}{position[2]:8.3f}  1.00  0.00          {element:>2}")
            symbols = [atom.GetSymbol() for atom in molecule.GetAtoms()]
            species = list(dict.fromkeys(symbols))
            hubbard = {"C": -.1492, "H": -.1857, "N": -.1535, "O": -.1575, "S": -.11,
                       "P": -.14, "F": -.1623, "Cl": -.0697, "Br": -.0573, "I": -.0433}
            sk_lines = []
            for first in species:
                for second in species:
                    path = alias / "sk" / f"{first}-{second}.skf"
                    if not path.is_file():
                        raise ValueError(f"missing parameter {path}")
                    sk_lines.append(f' {first.upper()} {second.upper()} "{path}"')
            charge = Chem.GetFormalCharge(molecule)
            header = [f" $CONTRL SCFTYP=RHF RUNTYP=ENERGY ICHARG={charge} MAXIT=200 $END",
                      " $BASIS GBASIS=DFTB $END", " $SCF SOSCF=.T. DIRSCF=.T. FDIFF=.F. $END",
                      " $SYSTEM TIMLIM=20 MWORDS=100 $END",
                      " $DFTB NDFTB=3 DAMPXH=.T. DAMPEX=4.0",
                      " HUBDER(1)=" + " ".join(str(hubbard[symbol]) for symbol in species) + " $END",
                      " $DFTBSK", *sk_lines, " $END",
                      " $QUANPO NSTEP=0 NFFTYP=30000 NFFFILE=3 IDOPOL=0",
                      " ITYPWAT=305 DFTBMM=0.70 LJQMMM=1 SWRA=10 SWRB=12 ISWITCH=2",
                      f" TOPAMIA='{alias}/qp/amino10.in'",
                      f" TOPNTER='{alias}/qp/aminont10.in'",
                      f" TOPCTER='{alias}/qp/aminoct10.in'",
                      f" PARFILE='{alias}/qp/parm10.dat'",
                      f" PARFIL2='{alias}/qp/frcmod.ff14SB'", " $END", " $DATA",
                      f"{args.target} {leg} endpoint {endpoint} finite shell diagnostic", "C1"]
            qm_lines = [f" {symbol} {atom.GetAtomicNum():.1f} " + " ".join(f"{value:.10f}" for value in position)
                        for symbol, atom, position in zip(symbols, molecule.GetAtoms(), ligand_xyz, strict=True)]
            directory = args.output / leg / f"endpoint_{endpoint:02d}"
            directory.mkdir(parents=True, exist_ok=True)
            deck = directory / f"{args.target.lower()}_{leg}_{endpoint:02d}.inp"
            deck.write_text("\n".join([*header, *qm_lines, " $END", " $FFPDB", *lines, " $END", ""]))
            (directory / "environment.pdb").write_text("\n".join([*lines, "END", ""]))
            environment_indices = [index for row in inventory for index in row["indices"]]
            np.savez(directory / "physical_coordinates.npz", ligand_xyz=ligand_xyz,
                     environment_xyz=centered[environment_indices],
                     environment_original_indices=saved[environment_indices])
            task = {"leg": leg, "endpoint": endpoint, "deck": str(deck.resolve()),
                    "deck_sha256": digest(deck), "trajectory": repeat["trajectory"],
                    "iteration": args.iteration, "replica": frame["replica"],
                    "topology": str(topology_path.resolve()), "topology_sha256": digest(topology_path),
                    "mapping_sha256": digest(root / "mapping.json"), "qm_atoms": len(symbols),
                    "mm_atoms": serial, "formal_charge": charge, "residues": inventory,
                    "box_angstrom": box.tolist(), "ligand_xyz_angstrom": ligand_xyz.tolist()}
            write_json(directory / "input_manifest.json", task)
            tasks.append(task)
            print(leg, endpoint, "QM", len(symbols), "MM", serial, flush=True)
    write_json(args.output / "plan.json", {"status": "INPUTS_BUILT", "target": args.target,
        "labels_used": False, "repeat_index": args.repeat,
        "frame_rule": "lexically indexed diagnostic repeat per endpoint, fixed iteration",
        "shell_angstrom": args.shell, "alias": str(alias), "tasks": tasks,
        "limitations": ["Finite solvent shell; not the periodic OpenFE Hamiltonian.",
                        "All protein residues retained; no receptor covalent cuts.",
                        "DFTB3/AMBER ff14SB QM/MM diagnostic, not a binding free energy or FEP correction.",
                        "QuanPol atom typing and energy closure must be audited before any interpretation."]})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--gamess", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--target", default="CDK8")
    parser.add_argument("--iteration", type=int, default=2000)
    parser.add_argument("--repeat", type=int, default=0)
    parser.add_argument("--shell", type=float, default=12.)
    args = parser.parse_args()
    if args.shell <= 0 or not np.isfinite(args.shell):
        parser.error("finite positive shell required")
    build(args)


if __name__ == "__main__":
    main()
