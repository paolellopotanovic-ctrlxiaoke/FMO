"""Explicit source-state tartrate model; never infer neutral CCD protonation."""
import json
from pathlib import Path

import numpy as np
from rdkit import Chem

from .component_inventory import inventory, sha


HEAVY = {"C1": "C", "C2": "C", "C3": "C", "C4": "C",
         "O1": "O", "O11": "O", "O2": "O", "O3": "O", "O4": "O", "O41": "O"}
BONDS = [("C1", "O1", 2), ("C1", "O11", 1), ("C1", "C2", 1),
         ("C2", "O2", 1), ("C2", "C3", 1), ("C3", "O3", 1),
         ("C3", "C4", 1), ("C4", "O4", 2), ("C4", "O41", 1)]


def source_tartrate(record):
    if record["id"].split(":")[1] != "TLA" or any(atom.get("altloc") for atom in record["atoms"]):
        raise ValueError("tartrate requires one resolved source TLA residue")
    atoms = record["atoms"]
    if any(not {"name", "element", "xyz"} <= set(atom) for atom in atoms):
        raise ValueError("tartrate source atom identity/coordinates are incomplete; do not delete or mutate")
    heavy = [atom for atom in atoms if atom["element"] != "H"]
    if len(heavy) != 10 or {atom["name"]: atom["element"] for atom in heavy} != HEAVY:
        raise ValueError("tartrate requires complete, unique named heavy atoms")
    positions = np.asarray([atom["xyz"] for atom in atoms])
    names = [atom["name"] for atom in atoms]
    indices = {atom["name"]: index for index, atom in enumerate(atoms) if atom["element"] != "H"}
    molecule = Chem.RWMol()
    for atom in atoms:
        current = Chem.Atom(atom["element"])
        current.SetNoImplicit(True)
        molecule.AddAtom(current)
    for first, second, order in BONDS:
        distance = float(np.linalg.norm(positions[indices[first]] - positions[indices[second]]))
        if not 1.10 <= distance <= 1.75:
            raise ValueError("tartrate source heavy bond geometry is not supported")
        molecule.AddBond(indices[first], indices[second], Chem.BondType.DOUBLE if order == 2 else Chem.BondType.SINGLE)
    attached = {name: [] for name in HEAVY}
    aliases = []
    for index, atom in enumerate(atoms):
        if atom["element"] != "H":
            continue
        candidates = sorted((float(np.linalg.norm(positions[index] - positions[partner])), name, partner)
                            for name, partner in indices.items())
        distance, parent, partner = candidates[0]
        if not .7 <= distance <= 1.25 or candidates[1][0] - distance < .2:
            raise ValueError("tartrate source hydrogen parent is ambiguous")
        attached[parent].append(index)
        names[index] = f"H{index + 1}"
        aliases.append({"source_index": index, "source_name": atom["name"], "prepared_name": names[index], "parent": parent})
        molecule.AddBond(index, partner, Chem.BondType.SINGLE)
    if any(len(attached[name]) != 1 for name in ("C2", "C3", "O2", "O3")):
        raise ValueError("tartrate requires explicit complete carbon and alcohol hydrogens")
    if any(attached[name] for name in ("C1", "C4", "O1", "O4")):
        raise ValueError("unsupported tartrate proton position")
    for name in ("O11", "O41"):
        if len(attached[name]) > 1:
            raise ValueError("invalid tartrate carboxyl hydrogen inventory")
        if not attached[name]:
            molecule.GetAtomWithIdx(indices[name]).SetFormalCharge(-1)
    molecule = molecule.GetMol()
    Chem.SanitizeMol(molecule)
    conformer = Chem.Conformer(len(atoms))
    for index, point in enumerate(positions):
        conformer.SetAtomPosition(index, point)
        molecule.GetAtomWithIdx(index).SetProp("_TriposAtomName", names[index])
    molecule.AddConformer(conformer)
    Chem.AssignStereochemistryFrom3D(molecule)
    if any(label == "?" for _, label in Chem.FindMolChiralCenters(molecule, includeUnassigned=True)):
        raise ValueError("tartrate stereochemistry unresolved")
    molecule.SetProp("_Name", record["id"])
    return molecule, names, aliases


def restore_tartrate(raw, prepared):
    raw, prepared = Path(raw), Path(prepared)
    records = [row for row in inventory(raw)["residues"] if row["kind"] not in {"standard_residue", "water"}]
    if not records or any(row["id"].split(":")[1] != "TLA" for row in records):
        raise ValueError("tartrate-only restoration cannot drop other special components")
    lines = [line for line in prepared.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM"))]
    molecules, models, connections = [], [], []
    for record in records:
        molecule, names, aliases = source_tartrate(record)
        chain, residue, number = record["id"].split(":")
        if not number.isdigit() or chain == "_":
            raise ValueError("tartrate insertion/blank-chain identity requires a validated mapping")
        if any(line[21:22].strip() == chain and line[22:27].strip() == number for line in lines):
            raise ValueError("tartrate residue identity collides with prepared receptor")
        offset = len(lines) + 1
        for atom, name in zip(record["atoms"], names, strict=True):
            point = atom["xyz"]
            lines.append(f"HETATM{len(lines)+1:5d} {name:>4} {residue:>3} {chain}{int(number):4d}    "
                         f"{point[0]:8.3f}{point[1]:8.3f}{point[2]:8.3f}  1.00  0.00          {atom['element']:>2}")
        for bond in molecule.GetBonds():
            connections.append(f"CONECT{offset+bond.GetBeginAtomIdx():5d}{offset+bond.GetEndAtomIdx():5d}")
        molecules.append(molecule)
        models.append({"id": record["id"], "formal_charge": Chem.GetFormalCharge(molecule), "aliases": aliases,
                       "source_state_preserved_not_assay_population": True})
    serial_lines = [line[:6] + f"{index:5d}" + line[11:] for index, line in enumerate(lines, 1)]
    text = "\n".join(serial_lines + connections + ["END"]) + "\n"
    sdf = prepared.with_suffix(".components.sdf")
    with Chem.SDWriter(str(sdf)) as writer:
        for molecule in molecules:
            writer.write(molecule)
    prepared.write_text(text)
    report = {"schema": "source-tartrate-model-v1", "raw_sha256": sha(raw), "prepared_sha256": sha(prepared),
        "components_sdf_sha256": sha(sdf), "models": models, "experimental_labels_used": False,
        "small_molecule_forcefield": "openff-2.2.0", "charge_model": "openff-gnn-am1bcc-1.0.0.pt",
        "hydrogen_policy": "preserve explicit source atoms and validate parents; only carboxyl O lacking source H gets -1; not a pKa prediction",
        "qm_representation_validated": False, "scope": "explicit source-state MM model; in-pocket QM unsupported"}
    prepared.with_suffix(".components.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def component_molecules(prepared):
    from openff.toolkit import Molecule
    from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
    prepared = Path(prepared)
    metadata = prepared.with_suffix(".components.json")
    if not metadata.exists():
        return [], []
    report = json.loads(metadata.read_text())
    sdf = prepared.with_suffix(".components.sdf")
    if report["prepared_sha256"] != sha(prepared) or report["components_sdf_sha256"] != sha(sdf):
        raise ValueError("registered receptor component model changed")
    molecules = []
    for molecule in Chem.SDMolSupplier(str(sdf), removeHs=False):
        if molecule is None:
            raise ValueError("unreadable receptor component model")
        current = Molecule.from_rdkit(molecule, allow_undefined_stereo=False, hydrogens_are_explicit=True)
        current.assign_partial_charges(partial_charge_method="openff-gnn-am1bcc-1.0.0.pt", toolkit_registry=NAGLToolkitWrapper())
        molecules.append(current)
    if len(molecules) != len(report["models"]):
        raise ValueError("component model membership differs")
    return molecules, report["models"]
