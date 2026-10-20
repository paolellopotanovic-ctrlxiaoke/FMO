"""Explicit retained-ion/phosphoresidue models; no generic metal guesses."""
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np

from .component_inventory import inventory


ION_TEMPLATES = {"NA": ("NA", "Na", 1), "CA": ("CA", "Ca", 2)}
PHOSPHO_RESIDUES = {"TPO", "SEP", "PTR"}


def parameter_paths():
    import openmmforcefields
    base = Path(openmmforcefields.__file__).parent / "ffxml/amber"
    return [base / name for name in ("ff14SB.xml", "phosaa14SB.xml", "tip3p_standard.xml")]


def validate_phosphoresidue_identity(record):
    name = record["id"].split(":")[1]
    root = ET.parse(parameter_paths()[1]).getroot()
    template = next((row for row in root.find("Residues") if row.attrib["name"] == name), None)
    if template is None:
        raise ValueError("unsupported phosphoresidue identity")
    elements = {row.attrib["name"]: row.attrib["element"].upper()
                for path in parameter_paths() for row in ET.parse(path).getroot().find("AtomTypes")}
    expected = {atom.attrib["name"]: elements[atom.attrib["type"]] for atom in template.findall("Atom")}
    aliases = {"1HG2": "HG21", "2HG2": "HG22", "3HG2": "HG23"} if name == "TPO" else {}
    actual = {aliases.get(atom["name"], atom["name"]): atom["element"] for atom in record["atoms"]}
    if len(actual) != len(record["atoms"]) or actual != expected:
        raise ValueError("complete atom/element-identical phosphoresidue template required")


def write_phospho_bonds(destination):
    source = ET.parse(parameter_paths()[1]).getroot()
    root = ET.Element("Residues")
    for residue in source.find("Residues"):
        if residue.attrib["name"] not in PHOSPHO_RESIDUES:
            continue
        target = ET.SubElement(root, "Residue", name=residue.attrib["name"])
        for bond in residue.findall("Bond"):
            ET.SubElement(target, "Bond", **{"from": bond.attrib["atomName1"], "to": bond.attrib["atomName2"]})
        ET.SubElement(target, "Bond", **{"from": "-C", "to": "N"})
        ET.SubElement(target, "Bond", **{"from": "C", "to": "+N"})
    ET.ElementTree(root).write(destination, encoding="unicode")
    return Path(destination)


def registered_phosphoresidue_model(path):
    from .component_inventory import sha
    model_path = Path(path).with_suffix(".model.json")
    if not model_path.exists():
        return None
    model = json.loads(model_path.read_text())
    if model.get("schema") != "retained-special-chemistry-model-v1":
        return None
    if model.get("prepared_sha256") != sha(path):
        raise ValueError("registered special receptor fingerprint changed")
    expected_parameters = {str(parameter): sha(parameter) for parameter in parameter_paths()}
    if model.get("parameter_sha256") != expected_parameters:
        raise ValueError("registered special receptor parameters changed")
    records = {row["id"]: row for row in inventory(path)["residues"]}
    templates = model.get("matched_special_templates", {})
    if not templates:
        raise ValueError("registered special receptor has no validated templates")
    for identifier, template in templates.items():
        record = records.get(identifier)
        if record is None:
            raise ValueError("registered special receptor component missing")
        name = identifier.split(":")[1]
        if name in PHOSPHO_RESIDUES:
            validate_phosphoresidue_identity(record)
            expected_charge = -2
        elif name in ION_TEMPLATES and len(record["atoms"]) == 1:
            expected_charge = ION_TEMPLATES[name][2]
        else:
            raise ValueError("unregistered special receptor template")
        if template["template"] != name or abs(template["charge"] - expected_charge) > 1e-6:
            raise ValueError("registered special receptor template charge changed")
    return model


def load_phosphoresidue_bonds():
    from openmm.app import Topology
    from tempfile import TemporaryDirectory
    with TemporaryDirectory(prefix="fmo-phospho-bonds-") as temporary:
        Topology.loadBondDefinitions(str(write_phospho_bonds(Path(temporary) / "bonds.xml")))


def assemble_retained(raw, prepared, output):
    from openmm.app import PDBFile, ForceField, Topology, NoCutoff
    from openmm import unit
    raw, prepared, output = Path(raw), Path(prepared), Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    graph = inventory(raw)
    records = [row for row in graph["residues"] if row["kind"] not in {"standard_residue", "water"}]
    additions = {}
    for record in records:
        elements = record["elements"]
        name = record["id"].split(":")[1]
        if record["kind"] == "metal" and len(record["atoms"]) == 1 and next(iter(elements)) in ION_TEMPLATES:
            element = next(iter(elements))
            residue_name, symbol, charge = ION_TEMPLATES[element]
            additions[record["id"]] = {"model": "Amber TIP3P fixed nonbonded ion; no artificial metal bonds",
                "element": symbol, "formal_charge": charge, "target_name": residue_name}
        elif name in PHOSPHO_RESIDUES:
            additions[record["id"]] = {"model": "Amber phosaa14SB source explicit proton state",
                "target_name": name, "formal_charge": -2,
                "proton_population_or_assay_state_verified": False}
        else:
            raise ValueError(f"no validated retained-component implementation for {record['id']}; do not discard it")
    rows = [line for line in prepared.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM"))]
    raw_rows = [line for line in raw.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM"))]
    retained = []
    for identifier, model in additions.items():
        chain, name, number = identifier.split(":")
        chain = "" if chain == "_" else chain
        group = [line for line in raw_rows if line[21:22].strip() == chain and line[17:20].strip() == name and line[22:27].strip() == number]
        if any(line[16:17].strip() for line in group):
            raise ValueError("retained component has unresolved alternate locations")
        if any(line[21:22].strip() == chain and line[22:27].strip() == number for line in rows):
            raise ValueError("retained component collides with existing prepared residue identity")
        if name in PHOSPHO_RESIDUES:
            if not number.isdigit():
                raise ValueError("phosphoresidue insertion requires explicit polymer insertion-code mapping")
            following = [index for index, line in enumerate(rows) if line[21:22].strip() == chain and line[22:26].strip().isdigit() and int(line[22:26]) > int(number) and line.startswith("ATOM")]
            if not following:
                raise ValueError("terminal phosphoresidue requires dedicated terminal template")
            aliases = {"1HG2": "HG21", "2HG2": "HG22", "3HG2": "HG23"} if name == "TPO" else {}
            renamed = []
            for position, line in enumerate(group):
                atom_name = line[12:16].strip()
                replacement = aliases.get(atom_name, atom_name)
                if replacement != atom_name:
                    renamed.append({"source": atom_name, "model": replacement})
                line = "ATOM  " + line[6:12] + f"{replacement:>4}" + line[16:]
                group[position] = line
            model["hydrogen_name_aliases"] = renamed
            rows[following[0]:following[0]] = group
        else:
            replacement = model["target_name"]
            line = group[0]
            line = "HETATM" + line[6:12] + f"{replacement:>4}" + line[16:17] + f"{replacement:>3}" + line[20:]
            rows.append(line)
        retained.append({"raw_id": identifier, **model, "n_raw_atoms_retained": len(group)})
    text = "\n".join(line[:6] + f"{serial:5d}" + line[11:] for serial, line in enumerate(rows, 1)) + "\nEND\n"
    output.write_text(text)
    definitions = write_phospho_bonds(output.with_suffix(".bonds.xml"))
    Topology.loadBondDefinitions(str(definitions))
    pdb = PDBFile(str(output))
    forcefield = ForceField(*map(str, parameter_paths()))
    templates = forcefield.getMatchingTemplates(pdb.topology)
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff, constraints=None, rigidWater=False)
    if system.getNumParticles() != len(rows):
        raise ValueError("retained component system particle count mismatch")
    charges = {}
    for residue, template in zip(pdb.topology.residues(), templates, strict=True):
        if residue.name in PHOSPHO_RESIDUES or len(list(residue.atoms())) == 1:
            charge = sum(atom.parameters["charge"] for atom in template.atoms)
            charges[f"{residue.chain.id.strip() or '_'}:{residue.name}:{residue.id}"] = {"template": template.name, "charge": charge}
    coordinates = np.asarray(pdb.positions.value_in_unit(unit.angstrom))
    if not np.isfinite(coordinates).all():
        raise ValueError("invalid prepared coordinates")
    result = {"schema": "retained-special-chemistry-model-v1", "raw_sha256": graph["source_sha256"],
        "baseline_prepared_sha256": hashlib.sha256(prepared.read_bytes()).hexdigest(),
        "prepared_sha256": hashlib.sha256(output.read_bytes()).hexdigest(), "retained": retained,
        "matched_special_templates": charges, "all_residues_template_matched": True,
        "forcefield_system_created": True, "n_particles": system.getNumParticles(),
        "parameter_sha256": {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in parameter_paths()},
        "experimental_labels_used": False, "coordination_bonds_added": False,
        "not_a_general_transition_metal_or_covalent_cofactor_model": True,
        "protonation_policy": "standard residues retain registered PDB2PQR states; source phosphorylated residue retained only if complete phosaa14SB template matches"}
    output.with_suffix(".model.json").write_text(json.dumps(result, indent=2) + "\n")
    return result
