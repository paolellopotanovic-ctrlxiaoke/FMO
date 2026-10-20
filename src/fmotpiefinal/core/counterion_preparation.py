"""Retain supported monatomic ions with complete Amber validation."""
import json
from pathlib import Path

from .component_inventory import inventory, sha
from .special_chemistry import parameter_paths


ION_CHARGES = {"F": -1, "CL": -1, "BR": -1, "I": -1, "NA": 1, "CA": 2}


def supported_monatomic_components(components):
    return bool(components) and all(len(row["atoms"]) == 1 and
        row["atoms"][0]["element"] in ION_CHARGES and not row["atoms"][0].get("altloc", "")
        for row in components)


def assemble_counterions(raw, prepared, output):
    from openmm.app import PDBFile, ForceField, NoCutoff
    raw, prepared, output = Path(raw), Path(prepared), Path(output)
    if output.exists():
        raise FileExistsError("counterion preparation requires a fresh receptor output")
    graph = inventory(raw)
    if graph["n_models"] != 1:
        raise ValueError("counterion preparation requires one registered source model")
    components = [row for row in graph["residues"] if row["kind"] not in {"standard_residue", "water"}]
    if not supported_monatomic_components(components):
        raise ValueError("counterion route only supports complete monatomic halides; other chemistry must not be dropped")
    ion_serials = {line[6:11].strip() for line in raw.read_text().splitlines()
        if line.startswith(("ATOM  ", "HETATM")) and any(
            (line[21:22].strip() or "_") + ":" + line[17:20].strip() + ":" + line[22:27].strip() == row["id"]
            for row in components)}
    ion_addresses = {(row["id"].split(":")[0], row["id"].split(":")[1], row["id"].split(":")[2]) for row in components}
    source_atoms = {line[6:11].strip(): line for line in raw.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM"))}
    coordination = {}
    for line in raw.read_text().splitlines():
        if line.startswith("CONECT"):
            serials = [line[start:start + 5].strip() for start in range(6, len(line), 5)]
            if serials and any(serials):
                for partner in filter(None, serials[1:]):
                    first = serials[0]
                    if not ion_serials.intersection({first, partner}):
                        continue
                    metal_serial, donor_serial = (first, partner) if first in ion_serials else (partner, first)
                    if metal_serial not in source_atoms or donor_serial not in source_atoms:
                        raise ValueError("explicit ion connectivity has unresolved atom identity")
                    metal, donor = source_atoms[metal_serial], source_atoms[donor_serial]
                    metal_element, donor_element = metal[76:78].strip().upper(), donor[76:78].strip().upper()
                    import math
                    distance = math.sqrt(sum((float(metal[start:start + 8]) - float(donor[start:start + 8])) ** 2 for start in (30, 38, 46)))
                    if metal_element not in {"NA", "CA"} or donor_element != "O" or not 1.8 <= distance <= 3.2:
                        raise ValueError("explicit ion connectivity outside registered Na/Ca oxygen-contact model")
                    coordination[(metal_serial, donor_serial)] = {"metal_serial": metal_serial, "donor_serial": donor_serial,
                        "metal_element": metal_element, "donor_element": donor_element,
                        "distance_angstrom": distance, "source_record": "CONECT",
                        "representation": "registered Amber nonbonded ion; no harmonic bond added"}
        if line.startswith("LINK  "):
            endpoints = [((line[chain:chain + 1].strip() or "_"), line[residue:residue + 3].strip(), line[number:number + 5].strip())
                         for chain, residue, number in ((21, 17, 22), (51, 47, 52))]
            if ion_addresses.intersection(endpoints):
                raise ValueError("explicit ion LINK requires a registered coordination model")
    lines = [line for line in prepared.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM"))]
    raw_lines = [line for line in raw.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM"))]
    retained = []
    for component in components:
        element = component["atoms"][0]["element"]
        chain, residue, number = component["id"].split(":")
        chain = "" if chain == "_" else chain
        matching = [line for line in raw_lines if line[21:22].strip() == chain and line[17:20].strip() == residue and line[22:27].strip() == number]
        if len(matching) != 1 or matching[0][16:17].strip():
            raise ValueError("ambiguous original counterion identity")
        if any(line[21:22].strip() == chain and line[22:27].strip() == number for line in lines):
            raise ValueError("counterion restoration collides with prepared residue identity")
        source = matching[0]
        lines.append("HETATM" + source[6:12] + f"{element:>4}" + source[16:17] + f"{element:>3}" + source[20:])
        retained.append({"raw_id": component["id"], "element": element, "template": element,
            "formal_charge": ION_CHARGES[element], "coordinates_retained": True,
            "model": "Amber TIP3P monatomic nonbonded ion; no inferred coordination bonds"})
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(line[:6] + f"{serial:5d}" + line[11:] for serial, line in enumerate(lines, 1)) + "\nEND\n")
    pdb = PDBFile(str(output))
    forcefield = ForceField(*map(str, parameter_paths()))
    templates = forcefield.getMatchingTemplates(pdb.topology)
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff, constraints=None, rigidWater=False)
    actual = {f"{residue.chain.id.strip() or '_'}:{residue.name}:{residue.id}": {"template": template.name,
        "charge": sum(atom.parameters["charge"] for atom in template.atoms)} for residue, template in
        zip(pdb.topology.residues(), templates, strict=True) if residue.name in ION_CHARGES}
    if len(actual) != len(retained) or any(abs(row["charge"] - ION_CHARGES[identifier.split(":")[1]]) > 1e-8 for identifier, row in actual.items()):
        raise ValueError("counterion template/charge validation failed")
    schema = "retained-monatomic-ion-model-v1" if any(row["element"] in {"NA", "CA"} for row in retained) else "retained-halide-model-v1"
    result = {"schema": schema, "raw_sha256": graph["source_sha256"],
        "baseline_prepared_sha256": sha(prepared), "prepared_sha256": sha(output), "retained": retained,
        "all_residues_template_matched": True, "forcefield_system_created": True, "n_particles": system.getNumParticles(),
        "matched_ion_templates": actual, "parameter_sha256": {str(path): sha(path) for path in parameter_paths()},
        "experimental_labels_used": False, "coordination_bonds_added": False,
        "source_coordination_contacts": list(coordination.values()),
        "not_a_general_cofactor_or_in_pocket_QM_model": True}
    output.with_suffix(".model.json").write_text(json.dumps(result, indent=2) + "\n")
    return result
