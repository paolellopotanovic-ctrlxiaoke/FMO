"""Retain chemically identified ACE/methylamide caps with full Amber validation."""
import json
import re
from pathlib import Path

import numpy as np

from .component_inventory import inventory, sha


ALIASES = {"ACE": {"1H": "HH31", "2H": "HH32", "3H": "HH33"},
           "NMA": {"CA": "CH3", "1HA": "HH31", "2HA": "HH32", "3HA": "HH33"}, "NME": {}}
CAP_ATOMS = {"ACE": {"CH3": "C", "C": "C", "O": "O", "HH31": "H", "HH32": "H", "HH33": "H"},
             "NME": {"N": "N", "H": "H", "CH3": "C", "HH31": "H", "HH32": "H", "HH33": "H"}}


def pdb_residue_fields(identifier):
    match = re.fullmatch(r"(-?\d+)([A-Za-z]?)", identifier)
    if match is None:
        raise ValueError("terminal cap requires a numeric PDB residue ID and optional insertion code")
    number = int(match[1])
    if not -999 <= number <= 9999:
        raise ValueError("terminal cap residue number exceeds PDB field width")
    return f"{number:4d}{match[2] or ' '}"


def cap_identity(record):
    name = record["id"].split(":")[1]
    if name == "UNK":
        return identify_unknown_cap(record)
    if name not in ALIASES:
        raise ValueError("unsupported cap identity; do not mutate other nonstandard residues")
    target = "NME" if name == "NMA" else name
    renamed = [ALIASES[name].get(atom["name"], atom["name"]) for atom in record["atoms"]]
    if (len(renamed) != len(set(renamed)) or
            dict(zip(renamed, [atom["element"] for atom in record["atoms"]], strict=True)) != CAP_ATOMS[target]):
        raise ValueError("complete atom/element-identical cap required; aliases do not authorize missing atoms")
    return target, dict(zip([atom["name"] for atom in record["atoms"]], renamed, strict=True))


def identify_unknown_cap(record):
    atoms = record["atoms"]
    if len(atoms) != 6 or len({atom["name"] for atom in atoms}) != 6:
        raise ValueError("unknown cap requires six complete, uniquely named atoms")
    if any("xyz" not in atom for atom in atoms):
        raise ValueError("unknown cap requires source geometry")
    positions = np.asarray([atom["xyz"] for atom in atoms], dtype=float)
    if not np.isfinite(positions).all():
        raise ValueError("unknown cap has non-finite coordinates")
    elements = [atom["element"] for atom in atoms]
    target = "ACE" if sorted(elements) == sorted(CAP_ATOMS["ACE"].values()) else "NME"
    if sorted(elements) != sorted(CAP_ATOMS[target].values()):
        raise ValueError("unknown component is not an atom-identical terminal cap")
    bounds = {("C", "C"): (1.35, 1.65), ("C", "O"): (1.1, 1.35),
              ("C", "N"): (1.3, 1.6), ("C", "H"): (.85, 1.25), ("H", "N"): (.7, 1.25)}
    neighbors = {index: set() for index in range(len(atoms))}
    for first in range(len(atoms)):
        for second in range(first + 1, len(atoms)):
            pair = tuple(sorted((elements[first], elements[second])))
            distance = float(np.linalg.norm(positions[first] - positions[second]))
            if distance < .65:
                raise ValueError("unknown cap contains overlapping atoms")
            lower, upper = bounds.get(pair, (0., 0.))
            if lower <= distance <= upper:
                neighbors[first].add(second)
                neighbors[second].add(first)
    signatures = {index: sorted(elements[other] for other in connected) for index, connected in neighbors.items()}
    expected = {"CH3": ("C", ["C", "H", "H", "H"]), "C": ("C", ["C", "O"]),
                "O": ("O", ["C"])} if target == "ACE" else {
                "CH3": ("C", ["H", "H", "H", "N"]), "N": ("N", ["C", "H"]), "H": ("H", ["N"])}
    mapping = {}
    for canonical, (element, connected) in expected.items():
        matches = [index for index in neighbors if elements[index] == element and signatures[index] == connected]
        if len(matches) != 1:
            raise ValueError("unknown cap connectivity is incomplete or ambiguous")
        mapping[matches[0]] = canonical
    methyl_hydrogens = sorted((index for index in neighbors if elements[index] == "H" and signatures[index] == ["C"]),
                             key=lambda index: atoms[index]["name"])
    if len(methyl_hydrogens) != 3:
        raise ValueError("unknown cap methyl hydrogen graph mismatch")
    mapping.update(zip(methyl_hydrogens, ("HH31", "HH32", "HH33"), strict=True))
    if len(mapping) != len(atoms):
        raise ValueError("unknown cap atom coverage mismatch")
    return target, {atoms[index]["name"]: canonical for index, canonical in mapping.items()}


def assemble_terminal_caps(raw, prepared, output):
    from openmm import unit
    from openmm.app import ForceField, NoCutoff, PDBFile
    raw, prepared, output = Path(raw), Path(prepared), Path(output)
    if output.exists():
        raise FileExistsError("terminal-cap repair requires a fresh output")
    graph = inventory(raw)
    if graph["n_models"] != 1:
        raise ValueError("one resolved source model required")
    specials = [row for row in graph["residues"] if row["kind"] not in {"standard_residue", "water"}]
    if not specials:
        raise ValueError("no source terminal caps to retain")
    models = {row["id"]: cap_identity(row) for row in specials}
    rows = [line for line in prepared.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM"))]
    raw_rows = [line for line in raw.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM"))]
    if any(line[16:17].strip() for line in rows + raw_rows):
        raise ValueError("resolve alternate locations before cap restoration")
    retained, changes, attachments = [], [], []
    for record in specials:
        source_chain, source_name, number = record["id"].split(":")
        source_chain = "" if source_chain == "_" else source_chain
        target, aliases = models[record["id"]]
        partner_name = "N" if target == "ACE" else "C"
        cap_atom = next(atom for atom in record["atoms"] if aliases[atom["name"]] == ("C" if target == "ACE" else "N"))
        source_matches = [(residue, atom) for residue in graph["residues"] if residue["kind"] == "standard_residue"
                          for atom in residue["atoms"] if atom["name"] == partner_name and
                          1.15 <= np.linalg.norm(np.array(cap_atom["xyz"]) - atom["xyz"]) <= 1.55]
        if len(source_matches) != 1:
            raise ValueError("ambiguous source cap attachment")
        source_terminal, source_partner = source_matches[0]
        matched = [line for line in rows if line.startswith("ATOM  ") and line[12:16].strip() == partner_name
                   and line[17:20].strip() == source_terminal["id"].split(":")[1] and
                   np.linalg.norm(np.array([float(line[start:start + 8]) for start in (30, 38, 46)]) - source_partner["xyz"]) <= .02]
        if len(matched) != 1:
            raise ValueError("cap has no unique source-coordinate-matched prepared polymer terminus")
        terminal = matched[0]
        chain = terminal[21:22].strip()
        chain_rows = [line for line in rows if line.startswith("ATOM  ") and line[21:22].strip() == chain
                      and line[17:20].strip() not in CAP_ATOMS]
        if terminal[21:27] != chain_rows[0 if target == "ACE" else -1][21:27]:
            raise ValueError("cap attachment is not a prepared chain terminus")
        terminal_key = terminal[21:27]
        terminal_name = terminal[17:20].strip()
        partners = [line for line in chain_rows if line[21:27] == terminal_key and line[12:16].strip() == partner_name]
        if len(partners) != 1:
            raise ValueError("ambiguous terminal cap attachment")
        partner = partners[0]
        point = np.array([float(partner[start:start + 8]) for start in (30, 38, 46)])
        distance = float(np.linalg.norm(point - np.array(cap_atom["xyz"])))
        if not 1.15 <= distance <= 1.55:
            raise ValueError("source cap is not covalently adjacent to prepared chain end")
        if np.linalg.norm(point - source_partner["xyz"]) > .02:
            raise ValueError("cap and prepared receptor reference frames differ")
        if any(line[21:22].strip() == chain and line[22:27].strip() == number for line in rows):
            used = {line[22:27].strip() for line in rows if line[21:22].strip() == chain}
            number = next((str(index) for index in range(1, 10000) if str(index) not in used), None)
            if number is None:
                raise ValueError("no available PDB residue identifier for retained cap")
        amide_hydrogen = []
        if target == "ACE":
            amide_hydrogen = [atom for atom in source_terminal["atoms"] if atom["element"] == "H" and
                .7 <= np.linalg.norm(point - atom["xyz"]) <= 1.25]
            expected_hydrogens = 0 if terminal_name == "PRO" else 1
            if len(amide_hydrogen) != expected_hydrogens:
                raise ValueError("raw capped N terminus needs an unambiguous amide proton state")
        revised = []
        for line in rows:
            if line[21:27] == terminal_key:
                name = line[12:16].strip()
                remove = name in ({"H", "H1", "H2", "H3"} if target == "ACE" else {"OXT"})
                if remove:
                    if target != "ACE" and any(atom["name"] == name for atom in source_terminal["atoms"]):
                        raise ValueError("terminal atom proposed for removal exists in capped raw source")
                    changes.append({"action": "RESTORE_CAPPED_AMIDE_H_STATE" if target == "ACE" else "REMOVE_UNCAPPED_TERMINAL_ADDITION",
                                    "residue": terminal_key.strip(), "atom": name})
                    continue
            revised.append(line)
        rows = revised
        for hydrogen in amide_hydrogen:
            raw_chain, raw_name, raw_number = source_terminal["id"].split(":")
            originals = [line for line in raw_rows if (line[21:22].strip() or "_") == raw_chain and
                         line[17:20].strip() == raw_name and line[22:27].strip() == raw_number and line[12:16].strip() == hydrogen["name"]]
            if len(originals) != 1:
                raise ValueError("source amide hydrogen mapping ambiguous")
            positions = [index for index, line in enumerate(rows) if line[21:27] == terminal_key]
            rows.insert(positions[-1] + 1, originals[0][:12] + "   H" + originals[0][16:21] + terminal_key + originals[0][27:])
            changes.append({"action": "RETAIN_SOURCE_AMIDE_HYDROGEN", "source_atom": hydrogen["name"], "model_atom": "H"})
        source_number = record["id"].split(":")[2]
        group = [line for line in raw_rows if line[21:22].strip() == source_chain and line[17:20].strip() == source_name and line[22:27].strip() == source_number]
        converted = ["ATOM  " + line[6:12] + f"{aliases[line[12:16].strip()]:>4}" + line[16:17] + f"{target:>3}" +
                     line[20:21] + (chain or ' ') + pdb_residue_fields(number) + line[27:] for line in group]
        indices = [index for index, line in enumerate(rows) if line.startswith("ATOM  ") and line[21:22].strip() == chain]
        insertion = indices[0] if target == "ACE" else indices[-1] + 1
        rows[insertion:insertion] = converted
        retained.append({"raw_id": record["id"], "prepared_id": f"{chain or '_'}:{target}:{number}",
                         "model": target, "n_raw_atoms_retained": len(group), "atom_aliases": aliases})
        attachments.append({"raw_id": record["id"], "partner_residue": terminal_key.strip(), "distance_angstrom": distance})
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(line[:6] + f"{index:5d}" + line[11:] for index, line in enumerate(rows, 1)) + "\nEND\n")
    pdb = PDBFile(str(output))
    forcefield = ForceField("amber14/protein.ff14SB.xml", "amber14/tip3p.xml")
    templates = forcefield.getMatchingTemplates(pdb.topology)
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff, constraints=None, rigidWater=False)
    cap_models = []
    for residue, template in zip(pdb.topology.residues(), templates, strict=True):
        if residue.name in CAP_ATOMS:
            charge = sum(atom.parameters["charge"] for atom in template.atoms)
            if template.name != residue.name or abs(charge) > 1e-6:
                raise ValueError("cap model template or neutral charge mismatch")
            cap_models.append({"chain": residue.chain.id, "residue": residue.id,
                               "insertion_code": residue.insertionCode,
                               "template": template.name, "charge": charge})
    if len(cap_models) != len(specials) or system.getNumParticles() != len(rows):
        raise ValueError("cap model coverage mismatch")
    positions = np.asarray(pdb.positions.value_in_unit(unit.angstrom))
    if not np.isfinite(positions).all():
        raise ValueError("non-finite retained cap coordinates")
    result = {"schema": "source-terminal-cap-retention-v1", "experimental_labels_used": False,
        "raw_sha256": sha(raw), "baseline_prepared_sha256": sha(prepared), "prepared_sha256": sha(output),
        "retained": retained, "terminal_state_corrections": changes, "attachments": attachments,
        "validated_cap_models": cap_models, "forcefield_system_created": True, "n_particles": system.getNumParticles(),
        "scope": "full MM topology validation; not an approved in-pocket FMO cap fragment model"}
    output.with_suffix(".model.json").write_text(json.dumps(result, indent=2) + "\n")
    return result
