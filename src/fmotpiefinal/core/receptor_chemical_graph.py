"""Packaged receptor chemistry and coordination audit; never silently edits."""
from collections import Counter, defaultdict
import hashlib
import math
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


CANONICAL = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
VARIANTS = {"ASH": "ASP", "GLH": "GLU", "HID": "HIS", "HIE": "HIS", "HIP": "HIS", "CYX": "CYS", "CYM": "CYS", "LYN": "LYS"}
WATERS = {"HOH", "WAT", "TIP"}
METALS = set("LI NA K RB CS MG CA SR BA ZN FE MN CU CO NI CD HG AL CR V MO W".split())
COORDINATION_LIMITS = {"LI": 2.8, "NA": 3.2, "K": 3.8, "MG": 2.8, "CA": 3.2, "ZN": 3., "FE": 3., "MN": 3., "CU": 3., "CO": 3., "NI": 3.}


def fingerprint(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def residue_id(residue):
    return f"{residue.chain.id.strip() or '_'}:{residue.name}:{residue.id.strip()}{residue.insertionCode.strip()}"


def atom_id(atom):
    return residue_id(atom.residue) + "/" + atom.name


def missing_skf_pairs(elements, directory):
    if directory is None:
        return None
    elements = sorted(set(elements))
    return [f"{first.title()}-{second.title()}.skf" for first in elements for second in elements
            if not (Path(directory) / f"{first.title()}-{second.title()}.skf").is_file()]


def load_structure(path):
    from openmm.app import PDBFile
    from openmm import unit
    text = Path(path).read_text()
    alternatives = [line[16:17] for line in text.splitlines() if line.startswith(("ATOM  ", "HETATM")) and line[16:17].strip()]
    if alternatives:
        raise ValueError("unresolved alternate-location records; resolve occupancy/identity before chemical graph")
    pdb = PDBFile(str(path))
    if pdb.getNumFrames() != 1:
        raise ValueError("multiple MODEL frames require explicit registered conformer selection")
    coordinates = np.asarray(pdb.positions.value_in_unit(unit.angstrom))
    atoms = list(pdb.topology.atoms())
    if len(set(atom_id(atom) for atom in atoms)) != len(atoms):
        raise ValueError("duplicate chain/residue/insertion/atom identity")
    if not np.isfinite(coordinates).all():
        raise ValueError("non-finite input coordinates")
    if any(atom.element is None for atom in atoms):
        raise ValueError("atom with unresolved element")
    return pdb, atoms, coordinates


def assign_template(residue, internal_bonds, forcefield):
    from openmm.app import Topology
    topology = Topology()
    chain = topology.addChain(residue.chain.id)
    target = topology.addResidue(residue.name, chain, id=residue.id, insertionCode=residue.insertionCode)
    mapping = {atom: topology.addAtom(atom.name, atom.element, target) for atom in residue.atoms()}
    for first, second in internal_bonds:
        topology.addBond(mapping[first], mapping[second])
    try:
        template = forcefield.getMatchingTemplates(topology, ignoreExternalBonds=True)[0]
    except Exception as error:
        return {"status": "UNMATCHED", "error": str(error)[:500], "hydrogens_or_topology_may_be_incomplete": True}
    charge = sum(atom.parameters.get("charge", 0.) for atom in template.atoms)
    if abs(charge - round(charge)) > 1e-4:
        return {"status": "UNRESOLVED_NONINTEGER_TEMPLATE_CHARGE", "template": template.name, "partial_charge_sum": charge}
    return {"status": "MATCHED", "template": template.name, "formal_charge_from_template_sum": int(round(charge)),
            "partial_charge_sum": charge, "matching_ignores_external_bonds": True}


def standard_roles(atom, template, neighbors):
    name, element = atom.name, atom.element.symbol.upper()
    base = VARIANTS.get(atom.residue.name, atom.residue.name)
    attached_hydrogens = [neighbor.name for neighbor in neighbors if neighbor.element.symbol == "H"]
    donor = bool(attached_hydrogens) and element in {"N", "O", "S"}
    acceptor = False
    uncertain = False
    if element == "O":
        acidic_oxygen = (base == "ASP" and name in {"OD1", "OD2"}) or (base == "GLU" and name in {"OE1", "OE2"})
        acceptor = not (acidic_oxygen and attached_hydrogens)
    elif element == "N":
        acceptor = (base == "HIS" and name in {"ND1", "NE2"} and not attached_hydrogens and "HIP" not in template)
        if base == "LYS" and name == "NZ":
            acceptor = "LYN" in template
    elif element == "S":
        if base == "MET" and name == "SD":
            acceptor = True
        elif base == "CYS" and name == "SG":
            acceptor = "CYX" not in template
            uncertain = "CYX" in template
    return {"donor": donor, "acceptor": acceptor, "attached_hydrogens": attached_hydrogens,
            "role_assignment": "canonical_template_and_observed_bonds", "role_uncertain": uncertain,
            "not_a_pKa_or_population_estimate": True}


def ccd_component(path, actual_atom_names):
    import gemmi
    from rdkit import Chem, RDConfig
    from rdkit.Chem import ChemicalFeatures
    block = gemmi.cif.read_file(str(path)).sole_block()
    rows = block.find(["_chem_comp_atom.atom_id", "_chem_comp_atom.type_symbol", "_chem_comp_atom.charge"])
    if not rows:
        raise ValueError("CCD atom IDs/elements/formal charges missing")
    molecule = Chem.RWMol()
    names, charges, mapping = [], {}, {}
    for row in rows:
        name, symbol, charge = [gemmi.cif.as_string(value) for value in row]
        if name in mapping or charge in {"?", "."}:
            raise ValueError("ambiguous CCD atom identity or charge")
        atom = Chem.Atom(symbol)
        atom.SetFormalCharge(int(charge))
        mapping[name] = molecule.AddAtom(atom)
        names.append(name)
        charges[name] = int(charge)
    bond_rows = block.find(["_chem_comp_bond.atom_id_1", "_chem_comp_bond.atom_id_2", "_chem_comp_bond.value_order"])
    orders = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE, "TRIP": Chem.BondType.TRIPLE, "AROM": Chem.BondType.AROMATIC}
    bonds = []
    for row in bond_rows:
        first, second, order = [gemmi.cif.as_string(value) for value in row]
        if order not in orders:
            raise ValueError("unsupported CCD bond order")
        molecule.AddBond(mapping[first], mapping[second], orders[order])
        bonds.append({"first": first, "second": second, "order": order})
    molecule = molecule.GetMol()
    Chem.SanitizeMol(molecule)
    unknown = sorted(set(actual_atom_names) - set(mapping))
    missing_heavy = sorted(name for name, index in mapping.items() if molecule.GetAtomWithIdx(index).GetAtomicNum() > 1 and name not in actual_atom_names)
    features = defaultdict(set)
    definition = Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"
    factory = ChemicalFeatures.BuildFeatureFactory(str(definition))
    for feature in factory.GetFeaturesForMol(molecule):
        for index in feature.GetAtomIds():
            features[names[index]].add(feature.GetFamily())
    return {"status": "CCD_MODEL_ANNOTATED" if not unknown and not missing_heavy else "CCD_INVENTORY_MISMATCH",
        "ccd_sha256": fingerprint(path), "feature_definition_sha256": fingerprint(definition), "unknown_input_names": unknown,
        "missing_heavy_atoms": missing_heavy, "formal_charge_of_dictionary_component": sum(charges.values()),
        "atoms": {name: {"formal_charge": charges[name], "feature_families": sorted(features[name])} for name in names},
        "bonds": bonds, "dictionary_protonation_not_verified_for_input": True,
        "forcefield_and_fmo_support_not_implied": True}


def audit(path, parameter_dir=None, ccd_cache=None):
    from openmm.app import ForceField, Modeller
    from openmm import __version__ as openmm_version
    from .special_chemistry import registered_phosphoresidue_model, load_phosphoresidue_bonds, parameter_paths
    special_model = registered_phosphoresidue_model(path)
    if special_model is not None:
        load_phosphoresidue_bonds()
    pdb, atoms, coordinates = load_structure(path)
    topology = pdb.topology
    residues = list(topology.residues())
    by_residue = {residue: list(residue.atoms()) for residue in residues}
    neighbors = {atom: set() for atom in atoms}
    internal = defaultdict(list)
    cross = []
    covalent = []
    for first, second in topology.bonds():
        covalent.append({"source": atom_id(first), "target": atom_id(second), "basis": "PDB_topology_inferred_or_explicit_not_independent_bond_validation"})
        neighbors[first].add(second)
        neighbors[second].add(first)
        if first.residue == second.residue:
            internal[first.residue].append((first, second))
        else:
            distance = float(np.linalg.norm(coordinates[first.index] - coordinates[second.index]))
            peptide = {first.name, second.name} == {"C", "N"}
            disulfide = first.element.symbol == "S" and second.element.symbol == "S"
            cross.append({"source": atom_id(first), "target": atom_id(second), "distance_angstrom": distance,
                          "kind": "peptide" if peptide else "disulfide" if disulfide else "explicit_covalent_link",
                          "geometry_supported": 1.15 <= distance <= 1.65 if peptide else 1.7 <= distance <= 2.4 if disulfide else None})
    inferred = []
    unresolved_hydrogen = []
    for atom in atoms:
        if atom.element.symbol != "H" or neighbors[atom]:
            continue
        candidates = sorted((float(np.linalg.norm(coordinates[atom.index] - coordinates[parent.index])), parent.index, parent)
                            for parent in by_residue[atom.residue] if parent.element.symbol in {"C", "N", "O", "S", "P"})
        if not candidates:
            unresolved_hydrogen.append(atom_id(atom))
            continue
        distance, _, parent = candidates[0]
        maximum = 1.55 if parent.element.symbol in {"S", "P"} else 1.3
        if distance < .6 or distance > maximum or (len(candidates) > 1 and candidates[1][0] - distance < .1):
            unresolved_hydrogen.append(atom_id(atom))
            continue
        neighbors[atom].add(parent)
        neighbors[parent].add(atom)
        internal[atom.residue].append((atom, parent))
        topology.addBond(atom, parent)
        covalent.append({"source": atom_id(atom), "target": atom_id(parent), "basis": "unique_H_parent_distance_inference"})
        inferred.append({"hydrogen": atom_id(atom), "parent": atom_id(parent), "distance_angstrom": distance,
                         "basis": "unique_same_residue_nearest_heavy_for_missing_topology_bond"})
    forcefield = ForceField("amber14/protein.ff14SB.xml", "amber14/tip3p.xml")
    from .tartrate import component_molecules
    registered_molecules, registered_models = component_molecules(path)
    registered_ids = {row["id"] for row in registered_models}
    registered_ion_charges = {}
    if special_model is not None:
        forcefield = ForceField(*map(str, parameter_paths()))
        registered_ids.update(special_model["matched_special_templates"])
        registered_ion_charges.update({identifier: template["charge"]
            for identifier, template in special_model["matched_special_templates"].items()})
    counterion_model = Path(path).with_suffix(".model.json")
    if counterion_model.exists():
        import json
        counterions = json.loads(counterion_model.read_text())
        if counterions.get("schema") in {"retained-halide-model-v1", "retained-monatomic-ion-model-v1"}:
            if counterions.get("prepared_sha256") != fingerprint(path):
                raise ValueError("registered counterion model changed")
            if any(not Path(parameter).is_file() or fingerprint(parameter) != digest
                   for parameter, digest in counterions.get("parameter_sha256", {}).items()):
                raise ValueError("registered ion forcefield parameters changed")
            from .counterion_preparation import ION_CHARGES
            for identifier, model in counterions["matched_ion_templates"].items():
                element = identifier.split(":")[1]
                if element not in ION_CHARGES or abs(model["charge"] - ION_CHARGES[element]) > 1e-8:
                    raise ValueError("registered ion identity or charge changed")
                registered_ion_charges[identifier] = model["charge"]
            registered_ids.update(counterions["matched_ion_templates"])
    try:
        ligand_residues = {residue for residue in residues if residue.name == "UNK" and residue.chain.id == "X"}
        if any(first.residue != second.residue and ((first.residue in ligand_residues) != (second.residue in ligand_residues)) for first, second in topology.bonds()):
            raise ValueError("covalent ligand requires an explicit joint topology")
        matching = Modeller(topology, pdb.positions)
        matching.delete(ligand_residues)
        if registered_molecules:
            from openmm.app import NoCutoff
            from openmmforcefields.generators import SystemGenerator
            generator = SystemGenerator(forcefields=["amber14/protein.ff14SB.xml", "amber14/tip3p.xml"],
                small_molecule_forcefield="openff-2.2.0", molecules=registered_molecules,
                nonperiodic_forcefield_kwargs={"nonbondedMethod": NoCutoff})
            generator.create_system(matching.topology)
            forcefield = generator.forcefield
        templates = forcefield.getMatchingTemplates(matching.topology)
        complete_templates = {}
        matching_residues = [residue for residue in residues if residue not in ligand_residues]
        for residue, template in zip(matching_residues, templates, strict=True):
            charge = sum(atom.parameters.get("charge", 0.) for atom in template.atoms)
            expected = next((model["formal_charge"] for model in registered_models if model["id"] == residue_id(residue)), None)
            if residue_id(residue) in registered_ion_charges:
                expected = registered_ion_charges[residue_id(residue)]
            if expected is not None and abs(charge - expected) > 1e-5:
                raise ValueError("registered component template charge disagrees with explicit source model")
            complete_templates[residue] = {"status": "MATCHED", "template": template.name,
                "formal_charge_from_template_sum": int(round(charge)), "partial_charge_sum": charge,
                "matching_ignores_external_bonds": False}
    except Exception:
        complete_templates = {}
    residue_records, atom_records = [], {}
    for residue in residues:
        group = by_residue[residue]
        base = VARIANTS.get(residue.name, residue.name)
        elements = [atom.element.symbol.upper() for atom in group]
        kind = ("ligand" if residue.name == "UNK" and residue.chain.id == "X" else "metal" if len(group) == 1 and elements[0] in METALS
                else "water" if residue.name in WATERS else "canonical_protein" if base in CANONICAL
                else "terminal_cap" if residue.name in {"ACE", "NME"} and residue in complete_templates
                else "registered_component" if residue_id(residue) in registered_ids and residue in complete_templates
                else "other_ion" if len(group) == 1 and elements[0] in {"F", "CL", "BR", "I"} else "nonstandard_component")
        record = {"residue": residue_id(residue), "kind": kind, "n_atoms": len(group), "elements": dict(Counter(elements))}
        if kind in {"canonical_protein", "water", "metal", "other_ion", "terminal_cap", "registered_component"}:
            template = complete_templates.get(residue) or assign_template(residue, internal[residue], forcefield)
            record["forcefield_template"] = template
        else:
            template = {"status": "NOT_ASSIGNED"}
        if kind == "nonstandard_component":
            component_path = Path(ccd_cache) / f"{residue.name}.cif" if ccd_cache else None
            if component_path is not None and component_path.is_file():
                try:
                    record["ccd"] = ccd_component(component_path, {atom.name for atom in group})
                except Exception as error:
                    record["ccd"] = {"status": "CCD_UNRESOLVED", "error": repr(error)}
            else:
                record["ccd"] = {"status": "CCD_REQUIRED", "component": residue.name}
        for atom in group:
            row = {"residue": residue_id(residue), "element": atom.element.symbol.upper(), "xyz": coordinates[atom.index].tolist(), "kind": kind}
            if template["status"] == "MATCHED" and kind == "canonical_protein":
                row["chemical_role"] = standard_roles(atom, template["template"], neighbors[atom])
            elif kind == "water" and template["status"] == "MATCHED" and atom.element.symbol == "O":
                row["chemical_role"] = {"donor": any(neighbor.element.symbol == "H" for neighbor in neighbors[atom]), "acceptor": True,
                                       "role_uncertain": False, "role_assignment": "matched_water_template"}
            elif kind == "nonstandard_component" and record["ccd"]["status"] == "CCD_MODEL_ANNOTATED":
                chemical = record["ccd"]["atoms"].get(atom.name)
                row["chemical_role"] = {"donor": "Donor" in chemical["feature_families"] if chemical else None,
                    "acceptor": "Acceptor" in chemical["feature_families"] if chemical else None, "role_uncertain": True,
                    "role_assignment": "CCD_dictionary_protonation_not_validated_for_input"}
            else:
                row["chemical_role"] = {"donor": None, "acceptor": None, "role_uncertain": True, "role_assignment": "unresolved_or_separate_ligand_model"}
            atom_records[atom_id(atom)] = row
        residue_records.append(record)
    coordinate_edges = []
    for atom in atoms:
        element = atom.element.symbol.upper()
        if element not in METALS:
            continue
        cutoff = COORDINATION_LIMITS.get(element, 3.)
        candidates = []
        for other in atoms:
            if other.residue == atom.residue or other.element.symbol not in {"N", "O", "S"}:
                continue
            distance = float(np.linalg.norm(coordinates[atom.index] - coordinates[other.index]))
            if distance <= cutoff:
                candidates.append({"donor_candidate": atom_id(other), "distance_angstrom": distance,
                                   "chemical_acceptor_role": atom_records[atom_id(other)]["chemical_role"]["acceptor"]})
        angles = []
        for first_index, first in enumerate(candidates):
            for second in candidates[first_index + 1:]:
                vector1 = np.asarray(atom_records[first["donor_candidate"]]["xyz"]) - coordinates[atom.index]
                vector2 = np.asarray(atom_records[second["donor_candidate"]]["xyz"]) - coordinates[atom.index]
                length = np.linalg.norm(vector1) * np.linalg.norm(vector2)
                angle = math.degrees(math.acos(float(np.clip(np.dot(vector1, vector2) / length, -1., 1.)))) if length > 1e-8 else None
                angles.append({"first": first["donor_candidate"], "second": second["donor_candidate"], "angle_degrees": angle})
        coordinate_edges.append({"metal": atom_id(atom), "element": element, "cutoff_angstrom": cutoff,
            "registered_nonbonded_mm_model": residue_id(atom.residue) in registered_ids and
                atom.residue in complete_templates and element in {"NA", "CA"},
            "n_distance_candidates": len(candidates), "distance_candidates": candidates, "angles": angles,
            "oxidation_state_assigned": False, "true_coordination_number_assigned": False,
            "production_HETATM_metal_fragment_supported": False,
            "note": "distance candidates are not confirmed coordination bonds or validated MD/QM parameters"})
    categories = Counter(record["kind"] for record in residue_records)
    geometry_screens = []
    limits = {("C", "C"): (1.1, 1.75), ("C", "N"): (1.1, 1.65), ("C", "O"): (1.05, 1.65),
              ("C", "S"): (1.55, 2.05), ("N", "N"): (1.1, 1.65), ("S", "S"): (1.7, 2.4)}
    for edge in covalent:
        first, second = atom_records[edge["source"]], atom_records[edge["target"]]
        distance = float(np.linalg.norm(np.asarray(first["xyz"]) - second["xyz"]))
        edge["distance_angstrom"] = distance
        pair = tuple(sorted((first["element"], second["element"])))
        if first["kind"] != "canonical_protein" or second["kind"] != "canonical_protein" or pair not in limits:
            continue
        minimum, maximum = limits[pair]
        if not minimum <= distance <= maximum:
            geometry_screens.append({**edge, "screen_bounds_angstrom": [minimum, maximum],
                "severe": distance < .85 or distance > 2.2,
                "basis": "fixed_broad_element_pair_screen_not_bond_order_specific_restraint_or_statistical_significance"})
    templates = Counter(record["forcefield_template"].get("template", record["forcefield_template"]["status"])
                        for record in residue_records if "forcefield_template" in record)
    import openmm.app
    definitions = [Path(openmm.app.__file__).parent / "data" / name for name in ("amber14/protein.ff14SB.xml", "amber14/tip3p.xml")]
    from .water_hydrogen_clash import hydrogen_clashes
    water_h_contacts = hydrogen_clashes(path)
    return {"schema": "receptor-chemical-graph-v1", "source_path": str(Path(path).resolve()), "source_sha256": fingerprint(path),
        "experimental_labels_used": False, "coordinates_or_scores_modified": False, "openmm_version": openmm_version,
        "forcefield_definition_sha256": {str(path): fingerprint(path) for path in definitions},
        "residues": residue_records, "atoms": atom_records, "covalent_edges": covalent,
        "cross_residue_bonds": cross, "inferred_H_bonds": inferred,
        "protein_bond_geometry_screens": geometry_screens,
        "water_hydrogen_clashes": water_h_contacts,
        "unresolved_hydrogen_parent_atoms": unresolved_hydrogen, "metal_coordination": coordinate_edges,
        "elements_without_all_ordered_SKF_pairs": missing_skf_pairs({atom.element.symbol.upper() for atom in atoms}, parameter_dir),
        "summary": {"n_atoms": len(atoms), "n_residues": len(residues), "residue_kind_counts": dict(categories),
                    "template_counts": dict(templates), "n_unmatched_canonical_residues": sum(record["kind"] == "canonical_protein" and record["forcefield_template"]["status"] != "MATCHED" for record in residue_records),
                    "n_protein_bond_geometry_screen_outliers": len(geometry_screens),
                    "n_water_hydrogen_clashes": len(water_h_contacts),
                    "n_invalid_peptide_or_disulfide_geometry": sum(edge["geometry_supported"] is False for edge in cross)}}


def raw_component_retention(raw_path, prepared_graph):
    prepared = prepared_graph["atoms"]
    by_element = defaultdict(list)
    for identifier, row in prepared.items():
        if row["element"] != "H":
            by_element[row["element"]].append((identifier, row["xyz"]))
    trees = {element: cKDTree([row[1] for row in values]) for element, values in by_element.items()}
    raw_residues = defaultdict(list)
    alternatives = []
    model = "1"
    for line in Path(raw_path).read_text().splitlines():
        if line.startswith("MODEL "):
            model = line[10:14].strip()
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        name, residue, chain = line[12:16].strip(), line[17:20].strip(), line[21:22].strip() or "_"
        number, alternative = line[22:27].strip(), line[16:17].strip()
        element = line[76:78].strip().upper()
        if not element:
            raise ValueError("raw inventory requires explicit element columns; do not infer metals from atom names")
        row = {"id": f"{chain}:{residue}:{number}/{name}", "element": element,
            "xyz": [float(line[start:start + 8]) for start in (30, 38, 46)], "altloc": alternative,
            "occupancy": float(line[54:60]) if line[54:60].strip() else None}
        if not np.isfinite(row["xyz"]).all():
            raise ValueError("raw inventory non-finite coordinates")
        if alternative:
            alternatives.append({"atom": row["id"], "model": model, "altloc": alternative, "occupancy": row["occupancy"]})
        if VARIANTS.get(residue, residue) in CANONICAL or residue in WATERS or (residue == "UNK" and chain == "X"):
            continue
        raw_residues[(model, f"{chain}:{residue}:{number}")].append(row)
    records = []
    for (model, identifier), group in raw_residues.items():
        matches = []
        for atom in group:
            element = atom["element"]
            if element == "H":
                continue
            distance, nearest = trees[element].query(atom["xyz"]) if element in trees else (float("inf"), None)
            matches.append({"raw_atom": atom["id"], "element": element, "altloc": atom["altloc"],
                            "prepared_atom_within_0.1A": by_element[element][int(nearest)][0] if distance <= .1 else None,
                            "nearest_distance_angstrom": float(distance) if np.isfinite(distance) else None})
        records.append({"raw_residue": identifier, "model": model, "n_heavy_atoms": len(matches), "n_matched": sum(row["prepared_atom_within_0.1A"] is not None for row in matches),
            "contains_metal": any(atom["element"] in METALS for atom in group), "atoms": matches,
            "identity_conclusion": "coordinate_inventory_diagnostic_only; requires provenance for intentional removal, mutation or changed frame"})
    return {"raw_source_sha256": fingerprint(raw_path), "prepared_source_sha256": prepared_graph["source_sha256"],
            "residue_numbering_independent_coordinate_match": True, "global_alignment_performed": False,
            "alternate_location_records": alternatives, "alternate_records_not_silently_selected": True, "components": records}
