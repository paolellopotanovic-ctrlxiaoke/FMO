"""Trace size trends to frozen structures and archived ligand-fragment PIEDA."""
from __future__ import annotations

import argparse
import gzip
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator, rdFMCS
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist

from build_model_comparison import digest, read_csv
from diagnose_molecular_weight import pearson


COMPONENTS = {"tpie": "value_kcal_mol", "ES": "electrostatic_kcal_mol",
              "EX": "exchange_repulsion_kcal_mol", "CT": "charge_transfer_kcal_mol",
              "DI": "dispersion_kcal_mol", "SOLV": "solvation_kcal_mol"}
TARGETS = ("OpenBind", "P38", "CHK1", "BACE", "CDK2", "TYK2", "FAAH")


def parse_deck(text):
    section = re.search(r"\$FMO\s+(.*?)\$END", text, re.S | re.I).group(1)
    names = re.search(r"FRGNAM\(1\)=(.*?)INDAT", section, re.S).group(1)
    names = [name.strip() for name in names.replace("\n", "").split(",") if name.strip()]
    values = [int(value) for value in re.search(r"INDAT\(1\)=(.*)", section, re.S).group(1).split()]
    groups = []
    current = []
    for value in values:
        if value == 0:
            if current:
                groups.append(current)
                current = []
        elif value > 0:
            current.append(value - 1)
        else:
            if not current or -value <= current[-1] + 1:
                raise ValueError("invalid INDAT range")
            current.extend(range(current[-1] + 1, -value))
    if current:
        raise ValueError("unterminated INDAT")
    xyz_section = re.search(r"\$FMOXYZ\s*\n(.*?)\$END", text, re.S | re.I).group(1)
    entries = [line.split() for line in xyz_section.splitlines() if line.strip()]
    atomic_numbers = np.asarray([int(float(entry[1])) for entry in entries])
    coordinates = np.asarray([[float(value) for value in entry[2:5]] for entry in entries])
    if len(groups) != len(names) or sorted(index for group in groups for index in group) != list(range(len(entries))):
        raise ValueError("fragment atom mapping is incomplete")
    return names, groups, atomic_numbers, coordinates


def ligand_surface(ligand_xyz, ligand_numbers, receptor_xyz, receptor_numbers, points=192):
    table = Chem.GetPeriodicTable()
    radii = np.asarray([table.GetRvdw(int(number)) + 1.4 for number in ligand_numbers])
    receptor_radii = np.asarray([table.GetRvdw(int(number)) + 1.4 for number in receptor_numbers])
    position = np.arange(points) + .5
    vertical = 1 - 2 * position / points
    radial = np.sqrt(1 - vertical ** 2)
    angle = position * math.pi * (3 - math.sqrt(5))
    sphere = np.column_stack([radial * np.cos(angle), radial * np.sin(angle), vertical])
    receptor_tree = cKDTree(receptor_xyz)
    free, bound = [], []
    for index, (coordinate, radius) in enumerate(zip(ligand_xyz, radii)):
        mesh = coordinate + radius * sphere
        distances = cdist(mesh, ligand_xyz)
        distances[:, index] = np.inf
        available = np.all(distances >= radii - 1e-8, axis=1)
        neighbors = receptor_tree.query_ball_point(coordinate, radius + float(receptor_radii.max()))
        occluded = np.any(cdist(mesh, receptor_xyz[neighbors]) < receptor_radii[neighbors] - 1e-8, axis=1) if neighbors else np.zeros(points, dtype=bool)
        area = 4 * math.pi * radius ** 2 / points
        free.append(float(available.sum() * area))
        bound.append(float((available & ~occluded).sum() * area))
    return np.asarray(free), np.asarray(bound)


def inspect_case(root, ligand, expected, sources, geometry_root=None):
    geometry_root = geometry_root or root
    directory = re.sub(r"[^A-Za-z0-9_.-]", "_", ligand)
    result_path = root / "qm" / directory / "run/result.json"
    deck_path = root / "qm" / directory / "input.inp"
    pdb_path = geometry_root / "minimized" / ligand / "minimized_complex.pdb"
    minimization_path = geometry_root / "minimized" / ligand / "minimization.json"
    for source in (result_path, deck_path, pdb_path, minimization_path):
        sources[str(source)] = digest(source)
    result = json.loads(result_path.read_text())
    if result["status"] != "SUCCEEDED" or abs(result["tpie_kcal_mol"] + expected) > 1e-7:
        raise ValueError(f"score identity mismatch: {result_path}")
    names, groups, numbers, coordinates = parse_deck(deck_path.read_text())
    ligand_fragment = result["ligand_fragment"]
    ligand_indices = np.asarray(groups[ligand_fragment - 1])
    ligand_indices = ligand_indices[numbers[ligand_indices] > 1]
    ligand_xyz, ligand_numbers = coordinates[ligand_indices], numbers[ligand_indices]
    pairs = result["pairs"]
    if len(pairs) != len(groups) - 1:
        raise ValueError("incomplete pair coverage")
    fragments = []
    seen = set()
    for pair in pairs:
        if ligand_fragment not in (pair["fragment_i"], pair["fragment_j"]):
            raise ValueError("pair is not a ligand interaction")
        other = pair["fragment_j"] if pair["fragment_i"] == ligand_fragment else pair["fragment_i"]
        if other == ligand_fragment or other in seen:
            raise ValueError("duplicate or self pair")
        seen.add(other)
        indices = np.asarray(groups[other - 1])
        indices = indices[numbers[indices] > 1]
        distance = float(cdist(coordinates[indices], ligand_xyz).min())
        fragments.append({"fragment": names[other - 1], "distance": distance,
                          **{key: float(pair[field]) for key, field in COMPONENTS.items()}})
    components = {key: sum(row[key] for row in fragments) for key in COMPONENTS}
    if abs(components["tpie"] - result["tpie_kcal_mol"]) > 1e-7:
        raise ValueError("PIEDA sum mismatch")
    if abs(components["tpie"] - sum(components[key] for key in COMPONENTS if key != "tpie")) > .02 * len(pairs):
        raise ValueError("PIEDA components do not close")
    pdb_atoms = []
    table = Chem.GetPeriodicTable()
    for line in pdb_path.read_text().splitlines():
        if line[:6].strip() not in ("ATOM", "HETATM"):
            continue
        element = line[76:78].strip().title()
        if element in ("H", "D"):
            continue
        pdb_atoms.append({"number": table.GetAtomicNumber(element), "xyz": tuple(float(line[start:start + 8]) for start in (30, 38, 46)),
                          "ligand": line[21] == "X" and line[17:20].strip() == "UNK",
                          "water": line[17:20].strip() in ("HOH", "WAT", "TIP"),
                          "residue": f"{line[21]}:{line[17:20].strip()}:{line[22:27].strip()}",
                          "protein": line[:6].strip() == "ATOM"})
    ligand_pdb = [atom for atom in pdb_atoms if atom["ligand"]]
    if [atom["number"] for atom in ligand_pdb] != ligand_numbers.tolist() or np.max(np.abs(np.asarray([atom["xyz"] for atom in ligand_pdb]) - ligand_xyz)) > .002:
        raise ValueError("minimized ligand does not match QM deck")
    receptor = [atom for atom in pdb_atoms if not atom["ligand"] and not atom["water"]]
    receptor_xyz = np.asarray([atom["xyz"] for atom in receptor])
    receptor_numbers = np.asarray([atom["number"] for atom in receptor])
    distances = cdist(ligand_xyz, receptor_xyz)
    nearest = distances.min(axis=1)
    contact_residues = sorted({atom["residue"] for atom, close in zip(receptor, np.any(distances <= 4, axis=0)) if close})
    free, bound = ligand_surface(ligand_xyz, ligand_numbers, receptor_xyz, receptor_numbers)
    minimization = json.loads(minimization_path.read_text())
    polar = np.isin(ligand_numbers, [7, 8, 16])
    polar_contact = polar[:, None] & np.isin(receptor_numbers, [7, 8, 16])[None, :] & (distances <= 3.5)
    water_atoms = [atom for atom in pdb_atoms if atom["water"]]
    water_close = int(np.sum(np.min(cdist(ligand_xyz, [atom["xyz"] for atom in water_atoms]), axis=0) <= 3.5)) if water_atoms else 0
    bins = {}
    for label, lower, upper in (("le4", 0, 4), ("4to6", 4, 6), ("gt6", 6, np.inf)):
        selected = [row for row in fragments if lower < row["distance"] <= upper]
        bins[label] = {key: sum(row[key] for row in selected) for key in COMPONENTS}
    return {"components": components, "fragment_rows": fragments, "distance_bins": bins,
        "n_ligand_heavy": len(ligand_xyz), "n_contact_atoms_4A": int(np.sum(nearest <= 4)),
        "n_noncontact_atoms_4A": int(np.sum(nearest > 4)), "n_contact_residues_4A": len(contact_residues),
        "heavy_contact_pairs_4A": int(np.sum(distances <= 4)), "min_heavy_distance": float(distances.min()),
        "polar_proximity_pairs_3_5A_not_HBonds": int(polar_contact.sum()), "water_oxygen_contacts_3_5A": water_close,
        "ligand_free_SASA_A2": float(free.sum()), "ligand_bound_SASA_A2": float(bound.sum()),
        "ligand_buried_SASA_A2": float((free - bound).sum()),
        "buried_fraction": float((free - bound).sum() / free.sum()),
        "k10_rmsd_all_atoms_A": minimization["ligand_rmsd_angstrom"],
        "ligand_xyz": ligand_xyz.tolist(), "ligand_atomic_numbers": ligand_numbers.tolist(),
        "atom_nearest_nonwater_distance": nearest.tolist(), "atom_buried_SASA_A2": (free - bound).tolist(),
        "source_ligand_sdf": minimization["ligand_sdf"], "ligand_name": ligand,
        "source_root": str(root)}


def chemical_molecule(record, sources, cache):
    path = Path(record["source_ligand_sdf"])
    if path not in cache:
        sources[str(path)] = digest(path)
        cache[path] = {molecule.GetProp("_Name").strip(): Chem.RemoveHs(molecule)
                       for molecule in Chem.SDMolSupplier(str(path), removeHs=False) if molecule is not None}
    molecule = cache[path][record["ligand_name"]]
    if [atom.GetAtomicNum() for atom in molecule.GetAtoms()] != record["ligand_atomic_numbers"]:
        raise ValueError("ligand source/deck heavy atom order differs")
    return molecule


def compare_pair(first, second, records, molecules):
    left, right = records[first["id"]], records[second["id"]]
    residues_left = {row["fragment"]: row for row in left["fragment_rows"]}
    residues_right = {row["fragment"]: row for row in right["fragment_rows"]}
    residuals = []
    for name in residues_left.keys() | residues_right.keys():
        values = {key: residues_right.get(name, {}).get(key, 0) - residues_left.get(name, {}).get(key, 0) for key in COMPONENTS}
        residuals.append({"fragment": name, "shared_fragment_name": name in residues_left and name in residues_right,
            "first_min_distance": residues_left.get(name, {}).get("distance"),
            "second_min_distance": residues_right.get(name, {}).get("distance"), **values})
    residuals.sort(key=lambda row: abs(row["tpie"]), reverse=True)
    all_residue_delta = sum(row["tpie"] for row in residuals)
    delta = {key: right["components"][key] - left["components"][key] for key in COMPONENTS}
    if abs(all_residue_delta - delta["tpie"]) > 1e-6:
        raise ValueError("pair residue decomposition does not close")
    return {"first": first["id"], "second": second["id"],
        "delta_MW": second["mw"] - first["mw"], "delta_label": second["label"] - first["label"],
        "delta_raw_favorable": second["raw_score"] - first["raw_score"],
        "delta_components_energy_second_minus_first": delta,
        "delta_geometry": {key: right[key] - left[key] for key in ("n_ligand_heavy", "n_contact_atoms_4A", "n_noncontact_atoms_4A",
            "heavy_contact_pairs_4A", "n_contact_residues_4A", "ligand_buried_SASA_A2", "ligand_bound_SASA_A2", "buried_fraction",
            "polar_proximity_pairs_3_5A_not_HBonds", "water_oxygen_contacts_3_5A", "k10_rmsd_all_atoms_A")},
        "top_fragment_deltas": residuals[:8],
        "sum_all_fragment_delta": all_residue_delta,
        "chemical_identity": [Chem.MolToSmiles(molecules[first["id"]]), Chem.MolToSmiles(molecules[second["id"]])],
        "same_Murcko_scaffold": first["scaffold"] == second["scaffold"]}


def polar_contacts(record):
    root = Path(record["source_root"])
    geometry = root if (root / "minimized").exists() else root.parents[1] / "qm/full_target_fmo"
    path = geometry / "minimized" / record["ligand_name"] / "minimized_complex.pdb"
    atoms = []
    for line in path.read_text().splitlines():
        if line[:6].strip() not in ("ATOM", "HETATM") or line[76:78].strip().upper() not in ("N", "O", "S"):
            continue
        atoms.append({"id": f"{line[21]}:{line[17:20].strip()}:{line[22:27].strip()}:{line[12:16].strip()}",
                      "ligand": line[21] == "X" and line[17:20].strip() == "UNK",
                      "xyz": [float(line[start:start + 8]) for start in (30, 38, 46)]})
    contacts = []
    for left in [atom for atom in atoms if atom["ligand"]]:
        for right in [atom for atom in atoms if not atom["ligand"]]:
            distance = float(np.linalg.norm(np.asarray(left["xyz"]) - right["xyz"]))
            if distance <= 3.8:
                contacts.append({"ligand_atom": left["id"], "receptor_atom": right["id"], "distance_A": distance})
    return sorted(contacts, key=lambda value: value["distance_A"])


def build(workspace, root):
    diagnostic_path = root / "docs/MOLECULAR_WEIGHT_DIAGNOSTICS.json"
    diagnostic = json.loads(diagnostic_path.read_text())
    run = workspace / "FMOtpiefinal_HIV1lab/runs/all_accelerated_robust_20261020"
    paired = read_csv(workspace / "runs/openbind/score_distribution_audit/paired_compounds.csv")
    sources_paired = workspace / "runs/openbind/score_distribution_audit/paired_compounds.csv"
    crystal_ids = {row["compound_index"]: [name.strip() for name in row["fragalysis_codes"].split(";")] for row in paired}
    sources = {str(diagnostic_path): digest(diagnostic_path)}
    sources[str(sources_paired)] = digest(sources_paired)
    groups = {}
    molecule_cache = {}
    fingerprint = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    for target in TARGETS:
        rows = [row for row in diagnostic["cases"] if row["cohort"] == target]
        records, molecules, fingerprints = {}, {}, {}
        for index, row in enumerate(rows):
            if target == "OpenBind":
                identifiers = crystal_ids[row["id"]]
                crystal_records = []
                for name in identifiers:
                    case_run = run / "OpenBind/cases" / name
                    verification = json.loads((case_run / "independent_verification.json").read_text())
                    sources[str(case_run / "independent_verification.json")] = digest(case_run / "independent_verification.json")
                    geometry_root = case_run / "run/qm/full_target_fmo"
                    selected_root = geometry_root if verification["selected_arm"] == "strict_default" else case_run / "run/rescues" / verification["selected_arm"]
                    score = verification["records"][0]["raw_score"]
                    inspected = inspect_case(selected_root, name, score, sources, geometry_root)
                    inspected["selected_arm"] = verification["selected_arm"]
                    crystal_records.append(inspected)
                record = dict(crystal_records[0])
                record["n_crystals"] = len(crystal_records)
                record["crystals"] = crystal_records
                for key in ("components", "distance_bins"):
                    if key == "components":
                        record[key] = {field: float(np.mean([entry[key][field] for entry in crystal_records])) for field in COMPONENTS}
                    else:
                        record[key] = {band: {field: float(np.mean([entry[key][band][field] for entry in crystal_records])) for field in COMPONENTS} for band in record[key]}
                for key in ("n_ligand_heavy", "n_contact_atoms_4A", "n_noncontact_atoms_4A", "n_contact_residues_4A", "heavy_contact_pairs_4A",
                            "polar_proximity_pairs_3_5A_not_HBonds", "water_oxygen_contacts_3_5A", "ligand_buried_SASA_A2", "ligand_bound_SASA_A2", "buried_fraction", "k10_rmsd_all_atoms_A"):
                    record[key] = float(np.mean([entry[key] for entry in crystal_records]))
            else:
                geometry_root = run / target / "run/qm/full_target_fmo"
                verification_path = run / target / "independent_verification.json"
                sources[str(verification_path)] = digest(verification_path)
                selected_arm = json.loads(verification_path.read_text())["selected_arm"]
                selected_root = geometry_root if selected_arm == "strict_default" else run / target / "run/rescues" / selected_arm
                record = inspect_case(selected_root, row["id"], row["raw_score"], sources, geometry_root)
                record["selected_arm"] = selected_arm
                record["n_crystals"] = 1
            if abs(record["components"]["tpie"] + row["raw_score"]) > 1e-7:
                raise ValueError("compound aggregation mismatch")
            records[row["id"]] = record
            molecule = chemical_molecule(record, sources, molecule_cache)
            record["prepared_ligand_formal_charge"] = Chem.GetFormalCharge(molecule)
            record["benchmark_descriptor_formal_charge"] = row["formal_charge"]
            molecules[row["id"]] = molecule
            fingerprints[row["id"]] = fingerprint.GetFingerprint(molecule)
            if index % 50 == 0:
                print(target, index + 1, "/", len(rows), flush=True)
        summary = {}
        labels = np.asarray([row["label"] for row in rows])
        masses = np.asarray([row["mw"] for row in rows])
        for field in ("n_contact_atoms_4A", "n_noncontact_atoms_4A", "heavy_contact_pairs_4A", "ligand_buried_SASA_A2", "ligand_bound_SASA_A2", "buried_fraction", "k10_rmsd_all_atoms_A"):
            values = np.asarray([records[row["id"]][field] for row in rows])
            summary[field] = {"label_pearson": pearson(labels, values), "MW_pearson": pearson(masses, values), "mean": float(values.mean())}
        for field in COMPONENTS:
            values = -np.asarray([records[row["id"]]["components"][field] for row in rows])
            summary["favorable_" + field] = {"label_pearson": pearson(labels, values), "MW_pearson": pearson(masses, values)}
        for band in ("le4", "4to6", "gt6"):
            values = -np.asarray([records[row["id"]]["distance_bins"][band]["tpie"] for row in rows])
            summary["favorable_distance_" + band] = {"label_pearson": pearson(labels, values), "MW_pearson": pearson(masses, values)}
        charge_changes = [{"id": row["id"], "benchmark_charge": row["formal_charge"],
            "prepared_charge": records[row["id"]]["prepared_ligand_formal_charge"]}
            for row in rows if row["formal_charge"] != records[row["id"]]["prepared_ligand_formal_charge"]]
        candidates = []
        for first_index, first in enumerate(rows):
            for second in rows[first_index + 1:]:
                if abs(first["mw"] - second["mw"]) > 10 or abs(first["label"] - second["label"]) < .5:
                    continue
                similarity = float(DataStructs.TanimotoSimilarity(fingerprints[first["id"]], fingerprints[second["id"]]))
                if similarity < .65:
                    continue
                correct = (first["label"] - second["label"]) * (first["raw_score"] - second["raw_score"]) > 0
                candidates.append({"first": first, "second": second, "similarity": similarity, "correct": correct})
        candidates.sort(key=lambda candidate: (-abs(candidate["first"]["label"] - candidate["second"]["label"]), candidate["first"]["id"], candidate["second"]["id"]))
        selected = []
        for correct in (False, True):
            eligible = [candidate for candidate in candidates if candidate["correct"] == correct
                        and records[candidate["first"]["id"]]["n_crystals"] == records[candidate["second"]["id"]]["n_crystals"] == 1]
            for candidate in eligible[:3]:
                first, second = candidate["first"], candidate["second"]
                comparison = compare_pair(first, second, records, molecules)
                comparison.update({"correct": correct, "Tanimoto": candidate["similarity"]})
                comparison["polar_contacts_not_validated_HBonds"] = {item["id"]: polar_contacts(records[item["id"]]) for item in (first, second)}
                mcs = rdFMCS.FindMCS([molecules[first["id"]], molecules[second["id"]]], timeout=3,
                    ringMatchesRingOnly=True, completeRingsOnly=True, bondCompare=rdFMCS.BondCompare.CompareOrderExact)
                comparison["MCS_canceled"] = bool(mcs.canceled)
                if not mcs.canceled:
                    query = Chem.MolFromSmarts(mcs.smartsString)
                    comparison["MCS_heavy_atoms"] = mcs.numAtoms
                    added = {}
                    for item in (first, second):
                        molecule = molecules[item["id"]]
                        matches = molecule.GetSubstructMatches(query, maxMatches=100)
                        match = matches[0]
                        noncore = sorted(set(range(molecule.GetNumAtoms())) - set(match))
                        record = records[item["id"]]
                        added[item["id"]] = {"noncore_indices_zero_based": noncore, "n_MCS_matches_up_to_100": len(matches),
                            "noncore_elements": [molecule.GetAtomWithIdx(index).GetSymbol() for index in noncore],
                            "noncore_buried_SASA_A2": sum(record["atom_buried_SASA_A2"][index] for index in noncore),
                            "noncore_contact_atoms_4A": sum(record["atom_nearest_nonwater_distance"][index] <= 4 for index in noncore)}
                    comparison["MCS_noncore"] = added
                selected.append(comparison)
        groups[target] = {"n": len(rows), "correlations": summary, "charge_representation_differences": charge_changes, "records": records,
            "eligible_pairs": [{"first": candidate["first"]["id"], "second": candidate["second"]["id"], "correct": candidate["correct"], "Tanimoto": candidate["similarity"]} for candidate in candidates],
            "representative_pairs": selected}
        print(target, "completed", len(candidates), "eligible pairs", flush=True)
    result = {"schema": "size-structure-interaction-diagnostic-v1", "date": "2026-10-21",
        "production_changed": False, "new_QM": False,
        "selection": "All seven requested diagnostic/control cohorts; all compounds/crystals; pairs within 10 Da, label gap >=0.5 in within-target label units and Morgan Tanimoto >=0.65. Representatives: top three largest-gap single-crystal pairs in each correct/incorrect category, no manual picking.",
        "geometry": "Frozen minimized ligand checked against exact QM deck; nonwater full receptor heavy atoms for contact/SASA; 192-point sphere, RDKit vdW radii +1.4 A; water excluded for receptor burial, reported separately. Polar proximity is NOT a validated H-bond. Fragment boundaries include caps/backbone reassignment; fragment name is not a pure residue atomic energy.",
        "limitations": "Posthoc descriptive selection uses labels, not a deployable score. No atomwise QM energy attribution. MCS symmetry chooses first match and records ambiguity. Residue comparison across OpenBind structures is nominal fragment-name correspondence, not verified structural equivalence. No entropy, isolated ligand strain energy, or binding desolvation computed.",
        "groups": groups, "sources_sha256": sources}
    evidence_path = workspace / "FMOtpiefinal_HIV1lab/runs/size_structure_interactions_20261021/evidence.json.gz"
    evidence_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(evidence_path, "wt", encoding="utf-8") as handle:
        json.dump(result, handle, ensure_ascii=False)
    result["full_evidence"] = {"path": str(evidence_path), "sha256": digest(evidence_path)}
    excluded = {"fragment_rows", "crystals", "ligand_xyz", "ligand_atomic_numbers", "atom_nearest_nonwater_distance", "atom_buried_SASA_A2"}
    for group in result["groups"].values():
        group["records"] = {identifier: {key: value for key, value in record.items() if key not in excluded}
                            for identifier, record in group["records"].items()}
    (root / "docs/SIZE_STRUCTURE_INTERACTIONS.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", required=True, type=Path)
    args = parser.parse_args()
    build(args.workspace.resolve(), Path(__file__).resolve().parents[1])
