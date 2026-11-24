"""Validate local chemistry and geometric explanations of archived misrankings."""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from rdkit import Chem
from scipy.spatial.distance import cdist

from build_model_comparison import digest, metrics, read_csv
from diagnose_molecular_weight import pearson


def atoms(path):
    result = []
    for line in path.read_text().splitlines():
        if line[:6].strip() in ("ATOM", "HETATM"):
            result.append({"chain": line[21], "residue": line[17:20].strip(), "number": line[22:27].strip(),
                "name": line[12:16].strip(), "element": line[76:78].strip(),
                "xyz": np.asarray([float(line[start:start + 8]) for start in (30, 38, 46)])})
    return result


def atom_key(atom):
    return atom["chain"], atom["residue"], atom["number"], atom["name"]


def kabsch(first, second):
    center_first, center_second = first.mean(axis=0), second.mean(axis=0)
    left, _, right = np.linalg.svd((first - center_first).T @ (second - center_second))
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(left @ right)
    rotation = left @ correction @ right
    translation = center_second - center_first @ rotation
    return rotation, translation


def angle(donor, hydrogen, acceptor):
    first, second = donor - hydrogen, acceptor - hydrogen
    return float(np.degrees(np.arccos(np.clip(np.dot(first, second) / np.linalg.norm(first) / np.linalg.norm(second), -1, 1))))


def local_hydrogen_bonds(record, sources):
    root = Path(record["source_root"])
    pdb_path = root / "minimized" / record["ligand_name"] / "minimized_complex.pdb"
    sdf_path = Path(record["source_ligand_sdf"])
    sources[str(pdb_path)] = digest(pdb_path)
    sources[str(sdf_path)] = digest(sdf_path)
    parsed = atoms(pdb_path)
    ligand = [atom for atom in parsed if atom["chain"] == "X" and atom["residue"] == "UNK"]
    receptor = [atom for atom in parsed if not (atom["chain"] == "X" and atom["residue"] == "UNK")]
    molecule = next(molecule for molecule in Chem.SDMolSupplier(str(sdf_path), removeHs=False)
                    if molecule.GetProp("_Name").strip() == record["ligand_name"])
    if len(ligand) != molecule.GetNumAtoms() or [atom["element"].upper() for atom in ligand] != [atom.GetSymbol().upper() for atom in molecule.GetAtoms()]:
        raise ValueError("all-atom identity mismatch")
    contacts = []
    for donor in molecule.GetAtoms():
        if donor.GetSymbol() not in ("N", "O", "S"):
            continue
        for hydrogen in donor.GetNeighbors():
            if hydrogen.GetAtomicNum() != 1:
                continue
            donor_atom, hydrogen_atom = ligand[donor.GetIdx()], ligand[hydrogen.GetIdx()]
            for acceptor in receptor:
                if acceptor["element"] not in ("N", "O", "S"):
                    continue
                distance = float(np.linalg.norm(donor_atom["xyz"] - acceptor["xyz"]))
                if distance > 3.5:
                    continue
                contacts.append({"donor": donor_atom["name"], "hydrogen": hydrogen_atom["name"],
                    "acceptor": ":".join(atom_key(acceptor)), "DA_A": distance,
                    "pre_k10_DA_A": float(np.linalg.norm(np.asarray(molecule.GetConformer().GetAtomPosition(donor.GetIdx())) - acceptor["xyz"])),
                    "HA_A": float(np.linalg.norm(hydrogen_atom["xyz"] - acceptor["xyz"])),
                    "DHA_degrees": angle(donor_atom["xyz"], hydrogen_atom["xyz"], acceptor["xyz"])})
    reverse = []
    for donor in receptor:
        if donor["element"] not in ("N", "O", "S"):
            continue
        hydrogens = [atom for atom in receptor if atom["element"] == "H" and
            atom_key(atom)[:3] == atom_key(donor)[:3] and np.linalg.norm(atom["xyz"] - donor["xyz"]) <= 1.3]
        for hydrogen in hydrogens:
            for acceptor in ligand:
                if acceptor["element"] not in ("N", "O", "S"):
                    continue
                distance = float(np.linalg.norm(donor["xyz"] - acceptor["xyz"]))
                if distance <= 3.5:
                    reverse.append({"donor": ":".join(atom_key(donor)), "hydrogen": hydrogen["name"],
                        "acceptor": acceptor["name"], "DA_A": distance,
                        "HA_A": float(np.linalg.norm(hydrogen["xyz"] - acceptor["xyz"])),
                        "DHA_degrees": angle(donor["xyz"], hydrogen["xyz"], acceptor["xyz"])})
    return {"ligand_donor_contacts": contacts, "receptor_donor_proximity_not_chemical_acceptor_validation": reverse,
        "ASP_states": [{"chain": chain, "number": number,
            "sidechain_H_names": [atom["name"] for atom in parsed if atom["chain"] == chain and atom["residue"] == "ASP" and atom["number"] == number and atom["name"].startswith("HD")]}
            for chain, number in sorted({(atom["chain"], atom["number"]) for atom in receptor if atom["residue"] == "ASP"})]}


def align_openbind(first, second, sources):
    paths = [Path(record["source_root"]) / "minimized" / record["ligand_name"] / "minimized_complex.pdb" for record in (first, second)]
    for path in paths:
        sources[str(path)] = digest(path)
    parsed_first, parsed_second = map(atoms, paths)
    first_ca = {atom_key(atom): atom for atom in parsed_first if atom["name"] == "CA" and atom["residue"] != "UNK"}
    second_ca = {atom_key(atom): atom for atom in parsed_second if atom["name"] == "CA" and atom["residue"] != "UNK"}
    keys = sorted(first_ca.keys() & second_ca.keys())
    if len(keys) < 30:
        raise ValueError("too few matched C-alpha atoms")
    first_xyz = np.asarray([first_ca[key]["xyz"] for key in keys])
    second_xyz = np.asarray([second_ca[key]["xyz"] for key in keys])
    rotation, translation = kabsch(second_xyz, first_xyz)
    deviations = np.linalg.norm(second_xyz @ rotation + translation - first_xyz, axis=1)
    residue_rmsd = {}
    for residue_number in ("21", "39", "88", "93", "129"):
        first_atoms = {atom_key(atom): atom["xyz"] for atom in parsed_first if atom["number"] == residue_number and atom["chain"] == "A" and atom["element"] != "H"}
        second_atoms = {atom_key(atom): atom["xyz"] for atom in parsed_second if atom["number"] == residue_number and atom["chain"] == "A" and atom["element"] != "H"}
        shared = sorted(first_atoms.keys() & second_atoms.keys())
        residue_rmsd[residue_number] = {"n_shared": len(shared), "RMSD_A": float(np.sqrt(np.mean([np.sum((second_atoms[key] @ rotation + translation - first_atoms[key]) ** 2) for key in shared])))}
    source_preservation = []
    for record, prepared in zip((first, second), (parsed_first, parsed_second)):
        source = Path(record["source_root"]).parents[1] / "structure_input/original.pdb"
        sources[str(source)] = digest(source)
        source_atoms = {atom_key(atom): atom for atom in atoms(source) if atom["number"] == "129" and atom["element"] != "H"}
        prepared_atoms = {atom_key(atom): atom for atom in prepared if atom["number"] == "129" and atom["element"] != "H"}
        shared = sorted(source_atoms.keys() & prepared_atoms.keys())
        source_preservation.append({"ligand": record["ligand_name"], "ASN129_source_to_prepared_RMSD_A":
            float(np.sqrt(np.mean([np.sum((source_atoms[key]["xyz"] - prepared_atoms[key]["xyz"]) ** 2) for key in shared])))})
    return {"n_CA": len(keys), "CA_RMSD_A": float(np.sqrt(np.mean(deviations ** 2))),
            "max_CA_deviation_A": float(deviations.max()), "local_residue_heavy_RMSD": residue_rmsd,
            "source_preservation": source_preservation}


def build(root, workspace):
    summary_path = root / "docs/SIZE_STRUCTURE_INTERACTIONS.json"
    summary = json.loads(summary_path.read_text())
    sources = {str(summary_path): digest(summary_path)}
    evidence = Path(summary["full_evidence"]["path"])
    if digest(evidence) != summary["full_evidence"]["sha256"]:
        raise ValueError("full evidence changed")
    with gzip.open(evidence, "rt") as handle:
        full = json.load(handle)
    original = json.loads((root / "docs/MOLECULAR_WEIGHT_DIAGNOSTICS.json").read_text())
    case_lookup = {(row["cohort"], row["id"]): row for row in original["cases"]}
    results = {"date": "2026-10-21", "new_QM": False, "production_changed": False,
               "local_geometry": {}, "pair_concentration": {}, "OpenBind_charge_audit": [], "sources_sha256": sources}
    for target, identifiers in (("P38", ("2ff", "2l")), ("BACE", ("CAT-4k", "CAT-4l"))):
        for identifier in identifiers:
            record = full["groups"][target]["records"][identifier]
            results["local_geometry"][target + "/" + identifier] = local_hydrogen_bonds(record, sources)
    results["OpenBind_308_318_receptor_alignment"] = align_openbind(full["groups"]["OpenBind"]["records"]["308"], full["groups"]["OpenBind"]["records"]["318"], sources)
    for target, group in full["groups"].items():
        ratios = []
        fragment_counts = {}
        for pair in group["eligible_pairs"]:
            if pair["correct"]:
                continue
            first, second = [group["records"][pair[key]] for key in ("first", "second")]
            if first["n_crystals"] != 1 or second["n_crystals"] != 1:
                continue
            left = {row["fragment"]: row["tpie"] for row in first["fragment_rows"]}
            right = {row["fragment"]: row["tpie"] for row in second["fragment_rows"]}
            changes = {key: right.get(key, 0) - left.get(key, 0) for key in left.keys() | right.keys()}
            dominant = max(changes, key=lambda key: abs(changes[key]))
            total = sum(abs(value) for value in changes.values())
            ratios.append(abs(changes[dominant]) / total if total else 0)
            fragment_counts[dominant] = fragment_counts.get(dominant, 0) + 1
        results["pair_concentration"][target] = {"n_wrong_single_crystal_pairs": len(ratios),
            "median_top_fragment_absolute_share": float(np.median(ratios)) if ratios else None,
            "dominant_fragment_counts": fragment_counts}
    group = full["groups"]["OpenBind"]
    changed_ids = {row["id"] for row in group["charge_representation_differences"]}
    paired_path = workspace / "runs/openbind/score_distribution_audit/paired_compounds.csv"
    paired = {row["compound_index"]: row for row in read_csv(paired_path)}
    for identifier in sorted(changed_ids):
        record = group["records"][identifier]
        source = Path(record["source_ligand_sdf"])
        molecule = next(molecule for molecule in Chem.SDMolSupplier(str(source), removeHs=False) if molecule.GetProp("_Name").strip() == record["ligand_name"])
        original_source = workspace / "runs/openbind/cases" / record["ligand_name"] / "openbind_ligands.sdf"
        original_molecule = next(iter(Chem.SDMolSupplier(str(original_source), removeHs=False)))
        sources[str(original_source)] = digest(original_source)
        case_path = original_source.parent / "case.json"
        raw_source = workspace / json.loads(case_path.read_text())["source_ligand"]
        sources[str(case_path)] = digest(case_path)
        sources[str(raw_source)] = digest(raw_source)
        raw_molecule = next(iter(Chem.SDMolSupplier(str(raw_source), removeHs=False)))
        charged_amide = Chem.MolFromSmarts("[CX3](=[OX1])[N+]")
        results["OpenBind_charge_audit"].append({"id": identifier, "ligand": record["ligand_name"],
            "charged_amide_N_matches": [list(match) for match in molecule.GetSubstructMatches(charged_amide)],
            "upstream_source_charge": Chem.GetFormalCharge(original_molecule),
            "raw_crystal_source_charge": Chem.GetFormalCharge(raw_molecule),
            "raw_crystal_source_smiles": Chem.MolToSmiles(Chem.RemoveHs(raw_molecule)),
            "prepared_charge": Chem.GetFormalCharge(molecule), "benchmark_charge": 0,
            "source_smiles": Chem.MolToSmiles(Chem.RemoveHs(molecule)),
            "original_source_smiles": Chem.MolToSmiles(Chem.RemoveHs(original_molecule))})
    metrics_by_subset = {}
    for label, identifiers in (("all490", set(group["records"])), ("unchanged484", set(group["records"]) - changed_ids)):
        ordered = sorted(identifiers)
        observed = [case_lookup["OpenBind", identifier]["label"] for identifier in ordered]
        metrics_by_subset[label] = {"n": len(ordered),
            "MW": metrics(observed, [case_lookup["OpenBind", identifier]["mw"] for identifier in ordered]),
            "FMO_frozen_softcap": metrics(observed, [case_lookup["OpenBind", identifier]["softcap_score"] for identifier in ordered])}
    results["OpenBind_charge_subset_sensitivity_not_corrected_performance"] = metrics_by_subset
    (root / "docs/SIZE_DIAGNOSIS_CONCLUSIONS.json").write_text(json.dumps(results, ensure_ascii=False, indent=2) + "\n")
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", required=True, type=Path)
    args = parser.parse_args()
    build(Path(__file__).resolve().parents[1], args.workspace.resolve())
