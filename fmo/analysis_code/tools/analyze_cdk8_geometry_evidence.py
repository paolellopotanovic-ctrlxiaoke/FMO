"""Geometry, pair attribution and FEP comparisons without score selection."""
from __future__ import annotations

import argparse
import csv
from itertools import combinations
import json
from pathlib import Path
import re
import sys

import numpy as np
from rdkit import Chem
from scipy.stats import spearmanr
from scipy.spatial import cKDTree

from run_cdk8_solution_cycle import array, digest, group, metrics, write_json


def xyz(deck):
    rows = group(deck, "FMOXYZ")[1].strip().splitlines()
    return [row.split()[0].upper() for row in rows], np.array([row.split()[2:] for row in rows], dtype=float)


def fragments(deck):
    content = group(deck, "FMO")[1]
    names = array(content, "FRGNAM", "INDAT")
    membership = list(map(int, array(content, "INDAT", r"(?!)")))
    result, current = [], []
    if membership[0] != 0:
        raise ValueError("unexpected fragment inventory")
    for index in membership[1:]:
        if index == 0:
            result.append(current)
            current = []
        elif index < 0:
            current.extend(range(current[-1] + 1, -index + 1))
        else:
            current.append(index)
    if current or len(names) != len(result):
        raise ValueError("incomplete fragment inventory")
    return names, result


def bond_audit(molecule, positions):
    records = []
    for bond in molecule.GetBonds():
        first, second = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        distance = float(np.linalg.norm(positions[first] - positions[second]))
        hydrogen = any(molecule.GetAtomWithIdx(index).GetAtomicNum() == 1 for index in (first, second))
        records.append({"atoms": [first, second], "distance_angstrom": distance,
            "flagged": distance < 0.55 or distance > (1.5 if hydrogen else 2.3)})
    return {"n_flagged": sum(record["flagged"] for record in records), "bonds": records,
            "warning": "Distance screen only, not a Hessian, protonation or bond-order validation"}


def map_receptor_atoms(path, positions, symbols):
    records = []
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            records.append({"atom_name": line[12:16].strip(), "residue": line[17:20].strip(),
                "chain": line[21], "residue_number": line[22:26].strip(),
                "element": line[76:78].strip().upper(),
                "xyz": [float(line[start:stop]) for start, stop in ((30,38),(38,46),(46,54))]})
    result = {}
    for symbol in set(symbols):
        candidates = [record for record in records if record["element"] == symbol]
        if not candidates:
            continue
        tree = cKDTree([record["xyz"] for record in candidates])
        for index, element in enumerate(symbols):
            if element != symbol:
                continue
            distance, nearest = tree.query(positions[index])
            if distance <= 0.002:
                result[index] = {**candidates[nearest], "mapping_distance": float(distance)}
    return result


def amide_torsions(molecule, positions):
    from rdkit.Chem import rdMolTransforms
    positions = np.asarray(positions, dtype=float)
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    for index, position in enumerate(positions):
        conformer.SetAtomPosition(index, position)
    result = []
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        oxygens = [neighbor for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 8
                   and molecule.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.BondType.DOUBLE]
        nitrogens = [neighbor for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 7]
        for oxygen in oxygens:
            for nitrogen in nitrogens:
                nitrogen_index = nitrogen.GetIdx()
                neighbor_indices = [neighbor.GetIdx() for neighbor in nitrogen.GetNeighbors()]
                angle_sum = None
                if len(neighbor_indices) == 3:
                    vectors = positions[neighbor_indices] - positions[nitrogen_index]
                    angles = [float(np.degrees(np.arccos(np.clip(np.dot(vectors[first], vectors[second]) /
                        (np.linalg.norm(vectors[first])*np.linalg.norm(vectors[second])), -1, 1))))
                        for first,second in combinations(range(3),2)]
                    angle_sum = sum(angles)
                for substituent in nitrogen.GetNeighbors():
                    if substituent.GetIdx() == atom.GetIdx():
                        continue
                    indices = [oxygen.GetIdx(), atom.GetIdx(), nitrogen.GetIdx(), substituent.GetIdx()]
                    angle = rdMolTransforms.GetDihedralDeg(conformer, *indices)
                    deviation = min(abs(angle), abs(180 - abs(angle)))
                    result.append({"atoms_zero_based": indices, "angle_deg": angle,
                        "departure_from_planarity_deg": deviation,
                        "nitrogen_neighbor_angle_sum_deg": angle_sum,
                        "warning": "Dihedral departure can reflect both amide twist and nitrogen pyramidalization"})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    args = parser.parse_args()
    workspace, reference, run = args.workspace.resolve(), args.reference.resolve(), args.run.resolve()
    sys.path.insert(0, str(workspace / "src"))
    from sophosqm_baseline.parser import parse_pie
    evaluation = json.loads((run / "evaluation.json").read_text())
    if evaluation.get("status") != "COMPLETE":
        raise ValueError("full cohort must be complete before comparative report")
    dataset = workspace / "vendor/public_binding_free_energy_benchmark"
    sdf = dataset / "fep_benchmark_inputs/structure_inputs/misc/cdk8_koehler_ligands.sdf"
    molecules = {mol.GetProp("_Name").strip(): mol for mol in Chem.SDMolSupplier(str(sdf), removeHs=False)}
    fep_path = dataset / "21_4_results/ligand_predictions/misc/cdk8_koehler_213_out.csv"
    with fep_path.open() as handle:
        fep = {row["Ligand name"]: row for row in csv.DictReader(handle)}
    previous = {row["ligand"]: row for row in json.loads((reference / "final_evaluation.json").read_text())["rows"]}
    sampled = {row["ligand"]: row for row in json.loads((reference / "solution_sampling4/final_evaluation.json").read_text())["rows"]}
    rows = evaluation["rows"]
    if set(fep) != {row["ligand"] for row in rows}:
        raise ValueError("FEP identity mismatch")
    geometry, attribution = {}, {}
    baseline_root = workspace / "FMOtpiefinal_HIV1lab/runs/all_accelerated_robust_20261020/CDK8/run/qm/full_target_fmo"
    for row in rows:
        name = row["ligand"]
        deck = (reference / f"pcm/{name}/complex/input.inp").read_text()
        symbols, initial = xyz(deck)
        fragment_names, membership = fragments(deck)
        ligand_size = len(membership[-1])
        summary = json.loads((run / name / "summary.json").read_text())
        final = np.array([line.split()[2:] for line in summary["optimized_coordinates"]], dtype=float)
        mol = molecules[name]
        if [atom.GetSymbol().upper() for atom in mol.GetAtoms()] != symbols[-ligand_size:]:
            raise ValueError("SDF/deck ligand order mismatch")
        heavy = np.array([atom.GetAtomicNum() > 1 for atom in mol.GetAtoms()])
        delta = final[-ligand_size:] - initial[-ligand_size:]
        original = np.array(mol.GetConformer().GetPositions())
        initial_ligand, final_ligand = initial[-ligand_size:], final[-ligand_size:]
        semantic = map_receptor_atoms(baseline_root / "cdk8_prepared.pdb", initial[:-ligand_size], symbols[:-ligand_size])
        polar_contacts = []
        nearest_heavy_contacts = []
        for receptor_index, atom in semantic.items():
            if atom["element"] != "H":
                for ligand_atom in mol.GetAtoms():
                    if ligand_atom.GetAtomicNum() == 1:
                        continue
                    ligand_index = ligand_atom.GetIdx()
                    nearest_heavy_contacts.append({"receptor": atom,
                        "ligand_atom_zero_based": ligand_index,
                        "k10_distance": float(np.linalg.norm(initial[receptor_index]-initial_ligand[ligand_index])),
                        "qm_distance": float(np.linalg.norm(final[receptor_index]-final_ligand[ligand_index]))})
            if atom["element"] not in ("N", "O", "S"):
                continue
            for ligand_atom in mol.GetAtoms():
                if ligand_atom.GetAtomicNum() not in (7, 8, 16):
                    continue
                ligand_index = ligand_atom.GetIdx()
                before = float(np.linalg.norm(initial[receptor_index] - initial_ligand[ligand_index]))
                after = float(np.linalg.norm(final[receptor_index] - final_ligand[ligand_index]))
                if min(before, after) < 4:
                    polar_contacts.append({"receptor": atom, "ligand_atom_zero_based": ligand_index,
                        "ligand_element": ligand_atom.GetSymbol(), "k10_distance": before, "qm_distance": after,
                        "warning": "Polar distance only; donor/acceptor roles and angles not implied"})
        closest = {}
        for label in ("ALYS52", "AALA100", "AARG356", "AASP173"):
            if label not in fragment_names:
                continue
            indices = membership[fragment_names.index(label)]
            receptor_indices = [index - 1 for index in indices if symbols[index - 1] != "H"]
            ligand_indices = np.where(heavy)[0]
            distances = {}
            for state, positions in (("k10", initial_ligand), ("qm", final_ligand)):
                contacts = [(float(np.linalg.norm(final[index] - positions[ligand_index])), index + 1,
                    int(ligand_index), symbols[index], mol.GetAtomWithIdx(int(ligand_index)).GetSymbol())
                    for index in receptor_indices for ligand_index in ligand_indices]
                distances[state] = sorted(contacts)[:5]
            closest[label] = distances
        geometry[name] = {"heavy_rmsd_qm_vs_k10_no_alignment": float(np.sqrt(np.mean(np.sum(delta[heavy]**2, axis=1)))),
            "heavy_rmsd_k10_vs_input_no_alignment": float(np.sqrt(np.mean(np.sum((initial_ligand[heavy]-original[heavy])**2, axis=1)))),
            "heavy_rmsd_qm_vs_input_no_alignment": float(np.sqrt(np.mean(np.sum((final_ligand[heavy]-original[heavy])**2, axis=1)))),
            "max_frozen_receptor_displacement": float(np.max(np.linalg.norm(final[:-ligand_size]-initial[:-ligand_size], axis=1))),
            "bond_audit": bond_audit(mol, final_ligand), "nearest_fragment_contacts": closest,
            "mapped_receptor_atoms": len(semantic), "polar_contacts": polar_contacts,
            "nearest_heavy_contacts_k10": sorted(nearest_heavy_contacts, key=lambda record: record["k10_distance"])[:10],
            "nearest_heavy_contacts_qm": sorted(nearest_heavy_contacts, key=lambda record: record["qm_distance"])[:10],
            "amide_geometry": {state: amide_torsions(mol, coordinates) for state, coordinates in (
                ("input", original), ("k10", initial_ligand), ("qm_bound", final_ligand))},
            "fragment_boundary_warning": "FACIO fragment membership is not identical to whole-residue atom membership"}
        components = {}
        for state, path in (("k10", baseline_root / f"qm/{name}/run/gamess.log"),
                            ("qm", run / name / "pieda/gamess.log")):
            pairs = [pair for pair in parse_pie(path) if len(fragment_names) in (pair.fragment_i, pair.fragment_j)]
            components[state] = {fragment_names[(pair.fragment_j if pair.fragment_i == len(fragment_names) else pair.fragment_i)-1]: pair.value_kcal_mol for pair in pairs}
        attribution[name] = components
        row["fep_archived_dg"] = float(fep[name]["Pred. dG (kcal/mol)"])
        if abs(row["experimental_dg"] - float(fep[name]["Exp. dG (kcal/mol)"])) > 1e-4:
            raise ValueError("FEP experimental labels mismatch")
        row["original_multistart_strain"] = sampled[name]["solution_min_strain"]
        row["original_three_state"] = previous[name]["three_state_pcm_bound"]
    arms = list(evaluation["metrics"]) + ["fep_archived_dg"]
    truth = [row["experimental_dg"] for row in rows]
    pair_rows = []
    for first, second in combinations(rows, 2):
        difference = first["experimental_dg"] - second["experimental_dg"]
        pair = {"first": first["ligand"], "second": second["ligand"], "experimental_delta": difference,
                "experimental_gap_at_least_0p3": abs(difference) >= 0.3}
        for arm in arms:
            predicted = first[arm] - second[arm]
            pair[arm] = {"delta": predicted, "direction_correct": bool(predicted * difference > 0)}
        pair_rows.append(pair)
    summary_metrics = {arm: metrics(truth, [row[arm] for row in rows]) for arm in arms}
    for arm in arms:
        estimates = [float(spearmanr([row["experimental_dg"] for offset, row in enumerate(rows) if offset != index],
            [row[arm] for offset, row in enumerate(rows) if offset != index]).statistic) for index in range(len(rows))]
        summary_metrics[arm]["leave_one_out_spearman_range_diagnostic_only"] = [min(estimates), max(estimates)]
        summary_metrics[arm]["correct_pairs_of_45"] = sum(pair[arm]["direction_correct"] for pair in pair_rows)
    evidence = {"n": len(rows), "metrics": summary_metrics, "rows": rows, "geometry": geometry,
        "fragment_contributions": attribution, "pairs": pair_rows,
        "experimental_range_kcal_mol": max(truth)-min(truth), "production_changed": False,
        "sources_sha256": {str(path): digest(path) for path in (sdf, fep_path, run / "evaluation.json")},
        "limitations": ["Development set, not external validation", "No measured experimental uncertainty supplied",
            "0.3 kcal gap is a diagnostic threshold, not an assumed experimental uncertainty",
            "Converged gradient is not Hessian-confirmed minimum", "PIEDA fragment attribution is not causal proof",
            "Input SDF poses are not assumed independent crystal structures for all ten ligands"]}
    write_json(run / "geometry_evidence.json", evidence)
    lines = ["# CDK8：同模型几何对照与FEP比较", "", "完整10例；不删例、不拟合；生产主流程未替换。", "",
        "|方案|Pearson|Spearman|Kendall|正确分子对/45|", "|---|---:|---:|---:|---:|"]
    lines.extend(f"|{arm}|{value['pearson']:.3f}|{value['spearman']:.3f}|{value['kendall']:.3f}|{value['correct_pairs_of_45']}|" for arm, value in summary_metrics.items())
    lines += ["", "## 逐例变化", "", "|配体|实验DG|原tPIE|QM tPIE|原应变|QM应变|重原子RMSD Å|", "|---|---:|---:|---:|---:|---:|---:|"]
    lines.extend(f"|{row['ligand']}|{row['experimental_dg']:.3f}|{row['original_tpie']:.3f}|{row['qm_bound_tpie']:.3f}|{row['original_multistart_strain']:.3f}|{row['qm_strain_kcal_mol']:.3f}|{geometry[row['ligand']]['heavy_rmsd_qm_vs_k10_no_alignment']:.3f}|" for row in rows)
    lines += ["", "逐分子对、残基片段贡献、接触距离及限制见geometry_evidence.json。",
        "较低能量或较小应变不自动等于更好的亲和力排序；FEP为既有归档结果，非本轮新跑。"]
    (run / "GEOMETRY_REPORT.md").write_text("\n".join(lines) + "\n")
    print(json.dumps(summary_metrics, indent=2))


if __name__ == "__main__":
    main()
