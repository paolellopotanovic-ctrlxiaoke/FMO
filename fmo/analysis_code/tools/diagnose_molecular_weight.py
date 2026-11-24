"""Describe size confounding without changing or fitting production scores."""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.stats import pearsonr

from build_model_comparison import COHORTS, check_metrics, digest, metrics, read_csv


def pearson(first, second):
    if len(first) < 3 or np.std(first) < 1e-10 or np.std(second) < 1e-10:
        return None
    return float(pearsonr(first, second).statistic)


def center(values, groups):
    result = np.asarray(values, dtype=float).copy()
    for group in set(groups):
        selected = np.asarray([value == group for value in groups])
        result[selected] -= result[selected].mean()
    return result


def partial(first, second, mass):
    design = np.column_stack([np.ones(len(mass)), mass])
    residuals = [values - design @ np.linalg.lstsq(design, values, rcond=None)[0]
                 for values in (first, second)]
    return pearson(*residuals)


def pair_order(labels, scores, mass, tolerance):
    first, second = np.triu_indices(len(labels), k=1)
    selected = np.abs(mass[first] - mass[second]) <= tolerance
    products = (labels[first] - labels[second]) * (scores[first] - scores[second])
    correct = int(np.sum(selected & (products > 0)))
    incorrect = int(np.sum(selected & (products < 0)))
    tied = int(np.sum(selected & (products == 0)))
    return {"pairs": int(selected.sum()), "concordant": correct, "discordant": incorrect,
            "ties": tied, "accuracy_excluding_ties": correct / (correct + incorrect) if correct + incorrect else None}


def bootstrap(first, second, repeats=3000):
    random = np.random.default_rng(20261021)
    draws = random.integers(0, len(first), size=(repeats, len(first)))
    left, right = first[draws], second[draws]
    left, right = left - left.mean(axis=1, keepdims=True), right - right.mean(axis=1, keepdims=True)
    denominator = np.sqrt(np.sum(left ** 2, axis=1) * np.sum(right ** 2, axis=1))
    valid = denominator > 1e-10
    correlations = np.sum(left[valid] * right[valid], axis=1) / denominator[valid]
    return [float(value) for value in np.quantile(correlations, [0.025, 0.975])]


def describe(rows):
    labels = np.asarray([row["label"] for row in rows])
    mass = np.asarray([row["mw"] for row in rows])
    scores = np.asarray([row["softcap_score"] for row in rows])
    raw = np.asarray([row["raw_score"] for row in rows])
    scaffolds = [row["scaffold"] for row in rows]
    counts = Counter(scaffolds)
    bins = np.searchsorted(np.unique(np.quantile(mass, [1 / 3, 2 / 3])), mass, side="right").tolist()
    centered_label = center(labels, bins)
    repeated = np.asarray([counts[scaffold] >= 3 and bool(scaffold) for scaffold in scaffolds])
    selected_scaffolds = [scaffold for scaffold, selected in zip(scaffolds, repeated) if selected]
    same_scaffold = None
    if repeated.sum() >= 3:
        same_scaffold = {"n": int(repeated.sum()), "groups": len(set(selected_scaffolds)),
                         "MW_centered_pearson": pearson(center(labels[repeated], selected_scaffolds), center(mass[repeated], selected_scaffolds)),
                         "FMO_centered_pearson": pearson(center(labels[repeated], selected_scaffolds), center(scores[repeated], selected_scaffolds))}
    scaffold_summaries = []
    for scaffold, count in counts.most_common():
        if count < 5 or not scaffold:
            continue
        selected = np.asarray([value == scaffold for value in scaffolds])
        scaffold_summaries.append({"scaffold": scaffold, "n": count,
            "MW_label_pearson": pearson(labels[selected], mass[selected]),
            "FMO_label_pearson": pearson(labels[selected], scores[selected])})
    leave_one = [pearson(np.delete(labels, index), np.delete(mass, index)) for index in range(len(rows))]
    leave_one = [value for value in leave_one if value is not None]
    descriptors = {}
    for name in ("heavy_atoms", "logp", "rotatable_bonds", "aromatic_rings", "formal_charge", "halogens"):
        values = np.asarray([row[name] for row in rows])
        descriptors[name] = {"label_pearson": pearson(labels, values), "MW_pearson": pearson(mass, values)}
    return {"n": len(rows), "MW_range": [float(mass.min()), float(mass.max())],
        "MW_IQR": float(np.subtract(*np.quantile(mass, [.75, .25]))),
        "MW_label": metrics(labels.tolist(), mass.tolist()),
        "MW_pearson_bootstrap_95": bootstrap(labels, mass),
        "MW_pearson_leave_one_out_range": [min(leave_one), max(leave_one)],
        "FMO_softcap_label": metrics(labels.tolist(), scores.tolist()),
        "FMO_raw_label": metrics(labels.tolist(), raw.tolist()),
        "MW_FMO_softcap_pearson": pearson(mass, scores),
        "FMO_softcap_partial_pearson_linear_MW": partial(labels, scores, mass),
        "FMO_raw_partial_pearson_linear_MW": partial(labels, raw, mass),
        "tertiles": {"counts": dict(Counter(bins)),
            "MW_centered_pearson": pearson(centered_label, center(mass, bins)),
            "FMO_centered_pearson": pearson(centered_label, center(scores, bins))},
        "within_10Da": {"MW": pair_order(labels, mass, mass, 10), "FMO": pair_order(labels, scores, mass, 10)},
        "within_0_05Da": {"MW": pair_order(labels, mass, mass, .05), "FMO": pair_order(labels, scores, mass, .05)},
        "distinct_scaffolds": len(counts), "largest_scaffold_n": max(counts.values()),
        "same_scaffold_centered": same_scaffold, "scaffold_groups_n_ge_5": scaffold_summaries,
        "descriptors": descriptors}


def build(workspace, root):
    baseline_path = root / "docs/MOLECULAR_WEIGHT_COMPARISON.json"
    baseline = json.loads(baseline_path.read_text())
    sources = {str(baseline_path): digest(baseline_path)}
    original = read_csv(workspace / "data/fmo87_cohort.csv")
    paired = read_csv(workspace / "runs/openbind/score_distribution_audit/paired_compounds.csv")
    molecules = {}
    for target, (_, stem, _) in COHORTS.items():
        if target in ("CDK2", "JNK1", "P38", "TYK2"):
            molecules[target] = {row["ligand_name"]: Chem.MolFromSmiles(row["canonical_isomeric_smiles"])
                                 for row in original if row["target"] == target}
        else:
            source = workspace / "vendor/public_binding_free_energy_benchmark/fep_benchmark_inputs/structure_inputs" / f"{stem}_ligands.sdf"
            sources[str(source)] = digest(source)
            molecules[target] = {molecule.GetProp("_Name").strip(): molecule
                                 for molecule in Chem.SDMolSupplier(str(source), removeHs=False)}
    molecules["OpenBind"] = {row["compound_index"]: Chem.MolFromSmiles(row["compound_smiles"]) for row in paired}
    run = workspace / "FMOtpiefinal_HIV1lab/runs/all_accelerated_robust_20261020"
    source = run / "evaluation.json"
    sources[str(source)] = digest(source)
    evaluation = json.loads(source.read_text())
    scores = {}
    for target in COHORTS:
        source = run / target / "independent_verification.json"
        sources[str(source)] = digest(source)
        scores[target] = {row["ligand"]: row for row in json.loads(source.read_text())["records"]}
    scores["OpenBind"] = {str(row["compound_index"]): row for row in evaluation["groups"]["OpenBind"]["compound_predictions"]}
    all_rows = []
    groups = {}
    for target in baseline["groups"]:
        rows = []
        for case in baseline["cases"]:
            if case["cohort"] != target:
                continue
            molecule = Chem.RemoveHs(molecules[target][case["id"]])
            if abs(Descriptors.MolWt(molecule) - case["mw"]) > .05:
                raise ValueError(f"MW identity mismatch: {target}/{case['id']}")
            prediction = scores[target][case["id"]]
            rows.append({**case, "raw_score": prediction["raw_score"], "softcap_score": prediction["softcap_score"],
                "heavy_atoms": molecule.GetNumHeavyAtoms(), "logp": Crippen.MolLogP(molecule),
                "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(molecule),
                "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(molecule),
                "formal_charge": Chem.GetFormalCharge(molecule),
                "halogens": sum(atom.GetAtomicNum() in (9, 17, 35, 53) for atom in molecule.GetAtoms()),
                "scaffold": MurckoScaffold.MurckoScaffoldSmiles(mol=molecule)})
        if len(rows) != baseline["groups"][target]["n"]:
            raise ValueError(f"coverage mismatch: {target}")
        groups[target] = describe(rows)
        check_metrics(groups[target]["MW_label"], baseline["groups"][target]["MW_only"])
        check_metrics(groups[target]["FMO_raw_label"], baseline["groups"][target]["FMO_raw"])
        all_rows.extend(rows)
    result = {"date": "2026-10-21", "production_modified": False,
        "warning": "Residualization uses labels only for descriptive diagnostics, not an unfitted deployable predictor. Pairs are dependent; bootstrap assumes independent compounds, not scaffold-family independence. Tertiles retain within-bin MW variation. No causal attribution or multiplicity-adjusted significance claim.",
        "groups": groups, "cases": all_rows, "sources_sha256": sources}
    (root / "docs/MOLECULAR_WEIGHT_DIAGNOSTICS.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", required=True, type=Path)
    arguments = parser.parse_args()
    build(arguments.workspace.resolve(), Path(__file__).resolve().parents[1])
