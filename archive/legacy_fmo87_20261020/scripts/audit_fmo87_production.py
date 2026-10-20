#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from sophosqm_baseline.parser import aggregate_tpie, find_selected_tie, parse_pie


EXPECTED_COUNTS = {"CDK2": 16, "JNK1": 21, "P38": 34, "TYK2": 16}
EXPECTED_COHORT_SHA256 = "19edfb8e2ccf3ab1029d848e4062d70cb35532d2b4fee63674f160b2a292997b"
READYINESS = "GAMESS_DECK_GENERATED_FULL_ACE_NME_SHIFTED_C_ALPHA_C"


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _prepared_water_residues(path: Path) -> int:
    residues: set[tuple[str, str, str]] = set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line[:6].strip() != "HETATM" or line[17:20].strip().upper() not in {"HOH", "WAT", "TIP"}:
                continue
            residues.add((line[21], line[22:27], line[17:20].strip().upper()))
    return len(residues)


def audit_preparation(
    root: Path,
    *,
    expected_policy: str,
    prepared_root: Path | None = None,
) -> dict:
    records = {}
    prepared_root = prepared_root or root.parent.parent / "prepared_proteins"
    for target in EXPECTED_COUNTS:
        log_path = next((
            candidate for candidate in (
                prepared_root / f"{target}_pdb2pqr.log",
                prepared_root / f"{target}_prepared.pdb2pqr.log",
            ) if candidate.exists()
        ), None)
        if log_path is None:
            raise FileNotFoundError(f"missing PDB2PQR log for {target}")
        log = log_path.read_text(encoding="utf-8", errors="replace")
        prepared_path = prepared_root / f"{target}_prepared.pdb"
        preparation_path = prepared_root / f"{target}_preparation.json"
        preparation = json.loads(preparation_path.read_text(encoding="utf-8")) if preparation_path.exists() else {}
        explicit_water_residues = _prepared_water_residues(prepared_path)
        records[target] = {
            "pdb2pqr_3_6_1": "PDB2PQR v3.6.1" in log,
            "propka": "Assigning titration states with PROPKA" in log,
            "ph_7_4": "Applying pKa values at a pH of 7.40" in log,
            "hydrogen_bond_optimization": "Optimizing hydrogen bonds" in log,
            "prepared_pdb_sha256": file_sha256(prepared_path),
            "explicit_water_residues": explicit_water_residues,
            "explicit_water_policy": preparation.get("explicit_water_policy"),
            "explicit_water_policy_matches": (
                expected_policy == "removed_before_preparation" and explicit_water_residues == 0
            ),
        }
    return records


def audit_compound(row: dict, root: Path) -> dict:
    key = f"{row['target']}_{row['ligand_name']}"
    geometry_dir = root / "minimized" / key
    job_dir = root / "fmo" / key
    record = {"key": key}
    try:
        geometry = json.loads((geometry_dir / "minimization.json").read_text(encoding="utf-8"))
        structure_manifest = json.loads((job_dir / "structure_manifest.json").read_text(encoding="utf-8"))
        record["geometry"] = {
            "pocket_cutoff_5_angstrom": structure_manifest["cutoff_angstrom"] == 5.0,
            "protein_coordinate_policy_free": geometry["protein_coordinate_policy"] == "free",
            "ligand_restraint_100": geometry["restraint_k_kj_mol_angstrom2"] == 100.0,
            "amber14": geometry["protein_force_field"] == "AMBER14 protein.ff14SB",
            "openff_2_2_0": geometry["ligand_force_field"] == "OpenFF 2.2.0",
        }
        record["geometry_metrics"] = {
            "explicit_water_residues": geometry.get("explicit_water_residues"),
            "ligand_rmsd_angstrom": geometry["ligand_rmsd_angstrom"],
            "ligand_max_displacement_angstrom": geometry["ligand_max_displacement_angstrom"],
        }
        fragment_manifest = json.loads((job_dir / "fragment_manifest.json").read_text(encoding="utf-8"))
        fragments = fragment_manifest["fragments"]
        preparation = fragment_manifest["metadata"]["preparation"]
        result = json.loads((job_dir / "run" / "result.json").read_text(encoding="utf-8"))
        deck = (job_dir / "job.inp").read_text(encoding="utf-8")
        deck_lines = deck.splitlines()
        log = (job_dir / "run" / "gamess.log").read_text(encoding="utf-8", errors="replace")
        pairs = parse_pie(job_dir / "run" / "gamess.log")
        ligand_fragment = fragment_manifest["metadata"]["ligand_fragment_id"]
        tpie = aggregate_tpie(pairs, ligand_fragment)
        official_tie = find_selected_tie(job_dir / "run" / "gamess.log", ligand_fragment)
        record["fmo"] = {
            "status_succeeded": result["status"] == "SUCCEEDED",
            "readiness_exact": fragment_manifest["metadata"]["readiness"] == READYINESS,
            "nme_replaces_oxt": "nme_replaced_oxt_atom_indices" in preparation,
            "full_terminal_caps": fragment_manifest["metadata"]["preparation"]["terminal_caps"] ==
                                   "full neutral ACE/NME residues at selected segment ends",
            "hop_count": len(fragment_manifest["metadata"]["preparation"]["hop_pairs_gamess_indices"]),
            "deck_lines_within_80": max(map(len, deck_lines)) <= 80,
            "dftb3_3ob": "$DFTB SCC=.TRUE. DFTB3=.TRUE. DAMPXH=.TRUE. DAMPEX=4.00" in deck and "$DFTBSK" in deck,
            "pieda_modorb": "IPIEDA=1" in deck and "MODORB=3" in deck,
            "all_dimers": "MODMOL=1" in deck and "NACUT=0" in deck,
            "water_pcm": "$PCM SOLVNT=WATER" in deck and "IFMO=-1" in deck,
            "hop_bonds": "$FMOBND" in deck and " HOP_C" in deck,
            "all_fragments_closed_shell": all(fragment["multiplicity"] == 1 for fragment in fragments),
            "normal_termination": "EXECUTION OF GAMESS TERMINATED NORMALLY" in log,
            "no_scf_failure": "SCF IS UNCONVERGED" not in log and "DIMER CALCULATION DID NOT CONVERGE" not in log,
            "no_geometry_failure": "Distance errors" not in log and "impossible!" not in log,
            "pair_count": len(pairs),
            "ligand_pair_count": sum(
                (pair.fragment_i == ligand_fragment) != (pair.fragment_j == ligand_fragment)
                for pair in pairs
            ),
            "fragment_count": len(fragments),
            "water_fragment_count": sum(fragment["kind"] == "water" for fragment in fragments),
            "selected_pie_pair_count_expected": len(pairs) == len(fragments) - 1,
            "ligand_pair_count_expected": sum(
                (pair.fragment_i == ligand_fragment) != (pair.fragment_j == ligand_fragment)
                for pair in pairs
            ) == len(fragments) - 1,
            "tpie_pair_sum": tpie,
            "official_tie": official_tie,
            "tie_reconciles": official_tie is not None and math.isclose(tpie, official_tie, abs_tol=0.005),
            "result_tpie": result.get("tpie_kcal_mol"),
            "input_sha_matches": result.get("input_sha256") == file_sha256(job_dir / "job.inp"),
        }
        record["passed"] = (
            all(record["geometry"].values()) and all(
                value for key, value in record["fmo"].items()
                if isinstance(value, bool)
            ) and record["fmo"]["tie_reconciles"] and record["fmo"]["input_sha_matches"]
        )
    except Exception as error:
        record["error"] = f"{type(error).__name__}: {error}"
        record["passed"] = False
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", default="data/fmo87_cohort.csv")
    parser.add_argument("--root", default="runs/fmo87/fmo_rebuild/targets")
    parser.add_argument("--output-dir", default="runs/fmo87/final")
    parser.add_argument("--require-complete", action="store_true")
    parser.add_argument(
        "--explicit-water-policy",
        choices=["retained", "removed_from_fmo_pocket", "removed_before_preparation"],
        default="retained",
    )
    parser.add_argument(
        "--prepared-proteins",
        type=Path,
        help="explicit prepared-protein root; defaults to two directories above the audit root",
    )
    args = parser.parse_args()
    cohort_path = Path(args.cohort)
    root = Path(args.root)
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    with cohort_path.open(encoding="utf-8", newline="") as handle:
        cohort = list(csv.DictReader(handle))
    counts = {target: sum(row["target"] == target for row in cohort) for target in EXPECTED_COUNTS}
    compound_records = []
    for row in cohort:
        compound_records.append(audit_compound(row, root / row["target"]))
    complete = len(compound_records) == 87 and counts == EXPECTED_COUNTS
    passed = sum(item["passed"] for item in compound_records)
    protein_preparation = audit_preparation(
        root,
        expected_policy=args.explicit_water_policy,
        prepared_root=args.prepared_proteins,
    )
    explicit_water_ablation = {
        "policy": args.explicit_water_policy,
        "all_prepared_proteins_zero_water": (
            all(item["explicit_water_residues"] == 0 for item in protein_preparation.values())
            if args.explicit_water_policy == "removed_before_preparation" else None
        ),
        "all_fmo_pockets_zero_water": (
            all(item.get("fmo", {}).get("water_fragment_count", 0) == 0 for item in compound_records)
            if args.explicit_water_policy != "retained" else None
        ),
    }
    explicit_water_complete = True
    if explicit_water_ablation["all_prepared_proteins_zero_water"] is not None:
        explicit_water_complete = (
            explicit_water_complete and explicit_water_ablation["all_prepared_proteins_zero_water"]
        )
    if explicit_water_ablation["all_fmo_pockets_zero_water"] is not None:
        explicit_water_complete = (
            explicit_water_complete and explicit_water_ablation["all_fmo_pockets_zero_water"]
        )
    audit = {
        "cohort": {
            "n": len(cohort), "counts": counts,
            "all_neutral": all(int(row["formal_charge"]) == 0 for row in cohort),
            "sha256": file_sha256(cohort_path), "sha256_matches_frozen": file_sha256(cohort_path) == EXPECTED_COHORT_SHA256,
        },
        "protein_preparation": protein_preparation,
        "n_compounds": len(compound_records), "n_passed": passed,
        "strict_complete": complete and passed == 87 and explicit_water_complete,
        "compounds": compound_records,
        "protocol_deviations": [
            "Ross/Chen reference poses replace unavailable Glide SP pose ensembles.",
            "OpenMM AMBER14/OpenFF 2.2.0 replaces unavailable MacroModel/OPLS3e.",
            "Water PCM is required by the self-developed GAMESS DFTB3 protocol for stable polar dimers.",
            "RDKit Crippen clogP is used only in secondary SophosQM-style fits, never raw Table 12 comparison.",
        ],
        "explicit_water_ablation": explicit_water_ablation,
    }
    (output / "production_audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    failures = [item for item in compound_records if not item["passed"]]
    report = [
        "# FMO87 production audit", "",
        f"- Cohort SHA-256: `{audit['cohort']['sha256']}`",
        f"- Compounds passed: {passed}/{len(compound_records)}",
        f"- Explicit-water policy: `{args.explicit_water_policy}`",
        f"- Strict complete: `{audit['strict_complete']}`", "",
    ]
    if failures:
        report.extend(["## Failures", ""])
        report.extend(f"- `{item['key']}`: {item.get('error', 'boolean gate failed')}" for item in failures)
    (output / "PRODUCTION_AUDIT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(json.dumps({"n_compounds": len(compound_records), "n_passed": passed,
                      "strict_complete": audit["strict_complete"],
                      "audit": str(output / "production_audit.json")}, indent=2))
    if args.require_complete and not audit["strict_complete"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
