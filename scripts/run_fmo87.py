#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys
import time
import shutil

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rdkit import Chem
from rdkit.Chem import Crippen

from sophosqm_baseline.jobs import prepare_job
from sophosqm_baseline.cohort import resolve_structure_path
from sophosqm_baseline.minimization import minimize_complex
from sophosqm_baseline.runner import run_gamess


def _existing_success(
    result_path: Path,
    job_dir: Path,
    minimized_dir: Path,
    explicit_water_policy: str,
) -> dict | None:
    if not result_path.exists():
        return None
    result = json.loads(result_path.read_text(encoding="utf-8"))
    if result.get("status") != "SUCCEEDED":
        return None
    structure_manifest = json.loads(
        (job_dir / "structure_manifest.json").read_text(encoding="utf-8")
    )
    expected_pocket_water = explicit_water_policy == "retained"
    if structure_manifest.get("include_explicit_waters") is not expected_pocket_water:
        return None
    if explicit_water_policy == "removed_before_preparation":
        geometry = json.loads(
            (minimized_dir / "minimization.json").read_text(encoding="utf-8")
        )
        if geometry.get("explicit_water_residues") not in {None, 0}:
            return None
    return result


def _validate_preparation_policy(prepared_proteins: Path, targets: list[str], policy: str) -> None:
    expected_keep_water = policy != "removed_before_preparation"
    for target in sorted(set(targets)):
        metadata_path = prepared_proteins / f"{target}_preparation.json"
        if not metadata_path.exists():
            raise FileNotFoundError(f"missing preparation provenance: {metadata_path}")
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if metadata.get("keep_water") is not expected_keep_water:
            raise ValueError(
                f"{target} preparation keep_water={metadata.get('keep_water')!r} "
                f"does not match explicit-water policy {policy!r}"
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", default="data/fmo87_cohort.csv")
    parser.add_argument(
        "--structure-root",
        default="data/fmo87_structures",
        help="directory containing the frozen cohort ligand SDF and protein PDB files",
    )
    parser.add_argument("--output", default="runs/fmo87/fmo_dftb3")
    parser.add_argument("--prepared-proteins", default="runs/fmo87/prepared_proteins")
    parser.add_argument("--gamess-dir", required=True)
    parser.add_argument("--library-dir", required=True)
    parser.add_argument("--dftb-parameter-dir", required=True)
    parser.add_argument("--timeout", type=float, default=1800)
    parser.add_argument("--stage", choices=["minimize", "fmo", "all"], default="all")
    parser.add_argument("--only-target")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--shard-index", type=int)
    parser.add_argument("--shard-count", type=int)
    parser.add_argument(
        "--explicit-water-policy",
        choices=["retained", "removed_from_fmo_pocket", "removed_before_preparation"],
        default="retained",
    )
    parser.add_argument("--drop-waters", action="store_true",
                        help="exclude explicit waters from the FMO pocket while retaining the minimized complex geometry")
    args = parser.parse_args()

    root = Path(args.output)
    params = Path(args.dftb_parameter_dir).resolve()
    with open(args.cohort, encoding="utf-8", newline="") as handle:
        records = list(csv.DictReader(handle))
    structure_root = Path(args.structure_root)
    for record in records:
        record["ligand_sdf"] = str(resolve_structure_path(
            record["ligand_sdf"], structure_root
        ))
        record["protein_pdb"] = str(resolve_structure_path(
            record["protein_pdb"], structure_root
        ))
    if args.only_target:
        records = [record for record in records if record["target"] == args.only_target]
    if args.limit:
        records = records[:args.limit]
    if (args.shard_index is None) != (args.shard_count is None):
        parser.error("--shard-index and --shard-count must be used together")
    if args.shard_count is not None:
        if args.shard_index < 0 or args.shard_count <= 0 or args.shard_index >= args.shard_count:
            parser.error("invalid shard index/count")
        records = [record for index, record in enumerate(records) if index % args.shard_count == args.shard_index]
    if args.drop_waters and args.explicit_water_policy == "retained":
        args.explicit_water_policy = "removed_from_fmo_pocket"
    if args.drop_waters and args.explicit_water_policy == "removed_before_preparation":
        parser.error("--drop-waters removes waters only from the FMO pocket, not before PDB2PQR/PROPKA")
    _validate_preparation_policy(
        Path(args.prepared_proteins), [record["target"] for record in records],
        args.explicit_water_policy,
    )
    rows = []

    for index, record in enumerate(records, start=1):
        key = f"{record['target']}_{record['ligand_name']}"
        minimized_dir = root / "minimized" / key
        complex_pdb = minimized_dir / "minimized_complex.pdb"
        job_dir = root / "fmo" / key
        row = dict(record)
        row.update({"job": key, "status": "PENDING"})
        try:
            if args.stage in {"minimize", "all"} and not complex_pdb.exists():
                protein = Path(args.prepared_proteins) / f"{record['target']}_prepared.pdb"
                minimize_complex(protein, record["ligand_sdf"], record["ligand_name"],
                                 minimized_dir)
            if args.stage in {"fmo", "all"} and complex_pdb.exists():
                job_dir.mkdir(parents=True, exist_ok=True)
                params_link = job_dir / "params"
                if not params_link.exists():
                    params_link.symlink_to(params, target_is_directory=True)
                run_dir = job_dir / "run"
                result_path = run_dir / "result.json"
                result = _existing_success(
                    result_path, job_dir, minimized_dir, args.explicit_water_policy
                )
                if result is None:
                    prepare_job(
                        complex_pdb, record["ligand_sdf"], record["ligand_name"], job_dir,
                        dftb_parameter_dir="../params",
                        minimized_complex=complex_pdb,
                        include_waters=not args.drop_waters,
                    )
                    if run_dir.exists():
                        attempts = job_dir / "attempts"
                        attempts.mkdir(parents=True, exist_ok=True)
                        archive = attempts / time.strftime("%Y%m%dT%H%M%S")
                        archive_index = 2
                        while archive.exists():
                            archive = attempts / f"{time.strftime('%Y%m%dT%H%M%S')}-{archive_index}"
                            archive_index += 1
                        shutil.move(str(run_dir), archive)
                    manifest = json.loads((job_dir / "job_manifest.json").read_text(encoding="utf-8"))
                    result = run_gamess(
                        job_dir / "job.inp", run_dir,
                        gamess_dir=args.gamess_dir, library_dir=args.library_dir,
                        ligand_fragment=manifest["metadata"]["ligand_fragment_id"],
                        timeout_seconds=args.timeout, ncores=1,
                    )
                row.update({
                    "status": result.get("status"),
                    "tpie_kcal_mol": result.get("tpie_kcal_mol"),
                    "wall_seconds": result.get("wall_seconds"),
                    "result_path": str(result_path),
                })
                if result.get("status") != "SUCCEEDED":
                    row["error"] = result.get("error", "")
            elif args.stage == "fmo" and not complex_pdb.exists():
                row.update({"status": "FAILED_MISSING_MINIMIZED_COMPLEX"})
        except Exception as error:
            row.update({"status": f"FAILED_{type(error).__name__}", "error": str(error)})
        if complex_pdb.exists() and (minimized_dir / "minimization.json").exists():
            minimization = json.loads((minimized_dir / "minimization.json").read_text(encoding="utf-8"))
            row.update({
                "ligand_rmsd_angstrom": minimization["ligand_rmsd_angstrom"],
                "explicit_water_residues": minimization.get("explicit_water_residues"),
                "ligand_max_displacement_angstrom": minimization["ligand_max_displacement_angstrom"],
                "minimization_json": str(minimized_dir / "minimization.json"),
            })
        molecule = Chem.MolFromSmiles(record["canonical_isomeric_smiles"])
        if molecule is None:
            raise ValueError(f"invalid cohort SMILES for {key}")
        row["rdkit_crippen_clogp"] = Crippen.MolLogP(molecule)
        rows.append(row)
        print(f"[{index}/{len(records)}] {key}: {row['status']}", flush=True)
        output = root / "fmo_features.csv"
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", encoding="utf-8", newline="") as handle:
            fieldnames = list(rows[0])
            for extra in {column for item in rows for column in item} - set(fieldnames):
                fieldnames.append(extra)
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

    summary = {
        "cohort": "fmo87",
        "n_records": len(rows),
        "n_succeeded": sum(row["status"] == "SUCCEEDED" for row in rows),
        "features": str(root / "fmo_features.csv"),
        "explicit_water_policy": args.explicit_water_policy,
        "protocol": "PDB2PQR/PROPKA + flexible-complex OpenMM/OpenFF ligand restraint + "
                    "full ACE/NME shifted FACIO + FMO2-DFTB3/water PCM PIEDA",
    }
    (root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
