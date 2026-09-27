#!/usr/bin/env python
"""Strict production audit for the frozen Boltz-2 four-target cohort."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import pandas as pd
import yaml


EXPECTED_COUNTS = {"CDK2": 16, "JNK1": 21, "P38": 34, "TYK2": 16}
EXPECTED_COHORT_SHA256 = "19edfb8e2ccf3ab1029d848e4062d70cb35532d2b4fee63674f160b2a292997b"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def audit_record(row: pd.Series, manifest_record: dict, provenance: dict,
                 checkpoint_hashes: dict[str, str]) -> dict:
    key = f"{row['target']}|{row['ligand_id']}"
    result = {"key": key}
    checks: dict[str, bool] = {}

    input_path = Path(row["input"])
    yaml_data = yaml.safe_load(input_path.read_text(encoding="utf-8"))
    sequence = yaml_data["sequences"][0]["protein"]
    ligand = yaml_data["sequences"][1]["ligand"]
    checks.update(
        input_exists=input_path.is_file(),
        input_hash_matches=sha256(input_path) == row["input_sha256"] == manifest_record["input_sha256"],
        manifest_input_matches=str(input_path) == manifest_record["input"],
        yaml_smiles_matches=ligand["smiles"] == manifest_record["smiles"],
        yaml_sequence_hash_matches=hashlib.sha256(
            sequence["sequence"].encode()).hexdigest() == manifest_record["protein_sequence_sha256"],
        msa_exists=Path(row["msa"]).is_file(),
        csv_affinity_json_matches=Path(row["affinity_json"]).is_file()
        and math.isclose(float(row["affinity_pred_value"]),
                         float(json.loads(Path(row["affinity_json"]).read_text())["affinity_pred_value"]),
                         rel_tol=0.0, abs_tol=1e-12),
        csv_probability_json_matches=Path(row["affinity_json"]).is_file()
        and math.isclose(float(row["affinity_probability_binary"]),
                         float(json.loads(Path(row["affinity_json"]).read_text())["affinity_probability_binary"]),
                         rel_tol=0.0, abs_tol=1e-12),
        checkpoint_conf_hash_matches=checkpoint_hashes["boltz2_conf.ckpt"] == provenance["checkpoint_sha256"]["boltz2_conf.ckpt"],
        checkpoint_aff_hash_matches=checkpoint_hashes["boltz2_aff.ckpt"] == provenance["checkpoint_sha256"]["boltz2_aff.ckpt"],
        ccd_hash_matches=checkpoint_hashes["ccd.pkl"] == provenance["checkpoint_sha256"]["ccd.pkl"],
    )
    result["checks"] = checks
    result["passed"] = all(checks.values())
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", default="data/fmo87_cohort.csv")
    parser.add_argument("--manifest", default="runs/fmo87/boltz_inputs/manifest.json")
    parser.add_argument("--root", default="runs/fmo87/boltz_predictions")
    parser.add_argument("--cache", default="assets/boltz/cache")
    parser.add_argument("--output-dir", default="runs/fmo87/final")
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()

    cohort_path = Path(args.cohort)
    manifest = json.loads(Path(args.manifest).read_text(encoding="utf-8"))
    manifest_records = {
        f"{record['target']}|{record['ligand_id']}": record
        for record in manifest["records"]
    }
    cohort = list(csv.DictReader(cohort_path.open(encoding="utf-8", newline="")))
    cohort_keys = {f"{row['target']}|{row['ligand_name']}" for row in cohort}
    records = []
    root = Path(args.root)
    cache = Path(args.cache)
    checkpoint_hashes = {
        name: sha256(cache / name)
        for name in ("boltz2_conf.ckpt", "boltz2_aff.ckpt", "ccd.pkl")
    }
    protocol_ok = True
    for target, expected_count in EXPECTED_COUNTS.items():
        target_dir = root / target
        frame = pd.read_csv(target_dir / "predictions.csv", dtype={"ligand_id": str})
        provenance = json.loads((target_dir / "provenance.json").read_text(encoding="utf-8"))
        target_protocol = {
            "model": provenance.get("model") == "boltz2",
            "seed_0": provenance.get("seed") == 0,
            "recycling_3": provenance.get("recycling_steps") == 3,
            "sampling_200": provenance.get("sampling_steps_affinity") == 200,
            "diffusion_5": provenance.get("diffusion_samples_affinity") == 5,
            "gpu_1": provenance.get("accelerator") == "gpu" and provenance.get("devices") == 1,
            "boltz_2_2_1": provenance.get("software", {}).get("boltz") == "2.2.1",
        }
        protocol_ok &= all(target_protocol.values())
        if len(frame) != expected_count or frame["target"].nunique() != 1:
            raise ValueError(f"invalid Boltz2 count for {target}: {len(frame)} != {expected_count}")
        if frame.duplicated(["target", "ligand_id"]).any():
            raise ValueError(f"duplicate Boltz2 identity in {target}")
        for row in frame.to_dict(orient="records"):
            key = f"{row['target']}|{row['ligand_id']}"
            record = audit_record(row, manifest_records[key], provenance, checkpoint_hashes)
            record["checks"]["status_succeeded"] = str(row["status"]).startswith("SUCCEEDED")
            record["checks"]["cohort_identity_exists"] = key in cohort_keys
            record["passed"] = all(record["checks"].values())
            records.append(record)

    passed = sum(record["passed"] for record in records)
    audit = {
        "cohort": {"n": len(cohort), "sha256": sha256(cohort_path),
                   "sha256_matches_frozen": sha256(cohort_path) == EXPECTED_COHORT_SHA256},
        "manifest": {"n_records": len(manifest_records),
                     "identities_match_cohort": set(manifest_records) == cohort_keys,
                     "defaults_match_production": manifest["defaults"].get("diffusion_samples_affinity") == 5},
        "protocol": {"all_targets_match": protocol_ok},
        "n_compounds": len(records), "n_passed": passed,
        "strict_complete": passed == 87 and protocol_ok and set(manifest_records) == cohort_keys,
        "records": records,
        "protocol_deviations": [
            "Boltz-2 is run zero-shot from canonical ligand SMILES and shared target MSAs, not Glide poses.",
            "The local Boltz-2 run uses seed 0, recycling 3, affinity sampling 200, and affinity diffusion samples 5.",
        ],
    }
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    path = output / "boltz_production_audit.json"
    path.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    failures = [record for record in records if not record["passed"]]
    report = [
        "# Boltz-2 FMO87 production audit", "",
        f"- Cohort SHA-256: `{audit['cohort']['sha256']}`",
        f"- Compounds passed: {passed}/{len(records)}",
        f"- Strict complete: `{audit['strict_complete']}`", "",
    ]
    if failures:
        report.extend(["## Failures", ""])
        report.extend(
            f"- `{record['key']}`: {[key for key, value in record['checks'].items() if not value]}"
            for record in failures
        )
    (output / "BOLTZ_PRODUCTION_AUDIT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(json.dumps({"n_compounds": len(records), "n_passed": passed,
                      "strict_complete": audit["strict_complete"], "audit": str(path)}, indent=2))
    if args.require_complete and not audit["strict_complete"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
