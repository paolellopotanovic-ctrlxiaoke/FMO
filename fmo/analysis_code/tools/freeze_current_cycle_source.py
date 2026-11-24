"""Snapshot accepted main-pipeline response decks for matched cycle regression."""
from __future__ import annotations

import argparse
import io
import json
from pathlib import Path

from fmotpiefinal.core.fmo_convergence import final_fmo_convergence
from fmotpiefinal.core.parser import parse_pie
from run_cdk8_solution_cycle import digest, write_json


def require_accepted_response(decision):
    if (decision.get("action") != "ACCEPT_LOCAL_REPAIR" or not decision.get("accepted")
            or decision.get("affinity_labels_used") is not False
            or decision.get("scores_used") is not False
            or not decision.get("checks") or not all(decision["checks"].values())):
        raise ValueError("verified label-free production acceptance required")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--response-root", type=Path, required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    frozen = json.loads(args.source.read_text())
    cohort = frozen["cohorts"][args.target]
    root = args.response_root.resolve()/args.target
    decision_path = root/"sidechain_response/decision.json"
    require_accepted_response(json.loads(decision_path.read_text()))
    predictions = {row["ligand"]: row for row in json.loads((root/"predictions.json").read_text())}
    if set(predictions) != set(cohort["names"]):
        raise ValueError("complete accepted cohort required")
    for path, fingerprint in cohort["input_sha256"].items():
        if digest(path) != fingerprint:
            raise ValueError("original registry evidence changed")
    for path in (decision_path, root/"predictions.json"):
        cohort["input_sha256"][str(path)] = digest(path)
    for name in cohort["names"]:
        source = root/"sidechain_response/qm"/name/"input.inp"
        log = source.parent/"run/gamess.log"
        runtime_path = source.parent/"run/result.json"
        geometry = root/"sidechain_response/minimized"/name/"minimized_complex.pdb"
        runtime = json.loads(runtime_path.read_text())
        prediction = predictions[name]
        if (runtime["status"] != "SUCCEEDED" or prediction["status"] != "SUCCEEDED"
                or digest(source) != prediction["input_sha256"] or digest(log) != prediction["log_sha256"]):
            raise ValueError("accepted response evidence mismatch")
        final_fmo_convergence(source.read_text(), io.StringIO(log.read_text()))
        pairs = [pair for pair in parse_pie(log) if runtime["ligand_fragment"] in (pair.fragment_i, pair.fragment_j)]
        if (len(pairs) != runtime["ligand_fragment"]-1
                or abs(sum(pair.value_kcal_mol for pair in pairs)-prediction["tpie"]) > 1.e-8):
            raise ValueError("accepted response PIEDA mismatch")
        cohort["sources"][name] = str(source)
        cohort["complexes"][name] = str(geometry)
        cohort["runtime"][name] = runtime
        cohort["baseline"][name] = prediction
        for path in (source, log, runtime_path, geometry):
            cohort["input_sha256"][str(path)] = digest(path)
    frozen["response_source"] = {"original_registry_sha256": digest(args.source),
        "accepted_target": args.target, "decision_sha256": digest(decision_path),
        "driver_sha256": digest(__file__), "affinity_labels_used": False,
        "warning": "Prepared-state snapshot, not a raw structure rerun"}
    if args.output.exists() and json.loads(args.output.read_text()) != frozen:
        raise ValueError("refusing to overwrite a different source snapshot")
    write_json(args.output, frozen)


if __name__ == "__main__":
    main()
