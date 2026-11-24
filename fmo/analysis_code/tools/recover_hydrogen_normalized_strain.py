"""Label-free QA rescue for unconverged constrained RFO optimizations only."""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

from report_local_sidechain_transfer import load_truth
from run_cdk8_solution_cycle import digest, write_json
from run_hydrogen_normalized_strain import evaluate, hydrogen_only_deck, worker
from run_local_sidechain_transfer import directory_key


def qa_recovery_deck(deck):
    if "METHOD=RFO NSTEP=500 OPTTOL=0.000100 TRMIN=0.001" not in deck:
        raise ValueError("unsupported original constrained optimizer profile")
    return deck.replace("METHOD=RFO NSTEP=500 OPTTOL=0.000100 TRMIN=0.001",
                        "METHOD=QA NSTEP=500 OPTTOL=0.000100 TRMIN=0.000010", 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--prior", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    args = parser.parse_args()
    source = args.run.resolve()
    frozen = copy.deepcopy(json.loads((source/"frozen_protocol.json").read_text()))
    original = json.loads((source/"evaluation.json").read_text())
    prior = json.loads(args.prior.read_text())
    if digest(args.prior) != frozen["evaluation_sha256"]:
        raise ValueError("source evaluation mismatch")
    names = [name for name, row in original["rows"].items() if row["status"] != "SUCCEEDED"]
    if not names:
        raise ValueError("no failed cases require recovery")
    output = source/"qa_optimizer_recovery"
    for name in names:
        task = frozen["tasks"][name]
        path = Path(task["input"])
        original_summary = json.loads((path.parent/"verified_summary.json").read_text())
        if (original_summary["status"] != "FAILED"
                or "isolated ligand geometry optimization did not converge" not in original_summary.get("error", "")
                or digest(path) != task["input_sha256"]):
            raise ValueError("only optimizer nonconvergence is eligible, not geometry/chemical/QM failures")
        parent = output/frozen["target"]/directory_key(name)
        parent.mkdir(parents=True, exist_ok=True)
        parameters = Path(task["input"]).parent.parent/"params"
        link = parent/"params"
        if not link.exists():
            link.symlink_to(parameters.resolve(), target_is_directory=True)
        new_input = parent/"hydrogen/input.inp"
        new_input.parent.mkdir(parents=True, exist_ok=True)
        deck = qa_recovery_deck(path.read_text())
        if new_input.exists() and new_input.read_text() != deck:
            raise ValueError("refusing changed recovery deck")
        new_input.write_text(deck)
        frozen["tasks"][name] = dict(task, input=str(new_input), input_sha256=digest(new_input))
    frozen["recovery"] = {"schema": "constrained-qa-unchanged-threshold-v1", "names": names,
        "driver_sha256": digest(__file__), "original_protocol_sha256": digest(source/"frozen_protocol.json"),
        "original_evaluation_sha256": digest(source/"evaluation.json"), "labels_used": False,
        "profile": "QA500 TRMIN1e-5 from the same original coordinates; unchanged OPTTOL1e-4, SCF, Hamiltonian, Cartesian frozen atoms and cavity"}
    path = output/"frozen_protocol.json"
    if path.exists() and json.loads(path.read_text()) != frozen:
        raise ValueError("frozen QA rescue changed")
    write_json(path, frozen)
    rows = original["rows"]
    for name in names:
        _, rows[name] = worker((frozen["tasks"][name], str(args.workspace.resolve())))
        print(name, rows[name]["status"], rows[name].get("error", ""), flush=True)
    result = evaluate(rows, frozen, prior, load_truth(args.workspace.resolve(), frozen["target"]))
    result.update(wall_seconds=original["wall_seconds"],
                  timing_warning="wall_seconds is original full run only; recovery task timing is retained per task")
    write_json(output/"evaluation.json", result)
    print(result["status"], result["metrics"], flush=True)


if __name__ == "__main__":
    main()
