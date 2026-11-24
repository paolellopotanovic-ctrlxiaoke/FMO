"""Freeze and supervise the six missing complete-cohort Boltz2 comparators."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import fcntl
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time


EXPECTED = {"CHK1": 13, "CDK8": 10, "BACE": 36, "HSP90": 11, "PTP1B": 23, "FAAH": 24}


def now():
    return datetime.now(timezone.utc).isoformat()


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def write_json(path, data, frozen=False):
    content = json.dumps(data, indent=2, ensure_ascii=False) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    if frozen and path.exists() and path.read_text() != content:
        raise ValueError(f"frozen input changed: {path}")
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(content)
    temporary.replace(path)


def read_csv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def protein_chains(path):
    from Bio.PDB import PDBParser, PPBuilder

    result = []
    for chain in PDBParser(QUIET=True).get_structure("receptor", str(path))[0]:
        peptides = PPBuilder().build_peptides(chain, aa_only=True)
        if not peptides:
            continue
        if len(peptides) != 1 or len(peptides[0]) < 50:
            raise ValueError("fragmented/short source protein requires explicit sequence review")
        sequence = str(peptides[0].get_sequence())
        if not set(sequence) <= set("ACDEFGHIKLMNPQRSTVWY"):
            raise ValueError("unsupported source sequence")
        result.append({"source_chain": chain.id, "id": chr(65 + len(result)), "sequence": sequence})
    if not result or len(result) > 20:
        raise ValueError("source protein chain inventory invalid")
    return result


def prepare(args):
    from rdkit import Chem
    from boltz.main import predict
    import yaml

    root, workspace, cache = args.root.resolve(), args.workspace.resolve(), args.cache.resolve(strict=True)
    source = workspace / "FMOtpiefinal_HIV1lab/runs/all_accelerated_robust_20261020/frozen_protocol.json"
    original = json.loads(source.read_text())
    reference = workspace / "runs/fmo87/boltz_external_predictions/Galectin3/provenance.json"
    previous = json.loads(reference.read_text())
    checkpoints = {name: digest(cache / name) for name in previous["checkpoint_sha256"]}
    if checkpoints != previous["checkpoint_sha256"]:
        raise ValueError("model assets differ from the previous Boltz comparator")
    software = {name: importlib.metadata.version(name) for name in previous["software"]}
    if software != previous["software"]:
        raise ValueError("Boltz runtime versions differ from the previous comparator")
    defaults = {parameter.name: parameter.default for parameter in predict.params}
    if defaults["sampling_steps"] != 200 or defaults["diffusion_samples"] != 1:
        raise ValueError("legacy structure sampling defaults differ from frozen comparator")
    groups = {}
    for target, expected in EXPECTED.items():
        command = original["jobs"][target]["command"]
        protein = Path(command[command.index("--protein") + 1])
        ligand_csv = Path(command[command.index("--ligands") + 1])
        ligands = read_csv(ligand_csv)
        names = [row["ligand_id"] for row in ligands]
        if (len(ligands) != expected or len(set(names)) != expected or
                any(set(row) != {"ligand_id", "smiles"} for row in ligands)):
            raise ValueError(f"label-free full cohort required: {target}")
        chains = protein_chains(protein)
        folder = root / target
        records = []
        for index, ligand in enumerate(ligands, 1):
            molecule = Chem.MolFromSmiles(ligand["smiles"])
            if molecule is None:
                raise ValueError(f"invalid ligand identity: {target}/{ligand['ligand_id']}")
            smiles = Chem.MolToSmiles(molecule, isomericSmiles=True)
            ligand_chain = chr(65 + len(chains))
            document = {"version": 1, "sequences": [
                {"protein": {"id": chain["id"], "sequence": chain["sequence"],
                             "msa": str(folder / "msa" / f"{target}_{chain['id']}.csv")}}
                for chain in chains] + [{"ligand": {"id": ligand_chain, "smiles": smiles}}],
                "properties": [{"affinity": {"binder": ligand_chain}}]}
            path = folder / "inputs" / f"{target}_{index:03d}.yaml"
            content = yaml.safe_dump(document, sort_keys=False)
            path.parent.mkdir(parents=True, exist_ok=True)
            if path.exists() and path.read_text() != content:
                raise ValueError("input YAML changed")
            path.write_text(content)
            records.append({"target": target, "ligand_id": ligand["ligand_id"], "input": str(path),
                            "input_sha256": digest(path), "msa": document["sequences"][0]["protein"]["msa"],
                            "smiles": smiles, "n_heavy_atoms": molecule.GetNumHeavyAtoms(), "ligand_chain": ligand_chain})
        group = {"n": expected, "protein": str(protein), "protein_sha256": digest(protein),
                 "ligands_csv": str(ligand_csv), "ligands_sha256": digest(ligand_csv),
                 "protein_chains": chains, "records": records}
        write_json(folder / "manifest.json", {"cohort": f"{target} complete missing Boltz2 comparator",
                   "n_compounds": expected, "records": records, "labels_used": False}, frozen=True)
        group["manifest_sha256"] = digest(folder / "manifest.json")
        groups[target] = group
    frozen = {"schema": "boltz2-missing-six-frozen-v1", "n_requested": 117, "groups": groups,
              "workspace": str(workspace), "cache": str(cache), "gpu": args.gpu,
              "checkpoint_sha256": checkpoints, "software": software, "reference_sha256": digest(reference),
              "source_fmo_protocol_sha256": digest(source), "runner_sha256": digest(Path(__file__)),
              "comparison_builder_sha256": digest(Path(__file__).with_name("build_model_comparison.py")),
              "legacy_runner": str(workspace / "scripts/run_boltz_fmo87.py"),
              "legacy_runner_sha256": digest(workspace / "scripts/run_boltz_fmo87.py"),
              "model": "boltz2", "seed": 0, "recycling_steps": 3, "sampling_steps": 200,
              "diffusion_samples": 1, "sampling_steps_affinity": 200, "diffusion_samples_affinity": 5,
              "labels_used": False, "source_coordinates_or_templates_supplied": False,
              "all_source_protein_chains_preserved": True, "max_same_protocol_attempts": 2,
              "msa_server_url": "https://api.colabfold.com", "msa_pairing_strategy": "greedy"}
    write_json(root / "frozen_protocol.json", frozen, frozen=True)
    write_json(root / "status.json", {"state": "PREPARED", "n_requested": 117, "n_verified": 0,
               "updated_utc": now(), "groups": {target: "PENDING" for target in EXPECTED}})
    return frozen


def check_frozen(root, frozen):
    if digest(Path(__file__)) != frozen["runner_sha256"] or digest(frozen["legacy_runner"]) != frozen["legacy_runner_sha256"]:
        raise ValueError("frozen execution code changed")
    if digest(Path(__file__).with_name("build_model_comparison.py")) != frozen["comparison_builder_sha256"]:
        raise ValueError("comparison builder changed after input freeze")
    for group in frozen["groups"].values():
        if digest(group["protein"]) != group["protein_sha256"] or digest(group["ligands_csv"]) != group["ligands_sha256"]:
            raise ValueError("source protein/ligand inputs changed")
        for record in group["records"]:
            if digest(record["input"]) != record["input_sha256"]:
                raise ValueError("input YAML changed")
        target = group["records"][0]["target"]
        if digest(root / target / "manifest.json") != group["manifest_sha256"]:
            raise ValueError("manifest changed")


def generate_msa(root, frozen, target):
    from boltz.main import compute_msa

    folder = root / target / "msa"
    folder.mkdir(parents=True, exist_ok=True)
    chains = frozen["groups"][target]["protein_chains"]
    data = {f"{target}_{chain['id']}": chain["sequence"] for chain in chains}
    if not all((folder / f"{name}.csv").exists() for name in data):
        compute_msa(data, target, folder, msa_server_url=frozen["msa_server_url"],
                    msa_pairing_strategy=frozen["msa_pairing_strategy"])
    provenance = {}
    for name, sequence in data.items():
        path = folder / f"{name}.csv"
        rows = read_csv(path)
        if len(rows) < 2 or rows[0]["sequence"] != sequence or not any(row["sequence"] != sequence for row in rows):
            raise ValueError("query-matched non-singleton MSA required; no silent fallback")
        provenance[str(path)] = {"sha256": digest(path), "n_sequences": len(rows), "query": sequence}
    write_json(root / target / "msa_provenance.json", provenance, frozen=True)


def logged_process(command, log, timeout):
    start = time.monotonic()
    with log.open("w") as handle:
        handle.write("$ " + " ".join(command) + "\n")
        handle.flush()
        process = subprocess.Popen(command, stdout=handle, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            returncode = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=20)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            returncode = 124
    return {"returncode": returncode, "wall_seconds": time.monotonic() - start, "log": str(log)}


def audit_target(root, frozen, target):
    from Bio.PDB import MMCIFParser
    import yaml

    group = frozen["groups"][target]
    folder = root / target
    msa = json.loads((folder / "msa_provenance.json").read_text())
    for path, proof in msa.items():
        if digest(path) != proof["sha256"]:
            raise ValueError("MSA fingerprint changed")
    provenance = json.loads((folder / "predictions/provenance.json").read_text())
    if (provenance["checkpoint_sha256"] != frozen["checkpoint_sha256"] or
            provenance["manifest_sha256"] != group["manifest_sha256"] or provenance["software"] != frozen["software"]):
        raise ValueError("inference checkpoints differ")
    for name in ("model", "seed", "recycling_steps", "sampling_steps_affinity", "diffusion_samples_affinity"):
        if provenance[name] != frozen[name]:
            raise ValueError("inference protocol differs")
    rows = read_csv(folder / "predictions/predictions.csv")
    mapped = {row["ligand_id"]: row for row in rows}
    if len(rows) != group["n"] or set(mapped) != {record["ligand_id"] for record in group["records"]}:
        raise ValueError("complete ordered cohort missing from inference outputs")
    cases = []
    for record in group["records"]:
        row = mapped[record["ligand_id"]]
        if not row["status"].startswith("SUCCEEDED") or int(row["returncode"]) != 0:
            raise ValueError(f"inference failed: {record['ligand_id']}")
        path = Path(record["input"])
        document = yaml.safe_load(path.read_text())
        ligand = document["sequences"][-1]["ligand"]
        if (digest(path) != record["input_sha256"] or row["input_sha256"] != record["input_sha256"] or
                ligand["smiles"] != record["smiles"]):
            raise ValueError("YAML identity changed")
        affinity = Path(row["affinity_json"])
        if affinity.name != f"affinity_{path.stem}.json" or not affinity.resolve().is_relative_to((folder / "predictions").resolve()):
            raise ValueError("affinity output identity/path differs")
        value = float(json.loads(affinity.read_text())["affinity_pred_value"])
        if not math.isfinite(value) or not math.isclose(value, float(row["affinity_pred_value"]), abs_tol=1e-12):
            raise ValueError("invalid/mismatched affinity score")
        model_path = affinity.with_name(f"{path.stem}_model_0.cif")
        model = MMCIFParser(QUIET=True).get_structure(target, str(model_path))[0]
        if set(model.child_dict) != {chain["id"] for chain in group["protein_chains"]} | {record["ligand_chain"]}:
            raise ValueError("predicted complex chain inventory differs")
        for chain in group["protein_chains"]:
            if len(model[chain["id"]]) != len(chain["sequence"]):
                raise ValueError("incomplete predicted protein chain")
        if sum(atom.element != "H" for atom in model[record["ligand_chain"]].get_atoms()) != record["n_heavy_atoms"]:
            raise ValueError("predicted ligand heavy atom inventory differs")
        log = folder / "predictions/logs" / f"{path.stem}.log"
        text = log.read_text()
        if "GPU available: True (cuda), used: True" not in text or "Number of failed examples: 0" not in text:
            raise ValueError("GPU execution or success proof missing")
        cases.append({"ligand": record["ligand_id"], "affinity_pred_value": value,
                      "score": -value, "input_sha256": digest(path), "affinity_json_sha256": digest(affinity),
                      "model_cif_sha256": digest(model_path), "log_sha256": digest(log)})
    result = {"status": "COMPLETE_INDEPENDENTLY_VERIFIED", "n": group["n"], "cases": cases,
              "msa_provenance_sha256": digest(folder / "msa_provenance.json"), "labels_read": False}
    write_json(folder / "independent_verification.json", result)
    return result


def supervise(args):
    root = args.root.resolve()
    frozen = json.loads((root / "frozen_protocol.json").read_text())
    with (root / "supervisor.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        check_frozen(root, frozen)
        start = time.monotonic()
        state = {"state": "RUNNING", "start_utc": now(), "pid": os.getpid(), "gpu": frozen["gpu"],
                 "n_requested": 117, "n_verified": 0, "groups": {}, "timings": {}}
        write_json(root / "status.json", state)
        successes, failures = {}, {}
        for target in EXPECTED:
            folder = root / target
            state["groups"][target] = "GENERATING_MSA"
            state["updated_utc"] = now()
            write_json(root / "status.json", state)
            timings = state["timings"][target] = {"msa_attempts": [], "inference_attempts": []}
            for attempt in range(1, 3):
                command = [sys.executable, str(Path(__file__).resolve()), "--root", str(root), "--msa-only", target]
                outcome = logged_process(command, folder / f"msa_attempt{attempt}.log", 1800)
                timings["msa_attempts"].append(outcome)
                if outcome["returncode"] == 0:
                    break
            if outcome["returncode"] != 0:
                failures[target] = "MSA_FAILED_NO_SINGLE_SEQUENCE_FALLBACK"
            else:
                state["groups"][target] = "RUNNING_GPU_INFERENCE"
                write_json(root / "status.json", state)
                for attempt in range(1, 3):
                    command = [sys.executable, frozen["legacy_runner"], "--manifest", str(folder / "manifest.json"),
                        "--out-dir", str(folder / "predictions"), "--python", sys.executable,
                        "--cache", frozen["cache"], "--cuda-visible-devices", str(frozen["gpu"]),
                        "--seed", "0", "--recycling-steps", "3", "--sampling-steps-affinity", "200",
                        "--diffusion-samples-affinity", "5", "--timeout", "3600"]
                    outcome = logged_process(command, folder / f"inference_attempt{attempt}.log", EXPECTED[target] * 3600 + 900)
                    timings["inference_attempts"].append(outcome)
                    predictions_csv = folder / "predictions/predictions.csv"
                    if predictions_csv.exists():
                        shutil.copy2(predictions_csv, folder / f"predictions_attempt{attempt}.csv")
                    try:
                        if outcome["returncode"] != 0:
                            raise ValueError(f"inference driver exit {outcome['returncode']}")
                        successes[target] = audit_target(root, frozen, target)
                        failures.pop(target, None)
                        break
                    except Exception as error:
                        failures[target] = repr(error)
            state["groups"][target] = successes.get(target, {}).get("status", "FAILED_REQUIRES_REVIEW")
            state["n_verified"] = sum(group["n"] for group in successes.values())
            state["failures"] = failures
            state["updated_utc"] = now()
            write_json(root / "status.json", state)
            print(target, state["groups"][target], state["n_verified"], "/117", flush=True)
        check_frozen(root, frozen)
        state.update(end_utc=now(), end_to_end_wall_seconds=time.monotonic() - start)
        if failures:
            state["state"] = "INCOMPLETE_NO_FULL_COMPARISON"
        else:
            for name, expected in frozen["checkpoint_sha256"].items():
                if digest(Path(frozen["cache"]) / name) != expected:
                    raise ValueError("checkpoint changed during inference")
            supplement = {"status": "COMPLETE_INDEPENDENTLY_VERIFIED", "n_verified": 117, "groups": successes,
                          "protocol_sha256": digest(root / "frozen_protocol.json"), "timing": state}
            write_json(root / "verified_results.json", supplement)
            from build_model_comparison import build

            comparison = build(Path(frozen["workspace"]), Path(__file__).resolve().parents[1], root / "verified_results.json")
            if comparison["missing_complete_Boltz2"]:
                raise ValueError("comparison still has missing Boltz2 cohorts")
            state["state"] = "COMPLETE_VERIFIED_COMPARISON_UPDATED"
        write_json(root / "status.json", state)
        print(json.dumps(state, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--workspace", type=Path)
    parser.add_argument("--cache", type=Path)
    parser.add_argument("--gpu", type=int, default=1)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--msa-only", choices=list(EXPECTED))
    args = parser.parse_args()
    if args.prepare_only:
        if args.workspace is None or args.cache is None:
            parser.error("preparation requires --workspace and --cache")
        frozen = prepare(args)
        print(json.dumps({"n_requested": frozen["n_requested"], "targets": list(frozen["groups"]), "gpu": frozen["gpu"]}))
    elif args.msa_only:
        generate_msa(args.root.resolve(), json.loads((args.root / "frozen_protocol.json").read_text()), args.msa_only)
    else:
        try:
            supervise(args)
        except Exception as error:
            path = args.root / "status.json"
            state = json.loads(path.read_text()) if path.exists() else {}
            state.update(state="SUPERVISOR_FAILED", error=repr(error), updated_utc=now())
            write_json(path, state)
            raise
