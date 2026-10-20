#!/usr/bin/env python
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import time


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default="runs/fmo87/boltz_inputs/manifest.json")
    parser.add_argument("--out-dir", default="runs/fmo87/boltz_predictions")
    parser.add_argument("--cache", default="assets/boltz/cache")
    parser.add_argument("--python", default="python")
    parser.add_argument("--accelerator", choices=["gpu", "cpu", "tpu"], default="gpu")
    parser.add_argument("--devices", type=int, default=1)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--recycling-steps", type=int, default=3)
    parser.add_argument("--sampling-steps-affinity", type=int, default=200)
    parser.add_argument("--diffusion-samples-affinity", type=int, default=5)
    parser.add_argument("--timeout", type=int, default=3600)
    parser.add_argument("--only-target")
    parser.add_argument("--cuda-visible-devices")
    args = parser.parse_args()

    manifest = json.loads(Path(args.manifest).read_text(encoding="utf-8"))
    records = manifest["records"]
    if args.only_target:
        records = [record for record in records if record["target"] == args.only_target]
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    boltz = shutil.which("boltz", path=str(Path(args.python).parent))
    if not boltz:
        raise RuntimeError(f"boltz executable not found beside {args.python}")

    rows = []
    def prediction_path(stem: str) -> Path:
        candidates = (
            out_dir / "boltz_results" / stem / "predictions" / stem / f"affinity_{stem}.json",
            out_dir / f"boltz_results_{stem}" / "predictions" / stem / f"affinity_{stem}.json",
            out_dir / "predictions" / stem / f"affinity_{stem}.json",
        )
        matches = [candidate for candidate in candidates if candidate.is_file()]
        if len(set(matches)) > 1:
            raise RuntimeError(f"multiple Boltz affinity outputs for {stem}: {matches}")
        return matches[0] if matches else candidates[-1]

    def wait_for_prediction(path: Path, timeout_seconds: float = 10.0) -> bool:
        deadline = time.monotonic() + timeout_seconds
        while time.monotonic() < deadline:
            if path.is_file():
                return True
            time.sleep(0.2)
        return path.is_file()

    def file_sha256(path: str | Path) -> str:
        digest = hashlib.sha256()
        with Path(path).open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        return digest.hexdigest()

    checkpoint_paths = {
        "boltz2_conf.ckpt": Path(args.cache) / "boltz2_conf.ckpt",
        "boltz2_aff.ckpt": Path(args.cache) / "boltz2_aff.ckpt",
        "ccd.pkl": Path(args.cache) / "ccd.pkl",
    }
    for checkpoint in checkpoint_paths.values():
        if not checkpoint.is_file():
            raise FileNotFoundError(checkpoint)
    provenance = {
        "cohort": manifest["cohort"],
        "model": "boltz2",
        "seed": args.seed,
        "accelerator": args.accelerator,
        "devices": args.devices,
        "cuda_visible_devices": args.cuda_visible_devices or os.environ.get("CUDA_VISIBLE_DEVICES"),
        "recycling_steps": args.recycling_steps,
        "sampling_steps_affinity": args.sampling_steps_affinity,
        "diffusion_samples_affinity": args.diffusion_samples_affinity,
        "checkpoint_sha256": {name: file_sha256(path) for name, path in checkpoint_paths.items()},
        "software": {
            package: importlib.metadata.version(package)
            for package in ("boltz", "torch", "pytorch-lightning", "cuequivariance",
                            "cuequivariance-ops-cu12", "cuequivariance-ops-torch-cu12",
                            "cuequivariance-torch", "triton", "antlr4-python3-runtime")
        },
        "trusted_checkpoint_load": "TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1",
    }
    (out_dir / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    for index, record in enumerate(records, start=1):
        stem = Path(record["input"]).stem
        prediction_json = prediction_path(stem)
        if prediction_json.exists():
            status = "SUCCEEDED_CACHED"
            affinity = json.loads(prediction_json.read_text(encoding="utf-8"))
            returncode = 0
            wall_seconds = 0.0
        else:
            command = [
                boltz, "predict", record["input"],
                "--out_dir", str(out_dir),
                "--cache", args.cache,
                "--model", "boltz2",
                "--accelerator", args.accelerator,
                "--devices", str(args.devices),
                "--seed", str(args.seed),
                "--recycling_steps", str(args.recycling_steps),
                "--sampling_steps_affinity", str(args.sampling_steps_affinity),
                "--diffusion_samples_affinity", str(args.diffusion_samples_affinity),
            ]
            started = time.monotonic()
            environment = os.environ.copy()
            environment["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "1"
            if args.cuda_visible_devices:
                environment["CUDA_VISIBLE_DEVICES"] = args.cuda_visible_devices
            process = subprocess.run(command, text=True, capture_output=True,
                                     timeout=args.timeout, check=False, env=environment)
            wall_seconds = time.monotonic() - started
            returncode = process.returncode
            log_path = out_dir / "logs" / f"{stem}.log"
            log_path.parent.mkdir(parents=True, exist_ok=True)
            log_path.write_text(
                "$ " + " ".join(command) + "\n" +
                process.stdout + "\n--- stderr ---\n" + process.stderr,
                encoding="utf-8",
            )
            if wait_for_prediction(prediction_json):
                status = "SUCCEEDED"
                affinity = json.loads(prediction_json.read_text(encoding="utf-8"))
            else:
                status = "FAILED_BOLTZ_RUNTIME"
                affinity = {}
        rows.append({
            "target": record["target"],
            "ligand_id": record["ligand_id"],
            "input": record["input"],
            "input_sha256": record["input_sha256"],
            "msa": record["msa"],
            "status": status,
            "returncode": returncode,
            "wall_seconds": wall_seconds,
            "affinity_json": str(prediction_json),
            "affinity_pred_value": affinity.get("affinity_pred_value"),
            "affinity_probability_binary": affinity.get("affinity_probability_binary"),
            "experimental_dg_kcal_mol": record["experimental_dg_kcal_mol"],
            "checkpoint_conf_sha256": provenance["checkpoint_sha256"]["boltz2_conf.ckpt"],
            "checkpoint_aff_sha256": provenance["checkpoint_sha256"]["boltz2_aff.ckpt"],
        })
        print(f"[{index}/{len(records)}] {stem}: {status}", flush=True)

    output_csv = out_dir / "predictions.csv"
    if rows:
        import csv
        with output_csv.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    print(json.dumps({"records": len(rows), "succeeded": sum(row["status"].startswith("SUCCEEDED") for row in rows),
                      "output_csv": str(output_csv)}, indent=2))


if __name__ == "__main__":
    main()
