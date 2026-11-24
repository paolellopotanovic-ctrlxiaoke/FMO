"""Use the GAMESS documented DFTB3/3ob PCM cavity recipe without label fitting."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import re
import time

from run_cdk8_solution_cycle import digest, group, metrics, write_json
from run_local_sidechain_transfer import directory_key
from run_pocket_boundary_control import execute


def documented_cavity(deck):
    pcm = group(deck, "PCM")
    content = pcm[1]
    for key in ("ICAV", "IDISP"):
        if len(re.findall(r"\b"+key+r"=\d+", content)) != 1:
            raise ValueError("explicit, unique original solvent controls required")
        content = re.sub(r"\b"+key+r"=\d+", key+"=0", content)
    deck = deck[:pcm.start(1)]+content+deck[pcm.end(1):]
    cavity = group(deck, "PCMCAV")
    if cavity[1].strip() != "RADII=VANDW":
        raise ValueError("unexpected original cavity; do not silently override custom radii")
    return deck[:cavity.start(1)]+"RADII=DFTB33OB ALPHA(1)=1.0 "+deck[cavity.end(1):]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--target", default="BACE")
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("positive workers required")
    start = time.time()
    source, output = args.source.resolve(), args.output.resolve()
    cohort = json.loads(source.read_text())["cohorts"][args.target]
    for path, expected in cohort["input_sha256"].items():
        if digest(path) != expected:
            raise ValueError("registered baseline evidence changed")
    jobs = {}
    for name in cohort["names"]:
        original = Path(cohort["sources"][name]).read_text()
        deck = documented_cavity(original)
        if group(deck, "FMOXYZ")[1] != group(original, "FMOXYZ")[1]:
            raise ValueError("cavity control changed geometry")
        parent = output / args.target / directory_key(name)
        path = parent / "input.inp"
        parent.mkdir(parents=True, exist_ok=True)
        if path.exists() and path.read_text() != deck:
            raise ValueError("refusing changed registered deck")
        path.write_text(deck)
        jobs[name] = {"input": str(path), "input_sha256": digest(path),
            "source_input_sha256": digest(cohort["sources"][name]),
            "nfrag": int(re.search(r"NFRAG=(\d+)", deck)[1]), "runtime": cohort["runtime"][name],
            "baseline_tpie": cohort["baseline"][name]["tpie"]}
    manual = Path(cohort["runtime"][cohort["names"][0]]["gamess_dir"]) / "docs-input.txt"
    frozen = {"schema": "documented-dftb33ob-pcm-cavity-v1", "target": args.target,
        "source_sha256": digest(source), "source": str(source), "jobs": jobs,
        "parameters": cohort["parameters"], "driver_sha256": digest(__file__),
        "executor_sha256": digest(Path(__file__).with_name("run_pocket_boundary_control.py")),
        "manual_sha256": digest(manual), "manual": str(manual),
        "recipe": "GAMESS $PCMCAV RADII=DFTB33OB ALPHA(1)=1.0; $PCM ICAV=0 IDISP=0 IEF=-10",
        "baseline": "VANDW / implicit ALPHA=1.2 / ICAV=1 IDISP=1",
        "labels_used_for_generation": False, "production_default_changed": False,
        "interpretation": "Whole documented cavity/nonpolar recipe; not a single-parameter radius fit or full binding free energy"}
    path = output / "frozen_protocol.json"
    if path.exists() and json.loads(path.read_text()) != frozen:
        raise ValueError("refusing changed frozen experiment")
    write_json(path, frozen)
    write_json(output / "status.json", {"status": "RUNNING", "n": len(jobs)})
    rows = {}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(execute, name, job, frozen["parameters"]) for name, job in jobs.items()]
        for future in as_completed(futures):
            name, row = future.result()
            rows[name] = row
            write_json(output / "progress.json", {"n_expected": len(jobs), "n_finished": len(rows),
                "n_succeeded": sum(item["status"] == "SUCCEEDED" for item in rows.values()), "rows": rows})
            print(name, row["status"], row.get("tpie", row.get("error")), flush=True)
    complete = len(rows) == len(jobs) and all(row["status"] == "SUCCEEDED" for row in rows.values())
    report = {"status": "COMPLETE_VERIFIED" if complete else "INCOMPLETE", "rows": rows,
              "wall_seconds": time.time()-start, "production_default_changed": False}
    if complete:
        from report_local_sidechain_transfer import load_truth
        truth = load_truth(args.workspace.resolve(), args.target)
        if set(truth) != set(rows):
            raise ValueError("full cohort evaluation identity mismatch")
        names = sorted(rows)
        observed = [truth[name] for name in names]
        report["baseline_metrics"] = metrics(observed, [-rows[name]["baseline_tpie"] for name in names])
        report["cavity_metrics"] = metrics(observed, [-rows[name]["tpie"] for name in names])
    write_json(output / "evaluation.json", report)
    write_json(output / "status.json", {key: value for key, value in report.items() if key != "rows"})


if __name__ == "__main__":
    main()
