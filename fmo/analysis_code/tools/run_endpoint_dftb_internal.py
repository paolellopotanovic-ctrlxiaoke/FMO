"""Matched gas-phase DFTB3 single points on frozen physical endpoint frames."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import lru_cache
import hashlib
import json
import multiprocessing
from pathlib import Path
import re
import sys
import tempfile
import time

from analyze_endpoint_internal_energy import frame_coordinates, physical_atom_maps, verify_bonds
from analyze_openfe_endpoints import read_topology
from run_matched_openfe import digest, write_json


FIXED_ITERATIONS = (400, 800, 1200, 1600, 2000)


@lru_cache(maxsize=512)
def stat_digest(path, size, modified_ns):
    return digest(path)


def checked_digest(path):
    path = Path(path).resolve(strict=True)
    status = path.stat()
    return stat_digest(str(path), status.st_size, status.st_mtime_ns)


def attach_reuse(plan, previous):
    old = json.loads((previous / "frozen_protocol.json").read_text())
    status = json.loads((previous / "status.json").read_text())
    if status["status"] != "COMPLETE" or old["method"] != plan["method"] or old["labels_used"]:
        raise ValueError("complete matching label-independent source required for reuse")
    identity = lambda task: tuple(task[key] for key in ("target", "leg", "endpoint", "trajectory", "iteration", "replica"))
    tasks = {identity(task): task for task in old["tasks"]}
    rows = {row["directory"]: row for row in status["results"]}
    for task in plan["tasks"]:
        candidate = tasks.get(identity(task))
        if candidate is None:
            continue
        if any(candidate[key] != task[key] for key in ("name", "elements", "xyz_angstrom", "formal_charge", "runtime", "parameters")):
            raise ValueError("same saved frame identity has changed inputs")
        row = rows[candidate["directory"]]
        task["reuse_candidate"] = {"directory": row.get("source_directory", candidate["directory"]),
                                    "short_parameter_path": candidate["short_parameter_path"]}
    plan["reuse_sources"] = {str(previous / name): digest(previous / name) for name in ("frozen_protocol.json", "status.json")}


def verified_cached_result(task, deck):
    from sophosqm_baseline.relative_total_energy import parse_full_energy

    candidate = task.get("reuse_candidate")
    if candidate is None:
        return None
    directory = Path(candidate["directory"])
    old_deck = deck.replace(task["short_parameter_path"], candidate["short_parameter_path"])
    try:
        source = directory / "result.json"
        record = json.loads(source.read_text())
        if record["status"] != "SUCCEEDED" or record["input_sha256"] != hashlib.sha256(old_deck.encode()).hexdigest():
            return None
        if record["executable_sha256"] != checked_digest(Path(task["runtime"]["gamess_dir"]) / "gamess.00.x"):
            return None
        parameter_hash = hashlib.sha256()
        paths = [line.split()[-1] for line in old_deck.splitlines()
                 if re.match(r"\s+[A-Za-z]+\s+[A-Za-z]+\s+\S+\.skf\s*$", line)]
        if not paths:
            return None
        for path in paths:
            parameter_hash.update(path.encode())
            parameter_hash.update(bytes.fromhex(checked_digest(directory / path)))
        if record["parameters_sha256"] != parameter_hash.hexdigest():
            return None
        log = directory / "gamess.log"
        if parse_full_energy(log.read_text(errors="replace"), input_text=old_deck) != record["energy_hartree"]:
            return None
        return {**record, "cache_reused": True, "source_directory": str(directory),
                "source_result_sha256": digest(source), "source_log_sha256": digest(log),
                "wall_seconds_this_invocation": 0.0}
    except (OSError, KeyError, ValueError):
        return None


def freeze_tasks(run, diagnostics_path, workspace, output, stride=400, *, full_density=False):
    from gufe import LigandAtomMapping
    from netCDF4 import Dataset
    from rdkit import Chem

    diagnostics = json.loads(diagnostics_path.read_text())
    if diagnostics["affinity_used_for_frame_selection"] or diagnostics["frozen_protocol_sha256"] != digest(run / "frozen_protocol.json"):
        raise ValueError("matching label-independent frozen diagnostics required")
    registry = workspace / "FMOtpiefinal_HIV1lab/runs/graph_sidechain_transfer_20261021/frozen_protocol.json"
    cohorts = json.loads(registry.read_text())["cohorts"]
    if stride not in (40, 400):
        raise ValueError("only registered 0.1 ns or 1 ns fixed schedules supported")
    iterations = tuple(range(400, 2001, stride)) if full_density else FIXED_ITERATIONS
    alias_file = output / "parameter_aliases.json"
    if alias_file.exists():
        aliases = json.loads(alias_file.read_text())
    else:
        alias_root = Path(tempfile.mkdtemp(prefix="fmo-endpoint-sk-"))
        aliases = {}
        for target in ("BACE", "CDK8"):
            link = alias_root / target
            link.symlink_to(Path(cohorts[target]["parameters"]).resolve(), target_is_directory=True)
            aliases[target] = str(link)
        write_json(alias_file, aliases)
    for target, alias in aliases.items():
        if Path(alias).resolve() != Path(cohorts[target]["parameters"]).resolve() or len(alias) > 100:
            raise ValueError("short parameter alias does not match frozen parameter set")
    tasks = []
    for target in ("BACE", "CDK8"):
        for leg in ("solution", "bound"):
            root = run / "production" / target / "edge_000" / leg
            mapping = LigandAtomMapping.from_json(root / "mapping.json")
            topology_path = sorted((root / "shared").glob("*/hybrid_system.pdb"))[0]
            indices, molecules = physical_atom_maps(read_topology(topology_path), mapping)
            bonds = {endpoint: verify_bonds(topology_path, selected, molecules[endpoint]) for endpoint, selected in indices.items()}
            for repeat in diagnostics["targets"][target][leg]["complete_repeats"]:
                endpoint = repeat["endpoint"]
                frames = [frame for frame in repeat["frames"] if frame["iteration"] in iterations]
                if tuple(frame["iteration"] for frame in frames) != iterations:
                    raise ValueError("all fixed frame times required")
                molecule = molecules[endpoint]
                name = diagnostics["targets"][target]["endpoint_names"][str(endpoint)]
                with Dataset(repeat["trajectory"]) as data:
                    for frame in frames:
                        coordinates = frame_coordinates(data, frame, endpoint, indices[endpoint], bonds[endpoint])
                        path = output / target / leg / Path(repeat["trajectory"]).parent.name / f'endpoint_{endpoint:02d}' / f'iteration_{frame["iteration"]:04d}'
                        tasks.append({"target": target, "leg": leg, "endpoint": endpoint, "name": name,
                            "iteration": frame["iteration"], "replica": frame["replica"],
                            "trajectory": repeat["trajectory"], "directory": str(path),
                            "elements": [atom.GetSymbol() for atom in molecule.GetAtoms()],
                            "xyz_angstrom": coordinates.tolist(), "formal_charge": Chem.GetFormalCharge(molecule),
                            "runtime": cohorts[target]["runtime"][name], "parameters": cohorts[target]["parameters"],
                            "short_parameter_path": aliases[target]})
    return {"schema": "frozen-endpoint-dftb-internal-v1", "tasks": tasks,
            "diagnostics_sha256": digest(diagnostics_path), "driver_sha256": digest(__file__),
            "registry_sha256": digest(registry), "labels_used": False,
            "fixed_iterations": list(iterations),
            "frame_density": "0.1ns" if full_density else "1ns",
            "method": "isolated complete ligand FMO1-DFTB3 gas phase, DC=.TRUE. IDCVER=4",
            "selection": f"All repeats from frozen diagnostics, production 1-5 ns every {stride * .0025:g} ns; no new MD.",
            "not_binding_free_energy": True, "main_pipeline_changed": False}


def run_task(payload):
    task, workspace, timeout = payload
    sys.path.insert(0, str(Path(workspace) / "src"))
    from fmotpiefinal.core.fmo_input import FMODeckOptions, render_fmo2_pieda
    from fmotpiefinal.core.fragmentation import Fragment
    from fmotpiefinal.core.structure import Atom, Ligand, Protein
    from sophosqm_baseline.relative_total_energy import execute_deck

    directory = Path(task["directory"])
    directory.mkdir(parents=True, exist_ok=True)
    atoms = tuple(Atom(index, symbol.upper(), tuple(xyz), record="HETATM")
                  for index, (symbol, xyz) in enumerate(zip(task["elements"], task["xyz_angstrom"], strict=True)))
    ligand = Ligand(task["name"], atoms, formal_charge=task["formal_charge"])
    fragment = Fragment(1, "LIGAND", "ligand", tuple(range(len(atoms))), charge=ligand.formal_charge)
    parameter_link = directory / "params"
    if not parameter_link.exists():
        parameter_link.symlink_to(Path(task["parameters"]).resolve(), target_is_directory=True)
    if parameter_link.resolve() != Path(task["parameters"]).resolve():
        raise ValueError("parameter link changed")
    deck = render_fmo2_pieda(Protein(()), ligand, [fragment], dftb_parameter_dir=task["short_parameter_path"],
                            full_energy=True, deck_options=FMODeckOptions(use_pcm=False))
    deck += " $DFT DC=.TRUE. IDCVER=4 $END\n"
    result = verified_cached_result(task, deck)
    if result is None:
        result = execute_deck(deck, directory, gamess_dir=task["runtime"]["gamess_dir"],
                              library_dir=task["runtime"]["library_dir"], ncores=1, timeout=timeout)
        result.update(cache_reused=False, wall_seconds_this_invocation=result.get("wall_seconds", 0.0))
    return {**{key: task[key] for key in ("target", "leg", "endpoint", "name", "iteration", "replica", "directory")}, **result}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--diagnostics", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--stride", type=int, choices=(40, 400), default=400)
    parser.add_argument("--reuse", type=Path)
    parser.add_argument("--full-density", action="store_true",
                        help="evaluate all 41 registered 0.1 ns frames, instead of five 1 ns frames")
    args = parser.parse_args()
    if args.workers < 1 or args.timeout <= 0:
        parser.error("positive workers and timeout required")
    started = time.monotonic()
    output = args.output.resolve()
    plan = freeze_tasks(args.run.resolve(), args.diagnostics.resolve(), args.workspace.resolve(), output,
                        args.stride, full_density=args.full_density)
    if args.reuse:
        attach_reuse(plan, args.reuse.resolve())
    frozen = output / "frozen_protocol.json"
    if frozen.exists() and json.loads(frozen.read_text()) != plan:
        raise ValueError("existing frozen experiment differs; use new output")
    write_json(frozen, plan)
    write_json(output / "status.json", {"status": "RUNNING", "n_requested": len(plan["tasks"]), "workers": args.workers})
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context("spawn")) as pool:
        futures = {pool.submit(run_task, (task, str(args.workspace.resolve()), args.timeout)): task for task in plan["tasks"]}
        for future in as_completed(futures):
            task = futures[future]
            try:
                row = future.result()
            except Exception as error:
                row = {**{key: task[key] for key in ("target", "leg", "endpoint", "iteration", "directory")},
                       "status": "FAILED", "error": f"{type(error).__name__}: {error}"}
            rows.append(row)
            write_json(output / "progress.json", rows)
            print(len(rows), len(plan["tasks"]), row["target"], row["status"], flush=True)
    succeeded = sum(row["status"] == "SUCCEEDED" for row in rows)
    final = {"status": "COMPLETE" if succeeded == len(rows) else "INCOMPLETE_WITH_FAILURES",
             "n_requested": len(plan["tasks"]), "n_succeeded": succeeded, "wall_seconds": time.monotonic() - started,
             "sum_singlepoint_wall_seconds": sum(row.get("wall_seconds", 0) for row in rows),
             "n_reused": sum(row.get("cache_reused", False) for row in rows),
             "sum_singlepoint_wall_seconds_this_invocation": sum(row.get("wall_seconds_this_invocation", 0) for row in rows),
             "results": rows, "binding_free_energy_claim_allowed": False}
    write_json(output / "status.json", final)
    if succeeded != len(plan["tasks"]):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
