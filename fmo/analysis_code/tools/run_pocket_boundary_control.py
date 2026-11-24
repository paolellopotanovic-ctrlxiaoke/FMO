"""Frozen full-cohort QM pocket expansion with unchanged bound poses and waters."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict
import io
import json
from pathlib import Path
import time

from fmotpiefinal.core.bound_geometry import load_bound_ligand, shared_pocket
from fmotpiefinal.core.fmo_convergence import final_fmo_convergence
from fmotpiefinal.core.fmo_input import render_fmo2_pieda
from fmotpiefinal.core.fragmentation import make_facio_like_fragments
from fmotpiefinal.core.parser import parse_pie
from fmotpiefinal.core.runner import run_gamess
from fmotpiefinal.core.structure import Protein, read_pdb, residue_key
from fmotpiefinal.pipeline.target_stages import receptor_at_minimized_geometry
from run_cdk8_solution_cycle import digest, group, metrics, write_json
from run_local_sidechain_transfer import directory_key


def physical_content(deck):
    names = ("CONTRL", "SYSTEM", "SCF", "BASIS", "DFTB", "DFTBSK", "FMOPRP",
             "DFT", "PCM", "PCMCAV", "TESCAV", "FMO", "FMOHYB", "FMOBND", "FMOXYZ")
    return {name: " ".join(group(deck, name)[1].split()) for name in names}


def expansion_keys(protein, ligands, baseline_keys, radius):
    if radius <= 5.0:
        raise ValueError("expanded radius must exceed the registered 5A baseline")
    baseline_keys = set(baseline_keys)
    known = {residue_key(atom) for atom in protein.atoms}
    if not baseline_keys <= known:
        raise ValueError("baseline residues missing from reference receptor")
    expanded = shared_pocket(protein, ligands, radius, include_waters=False)
    keys = baseline_keys | {residue_key(atom) for atom in expanded.atoms}
    waters = {residue_key(atom) for atom in protein.atoms
              if atom.residue_name in {"HOH", "WAT", "TIP"}}
    if keys & waters != baseline_keys & waters:
        raise ValueError("protein-only expansion changed explicit water inventory")
    return keys


def render(receptor, ligand, keys):
    pocket = Protein(tuple(atom for atom in receptor.atoms if residue_key(atom) in keys))
    capped, fragments, metadata = make_facio_like_fragments(
        pocket, ligand, reference_protein=receptor,
        merge_disulfides=True, repair_peptide_numbering=True)
    deck = render_fmo2_pieda(capped, ligand, fragments, dftb_parameter_dir="../params",
                           fragment_bonds=metadata["hop_pairs_gamess_indices"])
    return deck, fragments, metadata


def freeze(source, output, target, radius):
    cohort = json.loads(source.read_text())["cohorts"][target]
    for path, expected in cohort["input_sha256"].items():
        if digest(path) != expected:
            raise ValueError("registered baseline evidence changed")
    protein = read_pdb(cohort["protein"])
    ligands = {name: load_bound_ligand(cohort["sdf"], name, cohort["complexes"][name])
               for name in cohort["names"]}
    keys = expansion_keys(protein, list(ligands.values()), cohort["pocket_keys"], radius)
    jobs = {}
    for name, ligand in ligands.items():
        receptor = receptor_at_minimized_geometry(protein, Path(cohort["complexes"][name]))
        old, _, _ = render(receptor, ligand, set(cohort["pocket_keys"]))
        original = Path(cohort["sources"][name]).read_text()
        if physical_content(old) != physical_content(original):
            raise ValueError(f"5A baseline regeneration failed: {name}")
        deck, fragments, metadata = render(receptor, ligand, keys)
        parent = output / target / directory_key(name)
        path = parent / "input.inp"
        parent.mkdir(parents=True, exist_ok=True)
        if path.exists() and path.read_text() != deck:
            raise ValueError("refusing changed frozen deck")
        path.write_text(deck)
        write_json(parent / "fragment_manifest.json", {
            "fragments": [asdict(fragment) for fragment in fragments], "metadata": metadata})
        jobs[name] = {"input": str(path), "input_sha256": digest(path),
            "complex_sha256": digest(cohort["complexes"][name]),
            "source_input_sha256": digest(cohort["sources"][name]),
            "nfrag": len(fragments), "runtime": cohort["runtime"][name],
            "baseline_tpie": cohort["baseline"][name]["tpie"]}
    frozen = {"schema": "fixed-water-pocket-boundary-v1", "target": target,
        "radius_angstrom": radius, "source_sha256": digest(source), "source": str(source),
        "driver_sha256": digest(__file__), "parameters": cohort["parameters"], "jobs": jobs,
        "baseline_pocket_keys": sorted(cohort["pocket_keys"]), "expanded_keys": sorted(keys),
        "added_keys": sorted(keys - set(cohort["pocket_keys"])),
        "labels_used_for_generation": False, "production_default_changed": False,
        "primary_score": "Full ligand-receptor tPIE including all expanded fragments",
        "unchanged": ["bound ligand geometry", "receptor geometry", "water identities and coordinates",
                      "protonation", "DFTB3/D3/PCM controls", "strict SCF acceptance"],
        "interpretation": "Protein extent, cap boundaries and PCM cavity change together; not a pure polarization ablation"}
    path = output / "frozen_protocol.json"
    if path.exists() and json.loads(path.read_text()) != frozen:
        raise ValueError("refusing changed registered experiment")
    write_json(path, frozen)
    return frozen


def execute(name, job, parameters):
    path = Path(job["input"])
    try:
        if digest(path) != job["input_sha256"]:
            raise ValueError("frozen input changed")
        cached = path.parent / "run/result.json"
        if cached.exists():
            result = json.loads(cached.read_text())
        else:
            result = run_gamess(path, path.parent / "run", gamess_dir=job["runtime"]["gamess_dir"],
                library_dir=job["runtime"]["library_dir"], ligand_fragment=job["nfrag"],
                ncores=1, timeout_seconds=2400, parameter_dir=parameters)
        if result["status"] != "SUCCEEDED" or result["input_sha256"] != job["input_sha256"]:
            raise ValueError("QM failed or cached result identity mismatch")
        log = path.parent / "run/gamess.log"
        proof = final_fmo_convergence(path.read_text(), io.StringIO(log.read_text()))
        pairs = [pair for pair in parse_pie(log) if job["nfrag"] in (pair.fragment_i, pair.fragment_j)]
        inventory = {tuple(sorted((pair.fragment_i, pair.fragment_j))) for pair in pairs}
        expected = {(index, job["nfrag"]) for index in range(1, job["nfrag"])}
        if inventory != expected or len(pairs) != len(expected):
            raise ValueError("incomplete or duplicated ligand-pair inventory")
        row = {"status": "SUCCEEDED", "tpie": sum(pair.value_kcal_mol for pair in pairs),
            "baseline_tpie": job["baseline_tpie"], "input_sha256": digest(path),
            "log_sha256": digest(log), "convergence": proof,
            "pairs": [asdict(pair) for pair in pairs], "wall_seconds": result.get("wall_seconds")}
    except Exception as error:
        row = {"status": "FAILED", "error": f"{type(error).__name__}: {error}"}
    write_json(path.parent / "summary.json", row)
    return name, row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--target", default="BACE")
    parser.add_argument("--radius", type=float, default=7.)
    parser.add_argument("--workers", type=int, default=18)
    parser.add_argument("--register-only", action="store_true")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("positive workers required")
    start = time.time()
    output = args.output.resolve()
    frozen = freeze(args.source.resolve(), output, args.target, args.radius)
    write_json(output / "status.json", {"status": "REGISTERED", "n": len(frozen["jobs"])})
    if args.register_only:
        return
    rows = {}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(execute, name, job, frozen["parameters"])
                   for name, job in frozen["jobs"].items()]
        for future in as_completed(futures):
            name, row = future.result()
            rows[name] = row
            write_json(output / "progress.json", {"n_expected": len(futures), "n_finished": len(rows),
                "n_succeeded": sum(item["status"] == "SUCCEEDED" for item in rows.values()),
                "last": name, "rows": rows})
            print(name, row["status"], row.get("tpie", row.get("error")), flush=True)
    complete = len(rows) == len(frozen["jobs"]) and all(row["status"] == "SUCCEEDED" for row in rows.values())
    report = {"status": "COMPLETE_VERIFIED" if complete else "INCOMPLETE", "rows": rows,
              "wall_seconds": time.time()-start, "production_default_changed": False}
    if complete:
        from analyze_bace_mechanisms import load_truth
        truth = load_truth(args.workspace.resolve(), args.target)
        names = sorted(rows)
        if set(truth) != set(names):
            raise ValueError("full evaluation cohort identity mismatch")
        report["baseline_metrics"] = metrics([truth[name] for name in names],
                                             [-rows[name]["baseline_tpie"] for name in names])
        report["expanded_metrics"] = metrics([truth[name] for name in names],
                                             [-rows[name]["tpie"] for name in names])
    write_json(output / "evaluation.json", report)
    write_json(output / "status.json", {key: value for key, value in report.items() if key != "rows"})


if __name__ == "__main__":
    main()
