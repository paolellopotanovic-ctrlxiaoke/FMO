"""Exercise the installed main-pipeline response hook on verified prepared cohorts."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
from types import SimpleNamespace

from fmotpiefinal.core.bound_geometry import load_bound_ligand
from fmotpiefinal.core.fragmentation import make_facio_like_fragments
from fmotpiefinal.core.structure import Protein, read_pdb, residue_key
from fmotpiefinal.pipeline.sidechain_response import apply_response
from report_local_sidechain_transfer import load_truth
from run_cdk8_solution_cycle import digest, metrics, write_json


def execute(payload):
    cohort, output = payload
    names = cohort["names"]
    for path, expected in cohort["input_sha256"].items():
        if digest(path) != expected:
            raise ValueError("verified baseline evidence changed")
    target = cohort["target"]
    original_root = Path(cohort["sources"][names[0]]).parents[2]
    frozen = json.loads((original_root / "frozen_protocol.json").read_text())
    runtime = cohort["runtime"][names[0]]
    defaults = {"merge_disulfides": True, "repair_peptide_numbering": True, "use_pcm": True,
        "dampex": 4., "maxit": 200, "scf_solver": "soscf", "nguess": 2,
        "mconv_dimer": None, "scc_mixer": "broyden", "soscf_onset": .25, "cnvdmp": 50}
    options = {key: frozen.get(key, value) for key, value in defaults.items()}
    options.update(target_name=target, gamess_dir=Path(runtime["gamess_dir"]),
        library_dir=Path(runtime["library_dir"]), dftb_parameter_dir=Path(cohort["parameters"]),
        ncores=2, timeout=1800)
    args = SimpleNamespace(**options)
    protein = read_pdb(cohort["protein"])
    keys = set(cohort["pocket_keys"])
    pocket = Protein(tuple(atom for atom in protein.atoms if residue_key(atom) in keys))
    ligands = [load_bound_ligand(cohort["sdf"], name, cohort["complexes"][name]) for name in names]
    prepared, fragments, metadata = make_facio_like_fragments(pocket, ligands[0],
        include_waters=True, reference_protein=protein, use_hop=True,
        merge_disulfides=args.merge_disulfides, repair_peptide_numbering=args.repair_peptide_numbering)
    rows = [cohort["baseline"][name] for name in names]
    result = apply_response(Path(cohort["protein"]), Path(cohort["sdf"]), cohort["complexes"],
        protein, pocket, ligands, args, Path(output) / target, rows,
        (prepared, fragments[:-1], metadata["hop_pairs_gamess_indices"]))
    write_json(Path(output) / target / "predictions.json", result)
    return target, result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", required=True, type=Path)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    source, output = args.source.resolve(), args.output.resolve()
    cohorts = json.loads((source / "frozen_protocol.json").read_text())["cohorts"]
    write_json(output / "frozen_protocol.json", {"schema": "main-response-prepared-regression-v1",
        "source_sha256": digest(source / "frozen_protocol.json"), "targets": list(cohorts),
        "n": sum(len(cohort["names"]) for cohort in cohorts.values()),
        "scope": "Actual production response hook; reuse verified baseline preparation/QM, not a raw PDB-to-QM rerun",
        "labels_used_for_routing": False})
    results = {}
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures = [pool.submit(execute, (cohort, str(output))) for cohort in cohorts.values()]
        for future in as_completed(futures):
            target, records = future.result()
            truth = load_truth(args.workspace.resolve(), target)
            if set(truth) != {record["ligand"] for record in records}:
                raise ValueError("complete metric cohort required")
            decision = json.loads((output / target / "sidechain_response/decision.json").read_text())
            results[target] = {"n": len(records), "action": decision["action"],
                "metrics": metrics([truth[record["ligand"]] for record in records], [-record["tpie"] for record in records])}
            write_json(output / "results.json", results)
            print(target, results[target], flush=True)
    write_json(output / "final_status.json", {"status": "COMPLETE", "n_targets": len(results),
        "n_cases": sum(record["n"] for record in results.values())})


if __name__ == "__main__":
    main()
