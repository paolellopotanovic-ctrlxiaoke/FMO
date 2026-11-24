"""Matched CDK8 gas/PCM three-state experiment; never fit weights to labels."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import io
import json
import math
from pathlib import Path
import re
import sys
import time

from scipy.stats import kendalltau, pearsonr, spearmanr


HARTREE_TO_KCAL = 627.5094740631
ARMS = {
    "tpie_pcm": "Existing matched production PCM tPIE",
    "three_state_gas_bound": "Egas(PL)-Egas(P)-Egas(Lbound)",
    "three_state_pcm_bound": "Epcm(PL)-Epcm(P)-Epcm(Lbound)",
    "three_state_pcm_free_min": "Epcm(PL)-Epcm(P)-Epcm(Lminimum)",
    "tpie_plus_cycle_solvation": "tPIEpcm + (three_state_pcm_bound-three_state_gas_bound); exploratory hybrid, not a closed cycle",
    "tpie_plus_ligand_strain": "tPIEpcm + Epcm(Lbound)-Epcm(Lminimum); exploratory hybrid",
    "tpie_plus_solvation_and_strain": "Both unit-coefficient exploratory corrections; not a free energy",
}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def group(deck, name):
    matches = list(re.finditer(r"\$" + name + r"\s+(.*?)\$END", deck, re.S | re.I))
    if len(matches) != 1:
        raise ValueError(f"expected one {name} group")
    return matches[0]


def array(text, name, following):
    match = re.search(name + r"\(1\)=(.*?)(?=" + following + r"|$)", text, re.S)
    if not match:
        raise ValueError(f"missing {name}")
    return match[1].replace(",", " ").split()


def format_array(name, values):
    lines = [f"      {name}(1)="]
    for value in values:
        token = str(value) + (" " if name == "INDAT" else ", ")
        if len(lines[-1]) + len(token) > 72:
            lines.append("                ")
        lines[-1] += token
    return "\n".join(line.rstrip() for line in lines).rstrip(",")


def split_decks(source, parameters):
    """Preserve exact source atoms, caps, charges, HOPs and fragment membership."""
    if not parameters.is_dir():
        raise ValueError("parameter directory is missing")
    fragment_group = group(source, "FMO")
    content = fragment_group[1]
    count = int(re.search(r"NFRAG=(\d+)", content)[1])
    charges = array(content, "ICHARG", "FRGNAM")
    names = array(content, "FRGNAM", "INDAT")
    membership = list(map(int, array(content, "INDAT", r"(?!)")))
    coordinates = group(source, "FMOXYZ")[1].strip().splitlines()
    if len(charges) != count or len(names) != count or membership[0] != 0:
        raise ValueError("source fragment manifest mismatch")
    fragments, current = [], []
    for index in membership[1:]:
        if index == 0:
            fragments.append(current)
            current = []
        elif index < 0:
            if not current or -index <= current[-1]:
                raise ValueError("unsupported INDAT range")
            current.extend(range(current[-1] + 1, -index + 1))
        else:
            current.append(index)
    if current or len(fragments) != count:
        raise ValueError("incomplete fragment membership")
    if [index for fragment in fragments for index in fragment] != list(range(1, len(coordinates) + 1)):
        raise ValueError("experiment requires contiguous exact source atom ordering")
    ligand_start = fragments[-1][0] - 1
    if "$FMOBND" in source:
        for line in group(source, "FMOBND")[1].splitlines():
            fields = line.split()
            if fields and max(abs(int(fields[0])), int(fields[1])) > ligand_start:
                raise ValueError("covalently linked ligand needs a different cycle")

    def render(selected_charges, selected_names, selected_fragments, xyz, ligand=False):
        offsets = [0]
        atom_offset = 0
        for fragment in selected_fragments:
            size = len(fragment)
            offsets.extend([atom_offset + 1] + ([-(atom_offset + size)] if size > 1 else []) + [0])
            atom_offset += size
        nfrag = len(selected_fragments)
        replacement = (f"$FMO NFRAG={nfrag} NBODY={min(nfrag, 2)} NACUT=0\n"
            + "\n".join(format_array(name, values) for name, values in
                (("ICHARG", selected_charges), ("FRGNAM", selected_names), ("INDAT", offsets)))
            + "\n $END")
        deck = source[:fragment_group.start()] + replacement + source[fragment_group.end():]
        xyz_group = group(deck, "FMOXYZ")
        deck = deck[:xyz_group.start()] + "$FMOXYZ\n" + "\n".join(xyz) + "\n $END" + deck[xyz_group.end():]
        deck = deck.replace("IPIEDA=1", "IPIEDA=0")
        species = {line.split()[0].upper() for line in xyz}
        if "$DFTBSK" in deck:
            parameter_group = group(deck, "DFTBSK")
            selected = [line for line in parameter_group[1].splitlines() if line.strip()
                        and all(symbol.upper() in species for symbol in line.split()[:2])]
            deck = (deck[:parameter_group.start()] + "$DFTBSK\n" + "\n".join(selected)
                    + "\n $END" + deck[parameter_group.end():])
        data_group = group(deck, "DATA")
        data_lines = data_group[1].strip().splitlines()
        data_lines = data_lines[:2] + [line for line in data_lines[2:]
                                     if line.split()[0].split("-")[0].upper() in species]
        deck = deck[:data_group.start()] + "$DATA\n" + "\n".join(data_lines) + "\n $END" + deck[data_group.end():]
        if ligand:
            for name in ("FMOBND", "FMOHYB"):
                deck = re.sub(r"\s*\$" + name + r"\s+.*?\$END", "", deck, flags=re.S | re.I)
        return deck

    return {
        "complex": render(charges, names, fragments, coordinates),
        "receptor": render(charges[:-1], names[:-1], fragments[:-1], coordinates[:ligand_start]),
        "ligand": render(charges[-1:], names[-1:], fragments[-1:], coordinates[ligand_start:], True),
    }


def gas_deck(deck):
    for name in ("PCM", "PCMCAV", "TESCAV"):
        deck = re.sub(r"\s*\$" + name + r"\s+.*?\$END", "", deck, flags=re.S | re.I)
    return deck


def metrics(truth, predicted):
    if not all(math.isfinite(value) for value in truth + predicted):
        raise ValueError("nonfinite score")
    return {"pearson": float(pearsonr(truth, predicted).statistic),
            "spearman": float(spearmanr(truth, predicted).statistic),
            "kendall": float(kendalltau(truth, predicted).statistic)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--ncores", type=int, default=2)
    parser.add_argument("--pcm-idisp", type=int, choices=(0, 1), default=1)
    parser.add_argument("--skip-gas", action="store_true")
    parser.add_argument("--gas-only", action="store_true")
    args = parser.parse_args()
    if args.workers < 1 or args.ncores < 1:
        parser.error("positive workers and cores required")
    workspace, output = args.workspace.resolve(), args.output.resolve()
    sys.path.insert(0, str(workspace / "src"))
    from sophosqm_baseline.relative_total_energy import (
        execute_deck, build_ligand_optimization_deck, parse_converged_ligand_optimization)
    from fmotpiefinal.core.fmo_convergence import final_fmo_convergence
    from rdkit import Chem

    cohort = workspace / "FMOtpiefinal_HIV1lab/runs/all_accelerated_robust_20261020/CDK8"
    prepared = cohort / "run/qm/full_target_fmo"
    verification = json.loads((cohort / "independent_verification.json").read_text())
    if verification["status"] != "COMPLETE_INDEPENDENTLY_VERIFIED" or verification["n_verified"] != 10:
        raise ValueError("complete verified source required")
    runtime = json.loads((cohort.parent / "frozen_protocol.json").read_text())["jobs"]["CDK8"]["command"]
    gamess = Path(runtime[runtime.index("--gamess-dir") + 1])
    libraries = Path(runtime[runtime.index("--library-dir") + 1])
    parameters = Path(runtime[runtime.index("--dftb-parameter-dir") + 1])
    names = sorted(record["ligand"] for record in verification["records"])
    tasks, receptor_reference, sources = {}, None, {}
    for name in names:
        source = prepared / "qm" / name / "input.inp"
        text = source.read_text()
        record = next(record for record in verification["records"] if record["ligand"] == name)
        if digest(source) != record["input_sha256"]:
            raise ValueError("registered source deck changed")
        sources[str(source)] = digest(source)
        decks = split_decks(text, parameters)
        if args.pcm_idisp == 0:
            decks = {state: deck.replace("IDISP=1", "IDISP=0") for state, deck in decks.items()}
        data = group(decks["receptor"], "DATA")
        lines = data[1].strip().splitlines()
        decks["receptor"] = (decks["receptor"][:data.start()] + "$DATA\nshared receptor\n"
            + "\n".join(lines[1:]) + "\n $END" + decks["receptor"][data.end():])
        if receptor_reference is None:
            receptor_reference = decks["receptor"].split(" $CONTRL", 1)[1]
        elif decks["receptor"].split(" $CONTRL", 1)[1] != receptor_reference:
            raise ValueError("receptor must be identical across complete cohort")
        for state, deck in decks.items():
            identity = "shared" if state == "receptor" else name
            for model, rendered in (("pcm", deck), ("gas", gas_deck(deck))):
                if args.skip_gas and model == "gas":
                    continue
                if args.gas_only and model != "gas":
                    continue
                tasks[f"{model}/{identity}/{state}"] = rendered
        if not args.gas_only:
            tasks[f"optimized/{name}/ligand"] = build_ligand_optimization_deck(
                decks["ligand"], nstep=200, opttol=0.0001).replace(
                    "OPTTOL=0.000100 $END", "OPTTOL=0.000100 METHOD=RFO TRMIN=0.001 $END")
    frozen = {"schema": "cdk8-matched-solution-cycle-v1", "created_utc": datetime.now(timezone.utc).isoformat(),
        "n": 10, "names": names, "n_jobs": len(tasks), "arms": ARMS,
        "source_sha256": sources, "driver_sha256": digest(__file__),
        "deck_sha256": {key: hashlib.sha256(value.encode()).hexdigest() for key, value in tasks.items()},
        "experimental_labels_used_for_generation": False, "target_spearman": 0.6,
        "primary_arm": "three_state_pcm_free_min", "workers": args.workers, "ncores": args.ncores,
        "pcm_idisp": args.pcm_idisp, "skip_gas": args.skip_gas, "gas_only": args.gas_only,
        "scope": "Matched prepared geometry experiment, not a new raw-structure preparation run; no claim of complete binding free energy; fixed receptor and retained structural waters; single local solution minimum; no entropy or state populations",
        "selection_warning": "CDK8 is a development cohort; choosing the best of these arms is not blind validation"}
    frozen_path = output / "frozen_protocol.json"
    if frozen_path.exists():
        old = json.loads(frozen_path.read_text())
        for key in ("source_sha256", "driver_sha256", "deck_sha256", "arms"):
            if old[key] != frozen[key]:
                raise ValueError("frozen protocol changed; choose a new output directory")
    else:
        write_json(frozen_path, frozen)
        (output / "driver_snapshot.py").write_text(Path(__file__).read_text())
    for key in tasks:
        directory = output / key
        directory.parent.mkdir(parents=True, exist_ok=True)
        link = directory.parent / "params"
        if link.is_symlink() and link.resolve() != parameters.resolve():
            raise ValueError("registered parameter directory changed")
        if not link.exists():
            link.symlink_to(parameters.resolve(), target_is_directory=True)
    started = time.monotonic()

    def execute(key, deck):
        directory = output / key
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "input.inp").write_text(deck)
        result = execute_deck(deck, directory, gamess_dir=gamess, library_dir=libraries,
                              ncores=args.ncores, timeout=1800)
        text = (directory / "gamess.log").read_text(errors="replace")
        if key.startswith("optimized/"):
            reference = tasks[key.replace("optimized/", "pcm/")]
            try:
                coordinates, energy = parse_converged_ligand_optimization(text, ligand_deck=reference)
                if result.get("returncode") != 0:
                    raise ValueError("optimizer process did not succeed")
                result.update(status="SUCCEEDED", energy_hartree=energy, optimized_coordinates=coordinates)
                result.pop("error", None)
            except Exception as error:
                result.update(status="FAILED_OPTIMIZATION", error=repr(error))
        elif key.startswith("pcm/") and "NBODY=2" in deck and result["status"] == "SUCCEEDED":
            try:
                result["strict_final_convergence"] = final_fmo_convergence(deck, io.StringIO(text))
            except ValueError as error:
                result.update(status="FAILED_CONVERGENCE_AUDIT", error=repr(error))
        write_json(directory / "result.json", result)
        print(key, result["status"], round(result["wall_seconds"], 2), flush=True)
        return key, result

    results = {}
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [executor.submit(execute, key, deck) for key, deck in tasks.items()]
        for future in as_completed(futures):
            key, result = future.result()
            results[key] = result
            write_json(output / "progress.json", {"completed": len(results), "requested": len(tasks),
                "succeeded": sum(item["status"] == "SUCCEEDED" for item in results.values()),
                "elapsed_seconds": time.monotonic() - started})
    failed = {key: result for key, result in results.items() if result["status"] != "SUCCEEDED"}
    write_json(output / "execution.json", {"status": "INCOMPLETE" if failed else "COMPLETE",
        "failed": failed, "wall_seconds": time.monotonic() - started,
        "sum_job_wall_seconds": sum(result["wall_seconds"] for result in results.values())})
    baseline = {record["ligand"]: record["tpie"] for record in verification["records"]}
    sdf = workspace / "vendor/public_binding_free_energy_benchmark/fep_benchmark_inputs/structure_inputs/misc/cdk8_koehler_ligands.sdf"
    truth = {molecule.GetProp("_Name").strip(): float(molecule.GetProp("r_exp_dg"))
             for molecule in Chem.SDMolSupplier(str(sdf), removeHs=False)}
    if set(truth) != set(names):
        raise ValueError("complete exact label identity required")
    rows = []
    for name in names:
        row = {"ligand": name, "experimental_dg": truth[name], "tpie_pcm": baseline[name]}
        for model in ("gas", "pcm"):
            keys = [f"{model}/{name}/complex", f"{model}/shared/receptor", f"{model}/{name}/ligand"]
            if all(key in results and results[key]["status"] == "SUCCEEDED" for key in keys):
                energies = [results[key]["energy_hartree"] for key in keys]
                row[f"three_state_{model}_bound"] = (energies[0] - energies[1] - energies[2]) * HARTREE_TO_KCAL
        optimized, bound = results.get(f"optimized/{name}/ligand", {}), results.get(f"pcm/{name}/ligand", {})
        if optimized.get("status") == bound.get("status") == "SUCCEEDED":
            row["ligand_strain"] = (bound["energy_hartree"] - optimized["energy_hartree"]) * HARTREE_TO_KCAL
            row["tpie_plus_ligand_strain"] = baseline[name] + row["ligand_strain"]
            if "three_state_pcm_bound" in row:
                row["three_state_pcm_free_min"] = row["three_state_pcm_bound"] + row["ligand_strain"]
        if all(f"three_state_{model}_bound" in row for model in ("gas", "pcm")):
            row["cycle_solvation_increment"] = row["three_state_pcm_bound"] - row["three_state_gas_bound"]
            row["tpie_plus_cycle_solvation"] = baseline[name] + row["cycle_solvation_increment"]
            if "ligand_strain" in row:
                row["tpie_plus_solvation_and_strain"] = row["tpie_plus_cycle_solvation"] + row["ligand_strain"]
        rows.append(row)
    evaluated = {}
    for arm in ARMS:
        if all(arm in row for row in rows):
            evaluated[arm] = metrics([row["experimental_dg"] for row in rows], [row[arm] for row in rows])
    report = {"n": 10, "complete_jobs": not failed, "metrics": evaluated, "rows": rows,
        "labels_sha256": digest(sdf), "target_reached": {arm: values["spearman"] >= 0.6 for arm, values in evaluated.items()},
        "primary_arm": frozen["primary_arm"], "experimental_labels_used_for_fit": False,
        "warning": frozen["scope"] + "; " + frozen["selection_warning"]}
    write_json(output / "evaluation.json", report)
    lines = ["# CDK8 matched solution-cycle experiment", "", "All metrics use the complete ten ligands; lower score and lower experimental DG mean stronger binding.", "", "|Method|Pearson|Spearman|Kendall|", "|---|---:|---:|---:|"]
    lines.extend(f"|{arm}|{values['pearson']:.3f}|{values['spearman']:.3f}|{values['kendall']:.3f}|" for arm, values in evaluated.items())
    lines.extend(["", report["warning"], "", "Incomplete arms are not ranked. No coefficients are fitted. A target threshold is not a guaranteed result."])
    (output / "REPORT.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({"metrics": evaluated, "failed_jobs": list(failed)}, indent=2), flush=True)


if __name__ == "__main__":
    main()
