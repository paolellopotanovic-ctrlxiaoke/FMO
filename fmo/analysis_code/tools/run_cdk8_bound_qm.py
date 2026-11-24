"""Frozen-receptor QM geometry control, not a binding free-energy estimator."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import io
import json
from pathlib import Path
import re
import sys

import numpy as np

from run_cdk8_solution_cycle import digest, group, split_decks, write_json


def optimization_deck(deck, ligand_size, nstep):
    coordinates = group(deck, "FMOXYZ")[1].strip().splitlines()
    if not 0 < ligand_size < len(coordinates) or nstep < 1:
        raise ValueError("invalid active ligand or optimization budget")
    if "RUNTYP=ENERGY" not in deck or re.search(r"\b(?:MODMOL|MODGRD|NZVAR|IACTAT)\s*=", deck):
        raise ValueError("expected an unrestricted full-FMO energy reference")
    first, last = len(coordinates) - ligand_size + 1, len(coordinates)
    result = deck.replace("RUNTYP=ENERGY", "RUNTYP=OPTIMIZE", 1)
    control = group(result, "CONTRL")
    result = result[:control.end() - 4] + " NZVAR=0 $END" + result[control.end():]
    fragment = group(result, "FMO")
    result = result[:fragment.start()] + result[fragment.start():].replace(
        "NACUT=0", "NACUT=0 MODGRD=96 MODESP=0", 1)
    result += (f"\n $STATPT METHOD=RFO NSTEP={nstep} OPTTOL=0.000100\n"
               f" TRMIN=0.001 IACTAT(1)={first},-{last} $END\n")
    if any(len(line) > 80 for line in result.splitlines()):
        raise ValueError("GAMESS input exceeds 80 columns")
    return result


def converged_coordinates(text, reference, ligand_size):
    marker = "***** EQUILIBRIUM GEOMETRY LOCATED *****"
    if marker not in text or "EXECUTION OF GAMESS TERMINATED NORMALLY" not in text:
        raise ValueError("bound QM optimization did not converge and terminate normally")
    if "THE GEOMETRY SEARCH IS NOT CONVERGED" in text:
        raise ValueError("unconverged geometry is not a minimum")
    initial = group(reference, "FMOXYZ")[1].strip().splitlines()
    tail = text.rsplit(marker, 1)[1]
    table = tail.split("COORDINATES OF ALL ATOMS ARE (ANGS)", 1)[1].splitlines()[3:]
    final = []
    for original, row in zip(initial, table, strict=False):
        fields = row.split()
        if len(fields) != 5 or fields[0].upper() != original.split()[0].upper():
            raise ValueError("optimized atom identity or order changed")
        values = [float(value) for value in fields[1:]]
        if not np.isfinite(values).all():
            raise ValueError("nonfinite optimized coordinates")
        final.append(" ".join(fields[:2]) + " " + " ".join(f"{value:.10f}" for value in values[1:]))
    if len(final) != len(initial):
        raise ValueError("optimized coordinates are incomplete")
    initial_fixed = np.array([row.split()[2:] for row in initial[:-ligand_size]], dtype=float)
    final_fixed = np.array([row.split()[2:] for row in final[:-ligand_size]], dtype=float)
    if not np.allclose(initial_fixed, final_fixed, atol=2e-8, rtol=0):
        raise ValueError("supposedly frozen receptor moved")
    return final


def replace_coordinates(deck, coordinates):
    match = group(deck, "FMOXYZ")
    return deck[:match.start()] + "$FMOXYZ\n" + "\n".join(coordinates) + "\n $END" + deck[match.end():]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--names", nargs="+")
    parser.add_argument("--nstep", type=int, default=100)
    parser.add_argument("--timeout", type=int, default=7200)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    if min(args.nstep, args.timeout, args.workers) < 1:
        parser.error("budgets must be positive")
    workspace, reference, output = args.workspace.resolve(), args.reference.resolve(), args.output.resolve()
    sys.path.insert(0, str(workspace / "src"))
    from sophosqm_baseline.relative_total_energy import execute_deck
    from sophosqm_baseline.parser import parse_pie
    from fmotpiefinal.core.fmo_convergence import final_fmo_convergence
    names = json.loads((reference / "frozen_protocol.json").read_text())["names"]
    selected = args.names or names
    if len(set(selected)) != len(selected) or not set(selected) <= set(names):
        parser.error("names must be unique members of the frozen cohort")
    parameters = (reference / f"pcm/{names[0]}/params").resolve(strict=True)
    sources = {}
    for name in selected:
        for state in ("complex", "ligand"):
            path = reference / f"pcm/{name}/{state}/input.inp"
            sources[str(path)] = digest(path)
    frozen = {"driver_sha256": digest(__file__), "sources": sources, "names": selected,
        "complete_cohort": selected == names, "nstep": args.nstep, "timeout": args.timeout,
        "workers": args.workers, "ncores": 2, "labels_used_for_generation": False,
        "gradient": "MODGRD=96: analytic DFTB gradient, no global T/R projection",
        "arms": ["qm_bound_tpie", "qm_bound_three_state", "qm_bound_strain_hybrid"],
        "warning": "Exploratory geometry control; no entropy, no FEP, no fitted weights"}
    output.mkdir(parents=True, exist_ok=True)
    protocol = output / "frozen_protocol.json"
    if protocol.exists() and json.loads(protocol.read_text()) != frozen:
        raise ValueError("frozen protocol changed; use a new output directory")
    write_json(protocol, frozen)
    (output / "driver_snapshot.py").write_text(Path(__file__).read_text())

    def execute(name):
        deck = (reference / f"pcm/{name}/complex/input.inp").read_text()
        ligand = (reference / f"pcm/{name}/ligand/input.inp").read_text()
        ligand_size = len(group(ligand, "FMOXYZ")[1].strip().splitlines())
        runtime = json.loads((reference / f"pcm/{name}/complex/result.json").read_text())
        parent = output / name
        parent.mkdir(exist_ok=True)
        link = parent / "params"
        if not link.exists():
            link.symlink_to(parameters, target_is_directory=True)

        def run(text, stage):
            directory = parent / stage
            directory.mkdir(exist_ok=True)
            (directory / "input.inp").write_text(text)
            return execute_deck(text, directory, gamess_dir=runtime["gamess_dir"],
                library_dir=runtime["library_dir"], ncores=2, timeout=args.timeout)

        optimized = optimization_deck(deck, ligand_size, args.nstep)
        result = run(optimized, "optimization")
        try:
            if result.get("returncode") != 0:
                raise ValueError("optimization process failed")
            coordinates = converged_coordinates((parent / "optimization/gamess.log").read_text(), deck, ligand_size)
            final_deck = replace_coordinates(deck, coordinates)
            states = split_decks(final_deck, parameters)
            energies = {}
            for state in ("complex", "ligand"):
                checked = run(states[state], state)
                if checked["status"] != "SUCCEEDED":
                    raise ValueError(f"final {state} energy verification failed")
                if state == "complex":
                    final_fmo_convergence(states[state], io.StringIO((parent / state / "gamess.log").read_text()))
                energies[state] = checked["energy_hartree"]
            count = int(re.search(r"NFRAG=(\d+)", deck)[1])
            pairs = [pair for pair in parse_pie(parent / "complex/gamess.log")
                     if count in (pair.fragment_i, pair.fragment_j)]
            if len(pairs) != count - 1:
                raise ValueError("incomplete protein-ligand pair inventory")
            result = {"status": "SUCCEEDED", "energies_hartree": energies,
                "tpie_kcal_mol": sum(pair.value_kcal_mol for pair in pairs),
                "optimized_coordinates": coordinates}
        except (ValueError, IndexError) as error:
            result = {"status": "FAILED_GEOMETRY_CONTROL", "error": str(error), "execution": result}
        write_json(parent / "summary.json", result)
        print(name, result["status"], result.get("error", ""), flush=True)
        return name, result

    results = {}
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        for future in as_completed([executor.submit(execute, name) for name in selected]):
            name, result = future.result()
            results[name] = result
            write_json(output / "execution.json", results)


if __name__ == "__main__":
    main()
