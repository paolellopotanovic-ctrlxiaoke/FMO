"""Native subsystem decomposition matched to independently computed PCM states."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
from pathlib import Path
import subprocess
import sys
import uuid

from run_cdk8_solution_cycle import digest, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--native-name", default="native")
    args = parser.parse_args()
    workspace, root = args.workspace.resolve(), args.run.resolve()
    sys.path.insert(0, str(workspace / "scripts"))
    sys.path.insert(0, str(workspace / "src"))
    from validate_fmo_official_subsystem import render_official_subsystem
    from run_fmo87_native_decomposition import parse_native_terms
    from sophosqm_baseline.runner import classify_gamess_log
    frozen = json.loads((root / "frozen_protocol.json").read_text())
    decks = {name: render_official_subsystem((root / f"pcm/{name}/complex/input.inp").read_text())
             for name in frozen["names"]}
    native = root / args.native_name
    native.mkdir(parents=True, exist_ok=True)
    link = native / "params"
    if not link.exists():
        link.symlink_to((root / f"pcm/{frozen['names'][0]}/params").resolve(), target_is_directory=True)
    write_json(native / "plan.json", {"names": frozen["names"], "driver_sha256": digest(__file__),
        "labels_used": False, "arms": ["native_interaction", "native_interaction_plus_desolvation",
        "native_binding", "tpie_plus_native_desolvation"],
        "warning": "Exploratory decomposition after observed poor PCM three-state result, not blind validation"})

    def execute(name, deck):
        directory = native / name
        directory.mkdir(parents=True, exist_ok=True)
        job = "native_" + uuid.uuid4().hex[:12]
        (directory / "input.inp").write_text(deck)
        (directory / f"{job}.inp").write_text(deck)
        reference = json.loads((root / f"pcm/{name}/complex/result.json").read_text())
        environment = os.environ.copy()
        environment.update(LD_LIBRARY_PATH=reference["library_dir"] + ":" + environment.get("LD_LIBRARY_PATH", ""),
            OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OMPI_MCA_btl="self,tcp", OMPI_MCA_pml="ob1")
        if os.geteuid() == 0:
            environment.update(OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
        for subdirectory in ("tmp", "home"):
            (directory / subdirectory).mkdir(exist_ok=True)
        environment.update(TMPDIR=str(directory / "tmp"), HOME=str(directory / "home"))
        command = ["/bin/tcsh", str(Path(reference["gamess_dir"]) / "rungms"), job, "00", "2"]
        with (directory / "gamess.log").open("w") as handle:
            completed = subprocess.run(command, cwd=directory, env=environment, stdout=handle,
                stderr=subprocess.STDOUT, timeout=1800, check=False)
        text = (directory / "gamess.log").read_text(errors="replace")
        result = {"status": classify_gamess_log(completed.returncode, text), "input_sha256": digest(directory / "input.inp")}
        if result["status"] == "SUCCEEDED":
            result["terms"] = parse_native_terms(text)
            energies = [json.loads((root / key / "result.json").read_text())["energy_hartree"]
                for key in (f"pcm/{name}/complex", "pcm/shared/receptor", f"pcm/{name}/ligand")]
            independent = (energies[0] - energies[1] - energies[2]) * 627.5094740631
            result["native_minus_independent"] = result["terms"]["binding_kcal_mol"] - independent
            if abs(result["native_minus_independent"]) > 0.02:
                result["status"] = "FAILED_INDEPENDENT_RECONCILIATION"
        write_json(directory / "result.json", result)
        print(name, result, flush=True)
        return name, result

    results = {}
    with ThreadPoolExecutor(max_workers=4) as executor:
        for future in as_completed([executor.submit(execute, name, deck) for name, deck in decks.items()]):
            name, result = future.result()
            results[name] = result
            write_json(native / "execution.json", results)


if __name__ == "__main__":
    main()
