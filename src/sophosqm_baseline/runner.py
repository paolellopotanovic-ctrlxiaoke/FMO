"""Execute a configured GAMESS installation and persist auditable results."""
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import json
import math
import os
import re
import shutil
import subprocess
import time
import uuid

from .gamess import sha256
from .parser import aggregate_tpie, find_selected_tie, parse_pie


def classify_gamess_log(returncode: int, text: str) -> str:
    """Classify a GAMESS run before applying the PIEDA quality gate."""
    if "SCF IS UNCONVERGED" in text or "DIMER CALCULATION DID NOT CONVERGE" in text:
        return "FAILED_SCF_CONVERGENCE"
    if "TERMINATED -ABNORMALLY-" in text or returncode != 0:
        return "FAILED_GAMESS_RUNTIME"
    if "EXECUTION OF GAMESS TERMINATED NORMALLY" not in text:
        return "FAILED_GAMESS_RUNTIME"
    return "SUCCEEDED"


def run_gamess(input_path: str | Path, output_dir: str | Path, *,
               gamess_dir: str | Path, library_dir: str | Path,
               ligand_fragment: int, ncores: int = 1,
               timeout_seconds: float = 3600,
               max_abs_pair_kcal_mol: float = 1000.0) -> dict[str, object]:
    """Run one energy deck; exit status alone is not a GAMESS success check.

    The installation's rungms must already point scratch/restart paths inside
    the project. A unique job name avoids collisions in those shared paths.
    """
    if ncores < 1 or timeout_seconds <= 0 or ligand_fragment < 1 or max_abs_pair_kcal_mol <= 0:
        raise ValueError("ncores, timeout and ligand fragment must be positive")
    source = Path(input_path).resolve(strict=True)
    installation = Path(gamess_dir).resolve(strict=True)
    libraries = Path(library_dir).resolve(strict=True)
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=False)
    job_name = "fmo_" + uuid.uuid4().hex[:12]
    staged = output / f"{job_name}.inp"
    shutil.copyfile(source, staged)
    for directory in ("tmp", "home"):
        (output / directory).mkdir()
    bundled_parameters = installation / "auxdata" / "DFTB" / "3OB-3-1"
    if bundled_parameters.is_dir():
        (output / "tmp" / "params").symlink_to(bundled_parameters, target_is_directory=True)
    environment = os.environ.copy()
    environment.update({
        "LD_LIBRARY_PATH": f"{libraries}:{environment.get('LD_LIBRARY_PATH', '')}",
        "OPAL_PREFIX": str(installation / "runtime" / "mpi"),
        "PATH": f"{installation / 'runtime' / 'mpi' / 'bin'}:{environment.get('PATH', '')}",
        "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
        "OMPI_MCA_btl": "self,tcp", "OMPI_MCA_pml": "ob1",
        "TMPDIR": str(output / "tmp"), "HOME": str(output / "home"),
    })
    if os.geteuid() == 0:
        environment.update(OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
    command = ["/bin/tcsh", str(installation / "rungms"), job_name, "00", str(ncores)]
    result: dict[str, object] = {
        "status": "RUNNING", "command": command, "input_sha256": sha256(source),
        "executable_sha256": sha256(installation / "gamess.00.x"),
        "gamess_dir": str(installation), "library_dir": str(libraries),
        "ligand_fragment": ligand_fragment, "ncores": ncores,
    }
    result_path = output / "result.json"
    result_path.write_text(json.dumps(result, indent=2) + "\n")
    start = time.monotonic()
    log_path = output / "gamess.log"
    try:
        with log_path.open("w") as log:
            process = subprocess.Popen(command, cwd=output, env=environment,
                                       stdout=log, stderr=subprocess.STDOUT,
                                       start_new_session=True)
            try:
                returncode = process.wait(timeout=timeout_seconds)
            except subprocess.TimeoutExpired:
                import signal
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
                raise TimeoutError("GAMESS exceeded its wall-time limit")
        result["returncode"] = returncode
        text = log_path.read_text(errors="replace")
        log_status = classify_gamess_log(returncode, text)
        if log_status != "SUCCEEDED":
            result["status"] = log_status
            message = ("one or more FMO monomer/dimer SCF calculations did not converge"
                       if log_status == "FAILED_SCF_CONVERGENCE"
                       else "GAMESS did not terminate normally; inspect gamess.log")
            raise RuntimeError(message)
        records = parse_pie(log_path)
        outliers = [record for record in records
                    if abs(record.value_kcal_mol) > max_abs_pair_kcal_mol]
        if outliers:
            result.update(status="FAILED_UNPHYSICAL_PIEDA",
                          pair_quality_max_abs_kcal_mol=max_abs_pair_kcal_mol,
                          outlier_pairs=[asdict(record) for record in outliers[:20]])
            raise RuntimeError(
                f"{len(outliers)} PIEDA pairs exceed {max_abs_pair_kcal_mol:g} kcal/mol; "
                "structure/fragment preparation is not physically valid"
            )
        tpie = aggregate_tpie(records, ligand_fragment)
        selected_tie = find_selected_tie(log_path, ligand_fragment)
        if selected_tie is None:
            result.update(status="FAILED_MISSING_SELECTED_TIE")
            raise RuntimeError(f"GAMESS did not print selected TIE for fragment {ligand_fragment}")
        if not math.isclose(tpie, selected_tie, rel_tol=0.0, abs_tol=0.005):
            result.update(status="FAILED_TIE_RECONCILIATION", selected_tie_kcal_mol=selected_tie)
            raise RuntimeError(
                f"pair sum {tpie:.6f} does not reconcile with selected TIE {selected_tie:.6f}"
            )
        result.update(status="SUCCEEDED", pairs=[asdict(record) for record in records],
                      tpie_kcal_mol=tpie, selected_tie_kcal_mol=selected_tie)
        energies = re.findall(r"Total energy of the molecule:\s*(\S+\s*\(2\))=\s*([-+\d.]+)", text)
        result["fmo2_energies_hartree"] = {label: float(value) for label, value in energies}
    except Exception as error:
        if result.get("status") == "RUNNING":
            result["status"] = "FAILED"
        result["error"] = f"{type(error).__name__}: {error}"
    finally:
        result["wall_seconds"] = time.monotonic() - start
        result_path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    return result
