"""Execute a configured GAMESS installation and persist auditable results."""
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import json
import io
import math
import os
import re
import shutil
import subprocess
import time
import uuid

from .gamess import sha256
from .parser import aggregate_tpie, find_selected_tie, parse_pie
from .fmo_convergence import final_fmo_convergence
from .scf_stall import MonomerSCFBudget


def classify_gamess_log(returncode: int, text: str, *, deck: str | None = None) -> str:
    """Classify a GAMESS run before applying the PIEDA quality gate."""
    if "SCF IS UNCONVERGED" in text or "DIMER CALCULATION DID NOT CONVERGE" in text:
        if (deck is not None and returncode == 0 and
                "DIMER CALCULATION DID NOT CONVERGE" not in text):
            try:
                final_fmo_convergence(deck, io.StringIO(text))
            except (ValueError, IndexError, OverflowError):
                pass
            else:
                return "SUCCEEDED"
        return "FAILED_SCF_CONVERGENCE"
    if "TERMINATED -ABNORMALLY-" in text or returncode != 0:
        return "FAILED_GAMESS_RUNTIME"
    if "EXECUTION OF GAMESS TERMINATED NORMALLY" not in text:
        return "FAILED_GAMESS_RUNTIME"
    return "SUCCEEDED"


def pie_roundoff_tolerance(pair_count: int) -> float:
    if pair_count < 1:
        raise ValueError("PIEDA pair count must be positive")
    return max(0.005, 0.0005 * pair_count + 0.0001)


def charge_aware_pair_outliers(records, deck: str, max_abs_pair_kcal_mol: float):
    """Preserve strong opposite-charge salt bridges; flag unexplained extremes."""
    section = re.search(r"\$FMO\s+(.*?)\$END", deck, re.S | re.I)
    match = re.search(r"\bICHARG\(1\)\s*=\s*(.*?)\bFRGNAM\(1\)",
                      section.group(1), re.S | re.I) if section else None
    nfrag = re.search(r"\bNFRAG\s*=\s*(\d+)", section.group(1), re.I) if section else None
    charges = [int(value) for value in re.findall(r"[-+]?\d+", match.group(1))] if match else []
    if not nfrag or len(charges) != int(nfrag.group(1)):
        raise ValueError("cannot validate interaction energies without full fragment charges")
    outliers, ionic_pairs = [], []
    for record in records:
        if abs(record.value_kcal_mol) <= max_abs_pair_kcal_mol:
            continue
        first, second = charges[record.fragment_i - 1], charges[record.fragment_j - 1]
        if first * second < 0 and record.value_kcal_mol < 0 and (
                abs(record.value_kcal_mol) <= 2 * max_abs_pair_kcal_mol):
            ionic_pairs.append(record)
        else:
            outliers.append(record)
    return outliers, ionic_pairs


def run_gamess(input_path: str | Path, output_dir: str | Path, *,
               gamess_dir: str | Path, library_dir: str | Path,
               ligand_fragment: int, ncores: int = 1,
               timeout_seconds: float = 3600,
               max_abs_pair_kcal_mol: float = 100.0,
               parameter_dir: str | Path | None = None) -> dict[str, object]:
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
    if parameter_dir is not None:
        parameters = Path(parameter_dir).resolve(strict=True)
        link = output.parent / "params"
        if link.is_symlink() and link.resolve(strict=True) != parameters:
            raise ValueError("existing GAMESS parameter link has changed")
        if not link.exists():
            link.symlink_to(parameters, target_is_directory=True)
    job_name = "fmo_" + uuid.uuid4().hex[:12]
    staged = output / f"{job_name}.inp"
    shutil.copyfile(source, staged)
    for directory in ("tmp", "home"):
        (output / directory).mkdir()
    environment = os.environ.copy()
    environment.update({
        "LD_LIBRARY_PATH": f"{libraries}:{environment.get('LD_LIBRARY_PATH', '')}",
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
            with log_path.open(errors="replace") as progress:
                budget = MonomerSCFBudget()
                pending = ""
                while True:
                    try:
                        remaining = timeout_seconds - (time.monotonic() - start)
                        if remaining <= 0:
                            raise subprocess.TimeoutExpired(command, timeout_seconds)
                        returncode = process.wait(timeout=min(1.0, remaining))
                        break
                    except subprocess.TimeoutExpired:
                        import signal
                        if time.monotonic() - start >= timeout_seconds:
                            os.killpg(process.pid, signal.SIGKILL)
                            process.wait()
                            raise TimeoutError("GAMESS exceeded its wall-time limit")
                        pending += progress.read(4 * 1024 * 1024)
                        lines = pending.split("\n")
                        pending = lines.pop()
                        for line in lines:
                            if budget.consume(line) is not None:
                                os.killpg(process.pid, signal.SIGKILL)
                                process.wait()
                                result.update(status="FAILED_SCF_CONVERGENCE",
                                              numerical_budget=budget.evidence)
                                raise RuntimeError("repeated monomer SCF failure budget exhausted; numerical retry required")
        result["returncode"] = returncode
        text = log_path.read_text(errors="replace")
        log_status = classify_gamess_log(returncode, text, deck=source.read_text())
        if log_status != "SUCCEEDED":
            result["status"] = log_status
            message = ("one or more FMO monomer/dimer SCF calculations did not converge"
                       if log_status == "FAILED_SCF_CONVERGENCE"
                       else "GAMESS did not terminate normally; inspect gamess.log")
            raise RuntimeError(message)
        registered_deck = source.read_text()
        if ("$DFTB " in registered_deck and
                re.search(r"\$FMO\s+.*?\bNBODY=2\b", registered_deck, re.S | re.I) and
                re.search(r"\$PCM\s+.*?\bIFMO=-1\b", registered_deck, re.S | re.I)):
            try:
                result["final_fmo_convergence"] = final_fmo_convergence(registered_deck, io.StringIO(text))
            except (ValueError, IndexError, OverflowError) as error:
                result["status"] = "FAILED_SCF_CONVERGENCE"
                result["convergence_proof_error"] = f"{type(error).__name__}: {error}"
                raise RuntimeError("final FMO monomer/dimer convergence proof failed; numerical retry required") from error
        if "SCF IS UNCONVERGED" in text:
            result["final_fmo_convergence"] = final_fmo_convergence(source.read_text(), io.StringIO(text))
            result["recovered_early_monomer_scf_failures"] = True
        records = parse_pie(log_path)
        outliers, ionic_pairs = charge_aware_pair_outliers(
            records, source.read_text(), max_abs_pair_kcal_mol)
        if ionic_pairs:
            result["strong_ionic_pairs"] = [asdict(record) for record in ionic_pairs]
        if outliers:
            result.update(status="FAILED_UNPHYSICAL_PIEDA",
                          pair_quality_max_abs_kcal_mol=max_abs_pair_kcal_mol,
                          outlier_pairs=[asdict(record) for record in outliers[:20]])
            raise RuntimeError(
                f"{len(outliers)} non-ionic PIEDA pairs exceed {max_abs_pair_kcal_mol:g} "
                "kcal/mol (or exceed the ionic ceiling); inspect geometry"
            )
        tpie = aggregate_tpie(records, ligand_fragment)
        selected_tie = find_selected_tie(log_path, ligand_fragment)
        if selected_tie is None:
            result.update(status="FAILED_MISSING_SELECTED_TIE")
            raise RuntimeError(f"GAMESS did not print selected TIE for fragment {ligand_fragment}")
        pie_roundoff = pie_roundoff_tolerance(len(records))
        if not math.isclose(tpie, selected_tie, rel_tol=0.0, abs_tol=pie_roundoff):
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
