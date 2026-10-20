#!/usr/bin/env python
"""Fail-closed, label-blind alternate QM arms for a single new-target crystal."""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

from .geometry_audit import audit as audit_geometry
from .target_stages import directory_key


ROOT = Path.cwd()
RUNNER = ("-m", "fmotpiefinal.pipeline.target_stages")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def attempts(failure: str) -> list[tuple[str, list[str]]]:
    if failure == "FAILED_GAMESS_RUNTIME":
        return [("heavy_pocket", ["--pocket-heavy-only"]),
                ("heavy_pocket_diis", ["--pocket-heavy-only", "--scf-solver", "diis"])]
    if failure == "FAILED_SCF_CONVERGENCE":
        return [("fmo_soscf_strict_free_guess", ["--scf-solver", "fmo_soscf", "--scc-mixer", "gamess_soscf",
                                      "--soscf-onset", "1.0", "--nguess", "3", "--cnvdmp", "0"]),
                ("fmo_soscf_early", ["--scf-solver", "fmo_soscf", "--scc-mixer", "gamess_soscf",
                                      "--soscf-onset", "1.0"]),
                ("fmo_soscf", ["--scf-solver", "fmo_soscf", "--scc-mixer", "gamess_soscf"]),
                ("diis", ["--scf-solver", "diis"]),
                ("free_k10_diis", ["--repair-clashes", "--scf-solver", "diis"]),
                ("free_k10_soscf", ["--repair-clashes"]),
                ("free_k0p1_soscf", ["--repair-clashes", "--restraint-k", "0.1"]),
                ("fixed_k0p1_soscf", ["--restraint-k", "0.1"]),
                ("heavy_pocket_diis", ["--pocket-heavy-only", "--scf-solver", "diis"]),
                ("free_k10_diis_serial", ["--repair-clashes", "--scf-solver", "diis",
                                          "--minimization-threads", "1"])]
    if failure == "FAILED_POST_K10_GEOMETRY":
        return [("free_k10_soscf", ["--repair-clashes"]),
                ("free_k0p1_soscf", ["--repair-clashes", "--restraint-k", "0.1"])]
    return []


def strict_failure(case: Path, name: str, qm_path: Path | None = None) -> str | None:
    qm = qm_path or case / "qm"
    report = qm / "run_status.json"
    if not report.is_file():
        return None
    status = json.loads(report.read_text())
    if status.get("status") != "REVIEW_INCOMPLETE":
        return None
    result = qm / "full_target_fmo/qm" / directory_key(name) / "run/result.json"
    if result.is_file():
        failure = json.loads(result.read_text()).get("status")
        if attempts(failure):
            return failure
    if "missing or severe post-minimization pose" in status.get("failure", ""):
        return "FAILED_POST_K10_GEOMETRY"
    return None


def minimized_is_safe(output: Path, name: str, target_name: str | None = None) -> bool:
    geometry = audit_geometry(output / f"{(target_name or name).lower()}_prepared.pdb",
                              output / "ligands_no_properties.sdf",
                              minimized_root=output / "minimized")
    if (geometry.get("n_minimized_pose_audits") != 1 or
            geometry["cases"][0]["minimized_pose"]["n_severe_contacts"] != 0):
        return False
    minimization = json.loads((output / "minimized" / name / "minimization.json").read_text())
    return (minimization["ligand_rmsd_angstrom"] <= 2.0 and
            minimization["ligand_max_displacement_angstrom"] <= 4.0)


def rescue(name: str, failure: str, case: Path, library: Path, ncores: int,
           timeout: int, gamess: Path | None = None, params: Path | None = None,
           qm_path: Path | None = None, repair_peptide_numbering: bool = False) -> dict:
    if not attempts(failure):
        return {"structure": name, "status": "REVIEW_REQUIRED", "strict_failure": failure,
                "error": "no registered rescue arm for this failure"}
    qm = qm_path or case / "qm"
    if strict_failure(case, name, qm) != failure:
        return {"structure": name, "status": "REVIEW_REQUIRED", "strict_failure": failure,
                "error": "strict failure status changed; cannot rescue"}
    previous = qm / "full_target_fmo"
    protocol_path = previous / "frozen_protocol.json"
    frozen = json.loads(protocol_path.read_text()) if protocol_path.is_file() else {}
    target_name = frozen.get("target_name", name)
    prepared = previous / f"{target_name.lower()}_prepared.pdb"
    poses = previous / "ligands_no_properties.sdf"
    if not all(path.is_file() for path in (prepared, poses, protocol_path)):
        return {"structure": name, "status": "REVIEW_REQUIRED", "strict_failure": failure,
                "error": "missing frozen strict inputs"}
    strict_numbering = json.loads((previous / "frozen_protocol.json").read_text()).get(
        "repair_peptide_numbering", False)
    if strict_numbering != repair_peptide_numbering:
        raise ValueError("rescue peptide-numbering policy differs from strict QM input")
    existing = case / "rescues/diagnostic_report.json"
    prior = json.loads(existing.read_text()) if existing.exists() else {}
    if prior.get("status") == "RESCUED_ALTERNATE_PROTOCOL":
        if prior.get("strict_failure") != failure:
            raise ValueError("cached rescue was for a different strict failure")
        chosen = case / "rescues" / prior["selected_arm"] / "blind_predictions.json"
        if digest(chosen) != prior["arms"][-1]["prediction_sha256"]:
            raise ValueError("cached rescue prediction changed")
        prior["ranking_eligible"] = True
        prior["target"] = target_name
        prior["n_ligands"] = 1
        existing.write_text(json.dumps(prior, indent=2) + "\n")
        return prior
    if prior and prior.get("strict_failure") != failure:
        raise ValueError("prior rescue was for a different strict failure")
    reports = prior.get("arms", [])
    for arm, extra in attempts(failure):
        if arm in {report["arm"] for report in reports}:
            continue
        output = case / "rescues" / arm
        common = [sys.executable, *RUNNER, "--source", str(qm / "source"),
                  "--set-name", "newtarget", "--target-name", target_name, "--expected", "1",
                  "--output", str(output), "--prepared-protein", str(prepared),
                  "--ligand-sdf", str(poses), "--workers", "1", "--ncores", str(ncores),
                  "--timeout", str(timeout), "--merge-disulfides",
                  *(["--repair-peptide-numbering"] if repair_peptide_numbering else []), *extra]
        if "--repair-clashes" not in extra and "--restraint-k" not in extra:
            common += ["--poses-root", str(previous)]
            stages = ("prepare", "score")
        else:
            stages = ("prepare", "minimize", "score")
        result = {"arm": arm, "parameters": extra, "stages": {}}
        for stage in stages:
            command = [*common, "--stage", stage]
            if stage == "score":
                if stages == ("prepare", "minimize", "score") and not minimized_is_safe(output, name, target_name):
                    result["stages"][stage] = "BLOCKED_GEOMETRY_AUDIT"
                    break
                command += ["--gamess-dir", str(gamess or ROOT / "software/gamess"),
                            "--library-dir", str(library), "--dftb-parameter-dir",
                            str(params or ROOT / "software/gamess/auxdata/DFTB/3OB-3-1")]
            output.mkdir(parents=True, exist_ok=True)
            log = output / f"{stage}.stdout"
            with log.open("w") as handle:
                code = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT,
                                      check=False).returncode
            result["stages"][stage] = code
            if code:
                break
        prediction = output / "blind_predictions.json"
        if result["stages"].get("score") == 0 and prediction.exists():
            value = json.loads(prediction.read_text())
            if (value.get("n_requested") == 1 and value.get("n_succeeded") == 1 and
                    value["predictions"][0].get("status") == "SUCCEEDED" and
                    value["predictions"][0].get("ligand") == name):
                result["tpie"] = float(value["predictions"][0]["tpie"])
                result["prediction_sha256"] = digest(prediction)
        reports.append(result)
        if "tpie" in result:
            break
    report = {"structure": name, "target": target_name, "n_ligands": 1, "strict_failure": failure, "status":
              "RESCUED_ALTERNATE_PROTOCOL" if reports and "tpie" in reports[-1] else "REVIEW_REQUIRED",
              "arms": reports, "selected_arm": reports[-1]["arm"] if reports and
              "tpie" in reports[-1] else None}
    report["ranking_eligible"] = report["status"] == "RESCUED_ALTERNATE_PROTOCOL"
    existing.parent.mkdir(parents=True, exist_ok=True)
    existing.write_text(json.dumps(report, indent=2) + "\n")
    return report
