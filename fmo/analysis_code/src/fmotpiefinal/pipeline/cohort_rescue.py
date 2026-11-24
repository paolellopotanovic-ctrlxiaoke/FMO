#!/usr/bin/env python
"""Label-blind whole-cohort pocket rescue after a frozen strict QM failure."""
from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
from pathlib import Path


import numpy as np
from scipy.spatial import cKDTree

from .geometry_audit import minimized_atoms
from .target_intake import digest
from .single_case_rescue import ROOT, RUNNER
from .target_stages import directory_key
from fmotpiefinal.core.water_hydrogen_transfer import preserve_source_waters, valid_water, waters, xyz


RECOVERABLE = {"FAILED_GAMESS_RUNTIME", "FAILED_SCF_CONVERGENCE"}
NUMERICAL_ARMS = {"whole_cohort_fmo_soscf": 0.25, "whole_cohort_fmo_soscf_early": 1.0,
                  "whole_cohort_fmo_soscf_free_guess": 1.0,
                  "whole_cohort_fmo_soscf_strict_free_guess": 1.0}

FREE_GUESS_ARMS = {"whole_cohort_fmo_soscf_free_guess", "whole_cohort_fmo_soscf_strict_free_guess"}


def numerical_soscf_eligible(failure: dict) -> bool:
    failed = failure.get("failed_cases", {})
    return bool(failed) and all(state == "FAILED_SCF_CONVERGENCE" for state in failed.values())


def preparation_stages_verified(qm: Path, status: dict, prepared: Path) -> bool:
    stages = status.get("stage_completed")
    if stages == ["pdbfixer+geometry", "prepare", "minimize"]:
        return status.get("graph_preparation_sha256") is None
    if stages != ["pdbfixer+geometry", "chemical_graph+repair+reaudit", "prepare", "minimize"]:
        return False
    graph_path = qm / "graph_preparation/graph_preparation.json"
    if not graph_path.is_file() or status.get("graph_preparation_sha256") != digest(graph_path):
        return False
    graph = json.loads(graph_path.read_text())
    return (graph.get("status") in {"REPAIRED_AND_VERIFIED", "UNCHANGED_VERIFIED"}
            and graph.get("output_sha256") == digest(prepared))


def cohort_arm_controls(arm: str, excluded: str | None) -> list[str]:
    if arm in NUMERICAL_ARMS:
        return ["--scf-solver", "fmo_soscf", "--scc-mixer", "gamess_soscf",
                *(["--soscf-onset", "1.0"] if NUMERICAL_ARMS[arm] == 1.0 else []),
                *(["--nguess", "3"] if arm in FREE_GUESS_ARMS else []),
                *(["--cnvdmp", "0"] if arm == "whole_cohort_fmo_soscf_strict_free_guess" else [])]
    if arm == "whole_cohort_ligand_scf_diis":
        return ["--scf-solver", "diis"]
    return ["--exclude-pocket-residue", excluded] if excluded else ["--pocket-heavy-only"]


def failing_fragment(log: Path) -> str | None:
    fragments: dict[int, str] = {}
    current = None
    with log.open(errors="replace") as handle:
        for line in handle:
            fragment = re.match(r"^\s*(\d+)\s+([A-Za-z][A-Za-z]{3}\d+)\s+[-+]?\d+\s+\d+\s+\d+", line)
            if fragment:
                fragments[int(fragment[1])] = fragment[2].upper()
            monomer = re.search(r"\bRHF\s+monomer\s+(\d+)\s+corr=", line)
            if monomer:
                current = int(monomer[1])
            elif re.search(r"\bRHF\s+dimer\b", line):
                current = None
            if "SCF IS UNCONVERGED" in line:
                return fragments.get(current)
    return None


def diagnosed_peripheral_residue(qm: Path, failed: dict, pocket: dict) -> str | None:
    scf_failed = {name for name, state in failed.items() if state == "FAILED_SCF_CONVERGENCE"}
    if len(scf_failed) < 2 or any(state not in RECOVERABLE | {"FAILED_TIMEOUT"}
                                  for state in failed.values()):
        return None
    locations = set()
    for name in scf_failed:
        log = qm / "full_target_fmo/qm" / directory_key(name) / "run/gamess.log"
        if not log.is_file():
            return None
        locations.add(failing_fragment(log))
    if len(locations) != 1 or None in locations:
        return None
    fragment = locations.pop()
    matches = [entry for entry in pocket["membership_differences"]
               if entry["residue"].replace(":", "").upper() == fragment and
               entry.get("selected_all_atom") is True and
               entry.get("selected_heavy_atom") is False and
               entry.get("direct_contact_all_atom") is True and
               entry.get("direct_contact_heavy_atom") is False and
               entry.get("nearest_heavy_atom_angstrom", 0) > pocket["cutoff_angstrom"] and
               entry["residue"].split(":")[1] not in {"HOH", "WAT", "TIP"}]
    return matches[0]["residue"] if len(matches) == 1 else None


def eligible_failure(qm: Path, names: list[str]) -> dict | None:
    path = qm / "run_status.json"
    if not path.is_file():
        return None
    status = json.loads(path.read_text())
    if status.get("status") != "REVIEW_INCOMPLETE" or status.get("names") != names:
        return None
    results = {}
    for name in names:
        result = qm / "full_target_fmo/qm" / directory_key(name) / "run/result.json"
        if result.is_file():
            record = json.loads(result.read_text())
            state = record.get("status")
            if state == "FAILED" and record.get("error") == "TimeoutError: GAMESS exceeded its wall-time limit":
                state = "FAILED_TIMEOUT"
            results[name] = state
    failed = {name: state for name, state in results.items() if state in RECOVERABLE | {"FAILED_TIMEOUT"}}
    if not failed or any(state not in RECOVERABLE | {"FAILED_TIMEOUT", "SUCCEEDED"}
                         for state in results.values()):
        return None
    if not any(token in status.get("failure", "") for token in
               ("CalledProcessError", "GAMESS FMO output does not cover")):
        return None
    return {"status_sha256": digest(path), "failed_cases": failed}


def ligand_scf_failures(qm: Path, failed: dict) -> bool:
    if not failed or any(state != "FAILED_SCF_CONVERGENCE" for state in failed.values()):
        return False
    for name in failed:
        run = qm / "full_target_fmo/qm" / directory_key(name) / "run"
        result, log = run / "result.json", run / "gamess.log"
        if not result.is_file() or not log.is_file():
            return False
        ligand_id = json.loads(result.read_text()).get("ligand_fragment")
        if not isinstance(ligand_id, int):
            return False
        current = None
        found = False
        with log.open(errors="replace") as handle:
            for line in handle:
                monomer = re.search(r"\bRHF\s+monomer\s+(\d+)\s+corr=", line)
                if monomer:
                    current = int(monomer[1])
                elif re.search(r"\bRHF\s+dimer\b", line):
                    current = None
                if "SCF IS UNCONVERGED" in line:
                    if current != ligand_id:
                        return False
                    found = True
        if not found:
            return False
    return True


def source_water_rescue(qm: Path, prepared: Path, residue: str | None,
                        failure: dict) -> tuple[Path, dict] | None:
    source = qm / "input_target_water_treated.pdb"
    if (residue is None or not source.is_file() or
            len(failure["failed_cases"]) < 2 or not prepared.is_file()):
        return None
    original = waters(source)
    available = waters(prepared)
    if (not original or set(original) != set(available) or
            not all(valid_water(atoms) for atoms in original.values())):
        return None
    output = qm.parent / "rescues/validated_source_water_failed_fragment_cohort"
    output.mkdir(parents=True, exist_ok=True)
    candidate = output / "source_water_prepared.pdb"
    if not candidate.exists():
        shutil.copyfile(prepared, candidate)
        audit = preserve_source_waters(source, candidate)
    else:
        audit = json.loads(candidate.with_suffix(".water_hydrogens.json").read_text())
    if (audit.get("source_sha256") != digest(source) or
            audit.get("prepared_pretransfer_sha256") != digest(prepared) or
            audit.get("prepared_sha256") != digest(candidate) or
            audit.get("n_preserved_waters") != len(original) or
            audit.get("n_pqr_fallback_waters") != 0):
        raise ValueError("validated source-water rescue changed after freezing")
    return candidate, audit


def transferred_water_geometry(prepared: Path, poses_root: Path,
                               names: list[str]) -> dict:
    hydrogen_points = np.asarray([xyz(line) for atoms in waters(prepared).values()
                                  for name, line in atoms.items() if name in ("H1", "H2")])
    if not len(hydrogen_points):
        raise ValueError("source-water geometry audit has no hydrogens")
    per_ligand = []
    for ligand_name in names:
        pose = poses_root / "minimized" / ligand_name / "minimized_complex.pdb"
        _, ligand_atoms = minimized_atoms(pose)
        tree = cKDTree(np.asarray([atom["coordinates"] for atom in ligand_atoms]))
        distances, partners = tree.query(hydrogen_points, k=1)
        severe = int(sum(distance < (0.80 if ligand_atoms[index]["element"] == "H" else 1.15)
                         for distance, index in zip(distances, partners, strict=True)))
        per_ligand.append({"ligand": ligand_name, "n_severe_source_water_hydrogen_contacts": severe,
                           "nearest_ligand_distance_angstrom": float(min(distances))})
    return {"schema": "fmo-frozen-poses-source-water-contact-audit-v1",
            "prepared_sha256": digest(prepared),
            "minimized_complex_sha256": {name: digest(poses_root / "minimized" / name /
                                                     "minimized_complex.pdb") for name in names},
            "n_ligands": len(names), "n_water_hydrogens": len(hydrogen_points),
            "all_safe": all(not item["n_severe_source_water_hydrogen_contacts"]
                            for item in per_ligand), "cases": per_ligand}


def rescue_cohort(target: str, qm: Path, names: list[str], library: Path, ncores: int,
                  timeout: int, workers: int, gamess: Path, params: Path,
                  report_filename: str = "adaptive_status.json") -> dict:
    if report_filename not in {"adaptive_status.json", "fragment_repair_status.json"}:
        raise ValueError("unregistered cohort report filename")
    failure = eligible_failure(qm, names)
    if failure is None:
        return {"status": "REVIEW_REQUIRED", "reason": "unregistered strict failure",
                "ranking_eligible": False}
    strict = qm / "full_target_fmo"
    protocol_path = strict / "frozen_protocol.json"
    pocket_path = qm / "pocket_selection_audit.json"
    prepared = strict / f"{target.lower()}_prepared.pdb"
    poses = strict / "ligands_no_properties.sdf"
    if not all(path.is_file() for path in (protocol_path, pocket_path, prepared, poses)):
        raise ValueError("missing frozen strict receptor, poses or pocket audit")
    protocol = json.loads(protocol_path.read_text())
    pocket = json.loads(pocket_path.read_text())
    strict_status = json.loads((qm / "run_status.json").read_text())
    if protocol.get("reference_water_hydrogens_sha256") is not None:
        source = qm / "input_target_water_treated.pdb"
        if (not source.is_file() or
                protocol["reference_water_hydrogens_sha256"] != digest(source)):
            raise ValueError("source-water protocol fingerprint changed")
    if (strict_status.get("repair_peptide_numbering") is not True or
            strict_status.get("pocket_heavy_only") is not False or
            not preparation_stages_verified(qm, strict_status, prepared) or
            strict_status.get("pocket_selection_audit_sha256") != digest(pocket_path) or
            protocol.get("target_name") != target or protocol.get("ligand_names") != names or
            protocol.get("repair_peptide_numbering") is not True or
            protocol.get("pocket_heavy_only") or
            pocket.get("ligand_names") != names or
            pocket.get("selection_production") != "all_atom" or
            pocket.get("source_sha256", {}).get("receptor") != digest(prepared) or
            pocket.get("source_sha256", {}).get("ligands") != digest(poses) or
            pocket.get("source_sha256", {}).get("frozen_protocol") != digest(protocol_path) or
            not pocket.get("membership_differences") and not numerical_soscf_eligible(failure)):
        raise ValueError("strict target pocket or numbering audit does not permit heavy-pocket rescue")
    residue = diagnosed_peripheral_residue(qm, failure["failed_cases"], pocket)
    source_water = None if strict_status.get("preserve_input_water_hydrogens") else (
        source_water_rescue(qm, prepared, residue, failure))
    if source_water is not None:
        water_geometry = transferred_water_geometry(source_water[0], strict, names)
        report_file = source_water[0].parent / "source_water_post_k10_geometry.json"
        if report_file.is_file() and json.loads(report_file.read_text()) != water_geometry:
            raise ValueError("cached source-water contact audit changed")
        report_file.write_text(json.dumps(water_geometry, indent=2) + "\n")
        if not water_geometry["all_safe"]:
            source_water = None
    case = qm.parent
    report_path = case / report_filename
    if report_path.is_file():
        previous = json.loads(report_path.read_text())
        prediction = case / "rescues" / previous.get("selected_arm", "") / "blind_predictions.json"
        if previous.get("status") == "RESCUED_ALTERNATE_PROTOCOL":
            if (previous.get("strict_status_sha256") != failure["status_sha256"] or
                    previous.get("strict_protocol_sha256") != digest(protocol_path) or
                    not prediction.is_file() or
                    previous.get("blind_predictions_sha256") != digest(prediction)):
                raise ValueError("cached cohort rescue inputs or output changed")
            if previous.get("selected_arm") == "validated_source_water_failed_fragment_cohort":
                if source_water is None:
                    raise ValueError("cached source-water rescue is no longer eligible")
                water_report = source_water[0].with_suffix(".water_hydrogens.json")
                geometry = source_water[0].parent / "source_water_post_k10_geometry.json"
                if (not water_report.is_file() or not geometry.is_file() or
                        previous.get("arms", [{}])[-1].get("water_hydrogen_audit_sha256") !=
                        digest(water_report) or
                        previous["arms"][-1].get("water_geometry_audit_sha256") != digest(geometry)):
                    raise ValueError("cached source-water audit changed")
            return previous
        raise ValueError("previous incomplete cohort rescue must be reviewed before retry")
    arms = []
    if numerical_soscf_eligible(failure):
        if (strict / "cohort_failure_gate.json").is_file():
            arms.append(("whole_cohort_fmo_soscf_strict_free_guess", None, prepared))
            arms.append(("whole_cohort_fmo_soscf_free_guess", None, prepared))
        arms.append(("whole_cohort_fmo_soscf_early", None, prepared))
        arms.append(("whole_cohort_fmo_soscf", None, prepared))
    if source_water is not None:
        arms.append(("validated_source_water_failed_fragment_cohort", residue,
                     source_water[0]))
    if residue:
        arms.append(("failed_fragment_cohort", residue, prepared))
    if ligand_scf_failures(qm, failure["failed_cases"]):
        arms.append(("whole_cohort_ligand_scf_diis", None, prepared))
    if pocket.get("membership_differences"):
        arms.append(("heavy_pocket_cohort", None, prepared))
    valid = False
    attempts = []
    for arm, excluded, receptor in arms:
        output = case / "rescues" / arm
        common = [sys.executable, *RUNNER, "--source", str(qm / "source"),
                  "--set-name", "newtarget", "--target-name", target,
                  "--expected", str(len(names)), "--output", str(output),
                  "--prepared-protein", str(receptor), "--ligand-sdf", str(poses),
                  "--poses-root", str(strict), "--workers", str(workers),
                  "--ncores", str(ncores), "--timeout", str(timeout),
                  "--merge-disulfides", "--repair-peptide-numbering",
                  *cohort_arm_controls(arm, excluded)]
        stages = {}
        for stage in ("prepare", "score"):
            command = [*common, "--stage", stage]
            if stage == "score":
                command += ["--gamess-dir", str(gamess), "--library-dir", str(library),
                            "--dftb-parameter-dir", str(params)]
            output.mkdir(parents=True, exist_ok=True)
            with (output / f"{stage}.stdout").open("w") as handle:
                stages[stage] = subprocess.run(command, cwd=ROOT, stdout=handle,
                                               stderr=subprocess.STDOUT, check=False).returncode
            if stages[stage]:
                break
        prediction = output / "blind_predictions.json"
        if stages.get("score") == 0 and prediction.is_file():
            scores = json.loads(prediction.read_text())
            frozen = json.loads((output / "frozen_protocol.json").read_text())
            valid = (scores.get("n_requested") == len(names) == scores.get("n_succeeded") and
                     [row.get("ligand") for row in scores.get("predictions", [])] == names and
                     all(row.get("status") == "SUCCEEDED" for row in scores["predictions"]) and
                     frozen.get("pocket_heavy_only", False) is (arm == "heavy_pocket_cohort") and
                     frozen.get("excluded_pocket_residue") == excluded and
                     frozen.get("scf_solver", "soscf") ==
                     ("fmo_soscf" if arm in NUMERICAL_ARMS else
                      "diis" if arm == "whole_cohort_ligand_scf_diis" else "soscf") and
                     frozen.get("scc_mixer", "broyden") ==
                     ("gamess_soscf" if arm in NUMERICAL_ARMS else "broyden") and
                     frozen.get("soscf_onset", 0.25) == NUMERICAL_ARMS.get(arm, 0.25) and
                     frozen.get("nguess", 2) == (3 if arm in FREE_GUESS_ARMS else 2) and
                     frozen.get("cnvdmp", 50) == (0 if arm == "whole_cohort_fmo_soscf_strict_free_guess" else 50) and
                     frozen.get("repair_peptide_numbering") is True and
                     frozen.get("prepared_protein_override_sha256") == digest(receptor) and
                     frozen.get("poses_protocol_sha256") == digest(protocol_path))
        attempts.append({"arm": arm, "excluded_pocket_residue": excluded, "stages": stages,
                         "scf_solver": "fmo_soscf" if arm in NUMERICAL_ARMS else
                         "diis" if arm == "whole_cohort_ligand_scf_diis" else "soscf",
                         "soscf_onset": NUMERICAL_ARMS.get(arm, 0.25),
                         "nguess": 3 if arm in FREE_GUESS_ARMS else 2,
                         "cnvdmp": 0 if arm == "whole_cohort_fmo_soscf_strict_free_guess" else 50,
                         "complete": valid,
                         **({"water_hydrogen_audit_sha256": digest(
                             receptor.with_suffix(".water_hydrogens.json")),
                             "water_geometry_audit_sha256": digest(
                                 receptor.parent / "source_water_post_k10_geometry.json")}
                            if receptor != prepared else {})})
        if valid:
            break
    report = {"schema": "fmo-target-cohort-rescue-v1", "target": target,
              "n_ligands": len(names), "selected_arm": arm if valid else None,
              "diagnosed_failing_residue": residue,
              "strict_status_sha256": failure["status_sha256"],
              "strict_protocol_sha256": digest(protocol_path),
              "pocket_audit_sha256": digest(pocket_path),
              "failed_cases": failure["failed_cases"], "arms": attempts, "stages": stages,
              "status": "RESCUED_ALTERNATE_PROTOCOL" if valid else "REVIEW_REQUIRED",
              "ranking_eligible": valid}
    if valid:
        report["blind_predictions_sha256"] = digest(prediction)
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    return report
