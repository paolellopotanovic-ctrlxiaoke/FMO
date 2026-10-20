"""Fail-closed, label-free target-local tPIE soft-cap for completed QM cohorts."""
from __future__ import annotations

import csv
import hashlib
import json
import math
from importlib.resources import files
from pathlib import Path
from statistics import median


FACTOR = 3.0
MAD_CONVERSION = 1.4826
REFERENCE_SCHEMA = "label-free-target-local-tpie-softcap-v1"
BUILTIN = "references/retrospective_targets_20261016.json"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def selected_predictions(root: Path) -> tuple[dict, Path, list[dict]]:
    root = root.resolve()
    status_path = root / "adaptive_status.json"
    status = read_json(status_path)
    if status.get("ranking_eligible") is not True or status.get("status") not in (
            "BLIND_FMO_COMPLETE_REQUIRES_EXTERNAL_VALIDATION", "RESCUED_ALTERNATE_PROTOCOL"):
        raise ValueError("full-cohort adaptive QM coverage has not passed")
    if status["status"] == "RESCUED_ALTERNATE_PROTOCOL":
        selected = status.get("selected_arm")
        arms = [arm for arm in status.get("arms", []) if arm.get("arm") == selected]
        if (not isinstance(selected, str) or not selected or len(arms) != 1 or
                not (arms[0].get("complete") is True or "tpie" in arms[0])):
            raise ValueError("no complete registered rescue arm")
        source = root / "rescues" / selected / "blind_predictions.json"
        if not source.resolve().is_relative_to((root / "rescues").resolve()):
            raise ValueError("rescue arm escapes registered output")
        expected_sha256 = status.get("blind_predictions_sha256", arms[0].get("prediction_sha256"))
    else:
        source = root / "qm/full_target_fmo/blind_predictions.json"
        expected_sha256 = status.get("blind_predictions_sha256")
    if not source.is_file() or digest(source) != expected_sha256:
        raise ValueError("selected QM prediction fingerprint changed or is missing")
    predictions = read_json(source)
    rows = predictions.get("predictions", [])
    gate = read_json(root / "geometry/geometry_gate.json")
    if (gate.get("schema") != "fmo-new-target-geometry-gate-v1" or
            not isinstance(rows, list) or not rows or predictions.get("n_requested") != len(rows) or
            predictions.get("n_succeeded") != len(rows) or gate.get("n_ligands") != len(rows) or
            status.get("n_ligands", len(rows)) != len(rows) or
            status.get("target", gate.get("target")) != gate.get("target")):
        raise ValueError("incomplete QM prediction coverage")
    names = [row.get("ligand") for row in rows]
    if len(set(names)) != len(rows) or not all(isinstance(name, str) and name for name in names):
        raise ValueError("missing or duplicate ligand identities")
    for row in rows:
        if row.get("status") != "SUCCEEDED" or not isinstance(row.get("tpie"), (int, float)) or (
                not math.isfinite(row["tpie"])):
            raise ValueError("invalid FMO tPIE prediction")
    return status, source, rows


def parameters(scores: list[float]) -> dict | None:
    center = median(scores)
    spread = MAD_CONVERSION * median([abs(value - center) for value in scores])
    if not math.isfinite(spread) or spread <= 0:
        return None
    return {"n": len(scores), "median": center, "scaled_mad": spread}


def choose_reference(target: str, scores: list[float], path: Path | None) -> tuple[dict | None, dict]:
    if path is None:
        resource = files("fmotpiefinal").joinpath(BUILTIN)
        builtin = json.loads(resource.read_text(encoding="utf-8"))
        if target not in builtin["reference"]:
            chosen = parameters(scores)
            return chosen, {"mode": "cohort_local_descriptive_only", "factor": FACTOR,
                            "mad_conversion": MAD_CONVERSION, "target": target,
                            "reference": {target: chosen}, "prospective_reference": False}
        path = Path(str(resource))
        mode = "bundled_retrospective_target_reference"
    else:
        mode = "previously_frozen_target_reference"
    frozen = read_json(path)
    if (frozen.get("schema") != REFERENCE_SCHEMA or frozen.get("factor") != FACTOR or
            frozen.get("mad_conversion") != MAD_CONVERSION or target not in frozen.get("reference", {})):
        raise ValueError(f"no compatible frozen reference for target {target}")
    chosen = frozen["reference"][target]
    if (not isinstance(chosen, dict) or
            any(not isinstance(chosen.get(key), (int, float)) or
                not math.isfinite(chosen[key]) for key in ("median", "scaled_mad")) or
            chosen["scaled_mad"] <= 0):
        raise ValueError("invalid target median/MAD in frozen reference")
    return chosen, {"mode": mode, "target": target,
                    "reference_sha256": digest(path),
                    "prospective_reference": False,
                    "reference_development_set_retrospective": (
                        True if mode == "bundled_retrospective_target_reference" else None)}


def score_run(root: Path, reference: Path | None = None) -> dict:
    root = root.resolve()
    status, source, rows = selected_predictions(root)
    target = read_json(root / "geometry/geometry_gate.json").get("target")
    if not isinstance(target, str) or not target:
        raise ValueError("missing target identity in adaptive status")
    scores = [-float(row["tpie"]) for row in rows]
    settings, provenance = choose_reference(target, scores, reference)
    scored = []
    for row, value in zip(rows, scores, strict=True):
        transformed = None
        if settings is not None:
            width = FACTOR * settings["scaled_mad"]
            transformed = settings["median"] + width * math.tanh((value - settings["median"]) / width)
        scored.append({"ligand": row["ligand"], "tpie_kcal_mol": row["tpie"],
                       "raw_score": value, "softcap_score": transformed})
    ranked = sorted(scored, key=lambda entry: (-entry["raw_score"], entry["ligand"]))
    if settings is not None and any(
            first["raw_score"] > second["raw_score"] and
            first["softcap_score"] <= second["softcap_score"]
            for first, second in zip(ranked, ranked[1:])):
        raise ValueError("soft-cap lost ranking precision")
    for rank, row in enumerate(ranked, start=1):
        row["rank"] = rank
    report = {"schema": "fmotpiefinal-validated-tpie-ranking-v1", "target": target,
              "n": len(ranked), "status_sha256": digest(root / "adaptive_status.json"),
              "predictions_sha256": digest(source), "selected_predictions": str(source),
              "scoring": provenance, "higher_is_better": True,
              "softcap_defined": settings is not None, "ranked": ranked}
    from fmotpiefinal.diagnostics import audit
    diagnostic = audit(root, report, rows)
    destination = root / "ranking"
    destination.mkdir(parents=True, exist_ok=True)
    result = destination / "ranking.json"
    text = json.dumps(report, indent=2, ensure_ascii=False) + "\n"
    if result.exists() and result.read_text(encoding="utf-8") != text:
        raise ValueError("frozen ranking differs: use a fresh output directory")
    reference_file = destination / "local_reference.json"
    local_text = None
    if provenance["mode"] == "cohort_local_descriptive_only" and settings is not None:
        local = {"schema": REFERENCE_SCHEMA, "factor": FACTOR,
                 "mad_conversion": MAD_CONVERSION, "reference": {target: settings},
                 "score_inputs_sha256": {target: digest(source)},
                 "score_unit": "dimensionless ranking proxy, not binding free energy"}
        local_text = json.dumps(local, indent=2) + "\n"
        if reference_file.exists() and reference_file.read_text(encoding="utf-8") != local_text:
            raise ValueError("frozen local target reference differs: use a fresh output directory")
    result.write_text(text, encoding="utf-8")
    if local_text is not None:
        reference_file.write_text(local_text, encoding="utf-8")
    with (destination / "ranking.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("rank", "ligand", "tpie_kcal_mol",
                                                    "raw_score", "softcap_score"))
        writer.writeheader()
        writer.writerows(ranked)
    diagnostic_path = destination / "diagnostics.json"
    diagnostic_text = json.dumps(diagnostic, indent=2, ensure_ascii=False) + "\n"
    if diagnostic_path.exists() and diagnostic_path.read_text(encoding="utf-8") != diagnostic_text:
        raise ValueError("ranking diagnostics changed: use a fresh output directory")
    diagnostic_path.write_text(diagnostic_text, encoding="utf-8")
    from fmotpiefinal.diagnostics import pocket_sensitivity
    pocket_report = pocket_sensitivity(root, report)
    pocket_path = destination / "pocket_sensitivity.json"
    pocket_text = json.dumps(pocket_report, indent=2, ensure_ascii=False) + "\n"
    if pocket_path.exists() and pocket_path.read_text(encoding="utf-8") != pocket_text:
        raise ValueError("frozen pocket sensitivity audit differs: use a fresh output directory")
    pocket_path.write_text(pocket_text, encoding="utf-8")
    from fmotpiefinal.preparation_sensitivity import audit_selected_preparation
    preparation_report = audit_selected_preparation(root, report)
    preparation_path = destination / "preparation_sensitivity.json"
    preparation_text = json.dumps(preparation_report, indent=2, ensure_ascii=False) + "\n"
    if preparation_path.exists() and preparation_path.read_text(encoding="utf-8") != preparation_text:
        raise ValueError("frozen preparation sensitivity audit differs: use a fresh output directory")
    preparation_path.write_text(preparation_text, encoding="utf-8")
    return report
