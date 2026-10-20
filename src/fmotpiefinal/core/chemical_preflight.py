"""Persist label-free special-component plans before filtering or expensive QM."""
import json
from pathlib import Path
import time

from .component_inventory import inventory, sha
from .special_dispatch import choose_route


def component_plan(graph):
    special = [row for row in graph["residues"] if row["kind"] not in {"standard_residue", "water"}]
    try:
        route = choose_route(graph)
        error = None
    except ValueError as exception:
        route, error = "UNSUPPORTED_COMPONENT_COMBINATION", str(exception)
    return {"schema": "special-component-plan-v1", "source_sha256": graph["source_sha256"],
        "experimental_labels_used": False, "coordinates_modified": False,
        "components": [{key: row[key] for key in ("id", "model", "kind", "elements", "n_heavy_atoms")} for row in special],
        "available_mm_route": route, "blocking_reason": error,
        "status": "STANDARD_PATH" if not special and not error else "UNSUPPORTED_CHEMISTRY" if error else "VALIDATED_MM_PREPARATION_REQUIRED",
        "automatic_qm_approval": False,
        "required_followup": [] if not special and not error else [
            "Preserve original components and verify source-to-prepared mapping.",
            "Validate complete MM model before k10; do not delete unsupported atoms.",
            "Validate pocket membership, fragment charges and QM representation separately; MM support is not QM approval."],
        "ranking_eligible": False}


def write_component_plan(graph, output):
    output = Path(output)
    report = component_plan(graph)
    text = json.dumps(report, indent=2) + "\n"
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() and output.read_text() != text:
        raise ValueError("component plan changed; use a new run directory")
    output.write_text(text)
    return report


def audit_chemistry(raw, output, prepared=None, parameter_dir=None, ccd_cache=None):
    from .chemical_graph_policy import plan_actions
    from .receptor_chemical_graph import audit, raw_component_retention
    raw, output = Path(raw), Path(output)
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    graph = inventory(raw)
    (output / "raw_component_inventory.json").write_text(json.dumps(graph, indent=2) + "\n")
    plan = write_component_plan(graph, output / "component_model_plan.json")
    from .chemical_model_requirements import model_requirements
    requirements = model_requirements(graph, parameter_dir)
    (output / "chemical_model_requirements.json").write_text(json.dumps(requirements, indent=2) + "\n")
    receptor = Path(prepared) if prepared else raw
    report = {"schema": "toolkit-chemical-preflight-v1", "experimental_labels_used": False,
        "raw_sha256": graph["source_sha256"], "audited_sha256": sha(receptor),
        "component_plan": plan, "coordinates_modified": False, "ranking_eligible": False,
        "scope": "diagnostic only; no automatic chemical repair or QM model approval"}
    try:
        chemistry = audit(receptor, parameter_dir, ccd_cache)
        retention = raw_component_retention(raw, chemistry) if prepared else None
        actions = plan_actions(chemistry, retention)
        for name, value in (("receptor_chemical_graph", chemistry), ("chemical_action_plan", actions),
                            ("raw_component_retention", retention)):
            (output / f"{name}.json").write_text(json.dumps(value, indent=2) + "\n")
        report.update({"status": actions["quality_status"], "summary": actions["summary"]})
    except (ValueError, KeyError) as error:
        report.update({"status": "UNRESOLVED_CHEMICAL_GRAPH", "error": str(error)})
    if plan["status"] == "UNSUPPORTED_CHEMISTRY":
        report["status"] = "BLOCKED_UNSUPPORTED_CHEMISTRY"
    report["elapsed_seconds"] = time.monotonic() - started
    (output / "chemical_preflight.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
