"""Raw-receptor chemistry contract evaluated before chain or HETATM filtering."""
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path


STANDARD = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL ASH GLH HID HIE HIP CYX CYM LYN".split())
WATERS = {"HOH", "WAT", "TIP", "H2O"}
METALS = set("LI NA K RB CS MG CA SR BA ZN FE MN CU CO NI CD HG AL CR V MO W".split())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def inventory(source):
    groups = defaultdict(list)
    models = set()
    model = "1"
    for line in Path(source).read_text().splitlines():
        if line.startswith("MODEL "):
            model = line[10:14].strip()
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        models.add(model)
        chain, name, number = line[21:22].strip() or "_", line[17:20].strip(), line[22:27].strip()
        element = line[76:78].strip().upper()
        if not element:
            raise ValueError("raw receptor inventory requires explicit element columns")
        coordinates = [float(line[start:start + 8]) for start in (30, 38, 46)]
        if not all(math.isfinite(value) for value in coordinates):
            raise ValueError("non-finite raw receptor coordinate")
        groups[(model, f"{chain}:{name}:{number}")].append({"name": line[12:16].strip(), "element": element,
            "record": line[:6].strip(), "altloc": line[16:17].strip(), "xyz": coordinates})
    if not groups:
        raise ValueError("raw receptor has no atoms")
    records = []
    for (model, identifier), atoms in groups.items():
        name = identifier.split(":")[1]
        kind = ("standard_residue" if name in STANDARD else "water" if name in WATERS
                else "metal" if any(atom["element"] in METALS for atom in atoms)
                else "ion" if len(atoms) == 1 and atoms[0]["element"] in {"F", "CL", "BR", "I"}
                else "nonstandard_component")
        records.append({"id": identifier, "model": model, "kind": kind, "elements": dict(Counter(atom["element"] for atom in atoms)),
            "n_heavy_atoms": sum(atom["element"] not in {"H", "D"} for atom in atoms), "atoms": atoms})
    return {"schema": "fmo-raw-component-inventory-v1", "source_sha256": sha(source),
        "experimental_labels_used": False, "n_models": len(models), "residues": records}


def enforce_inventory(source, output, decisions=None):
    graph = inventory(source)
    from .chemical_preflight import write_component_plan
    write_component_plan(graph, Path(output).with_name("component_model_plan.json"))
    policy = json.loads(Path(decisions).read_text()) if decisions else {"components": []}
    if decisions and (policy.get("source_sha256") != graph["source_sha256"] or policy.get("experimental_labels_used") is not False):
        raise ValueError("component decisions must match raw source and be label-free")
    choices = policy.get("components", [])
    keys = [(row.get("model", "1"), row["id"]) for row in choices]
    if len(set(keys)) != len(keys):
        raise ValueError("duplicate component decisions")
    by_key = dict(zip(keys, choices, strict=True))
    present = {(row["model"], row["id"]) for row in graph["residues"]}
    if not set(keys) <= present:
        raise ValueError("component decision does not match an actual raw residue")
    tasks = []
    from .special_dispatch import choose_route
    try:
        registered_route = choose_route(graph)
        if registered_route == "PHOSPHORESIDUE_ALKALI_ALKALINE_EARTH_MM":
            from .special_chemistry import PHOSPHO_RESIDUES, validate_phosphoresidue_identity
            for row in graph["residues"]:
                if row["id"].split(":")[1] in PHOSPHO_RESIDUES:
                    validate_phosphoresidue_identity(row)
    except ValueError:
        registered_route = None
    for row in graph["residues"]:
        if row["kind"] in {"standard_residue", "water"}:
            continue
        decision = by_key.get((row["model"], row["id"]))
        if registered_route is not None and decision is None:
            continue
        if decision is None:
            tasks.append({"component": row["id"], "model": row["model"], "kind": row["kind"],
                "action": "RETAIN_AND_RESOLVE_MODEL", "reason": "special component cannot silently enter a filtering-only protein path"})
            continue
        if not isinstance(decision.get("reason"), str) or len(decision["reason"].strip()) < 10:
            raise ValueError("component decision needs a substantive structural justification")
        if (decision.get("action") == "EXCLUDE" and row["kind"] in {"ion", "nonstandard_component"}
                and all(atom["record"] == "HETATM" for atom in row["atoms"])):
            if not decision.get("evidence_sha256") or len(decision["evidence_sha256"]) != 64:
                raise ValueError("exclusion requires a pinned structural evidence artifact")
            evidence = Path(decision["evidence_path"])
            if sha(evidence) != decision["evidence_sha256"]:
                raise ValueError("component exclusion evidence changed")
        else:
            tasks.append({"component": row["id"], "model": row["model"], "kind": row["kind"],
                "action": "VALIDATED_SPECIAL_CHEMISTRY_ROUTE_REQUIRED", "requested_decision": decision,
                "reason": "generic protein path cannot validate a metal model or nonstandard conversion"})
    if graph["n_models"] != 1:
        tasks.append({"action": "REGISTER_SINGLE_MODEL_BEFORE_PREPARATION"})
    report = {**graph, "status": "BLOCKED_COMPONENT_MODEL_REQUIRED" if tasks else "RAW_COMPONENT_CONTRACT_PASS",
        "tasks": tasks, "decisions": choices, "decision_sha256": sha(decisions) if decisions else None,
        "coordinates_modified": False, "scope": "pre_filter_generic_protein_path",
        "registered_mm_route": registered_route}
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(report, indent=2) + "\n"
    if output.exists() and output.read_text() != text:
        raise ValueError("raw component audit changed; use a new run directory")
    output.write_text(text)
    if tasks:
        raise ValueError(f"raw receptor contains unresolved chemistry; retained inventory and required actions: {output}")
    return report
