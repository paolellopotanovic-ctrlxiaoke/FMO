"""Audit actual PCM energy composition before constructing a three-state score."""
from __future__ import annotations

import math
import re


HARTREE_TO_KCAL_MOL = 627.5094740631


def pcm_energy_inventory(text: str) -> dict:
    if ("EXECUTION OF GAMESS TERMINATED NORMALLY" not in text or
            "TERMINATED -ABNORMALLY-" in text or "DIMER CALCULATION DID NOT CONVERGE" in text):
        raise ValueError("normal completed PCM calculation required")
    patterns = {
        "electrostatic": r"ELECTROSTATIC INTERACTION\s*=\s*([-+\d.]+)\s*KCAL/MOL",
        "cavitation": r"PIEROTTI CAVITATION ENERGY\s*=\s*([-+\d.]+)\s*KCAL/MOL",
        "dispersion": r"DISPERSION FREE ENERGY\s*=\s*([-+\d.]+)\s*KCAL/MOL",
        "repulsion": r"REPULSION FREE ENERGY\s*=\s*([-+\d.]+)\s*KCAL/MOL",
        "total": r"TOTAL INTERACTION\s*=\s*([-+\d.]+)\s*KCAL/MOL",
        "free_hartree": r"Free unco\+D energy in solvent=\s*([-+\d.]+)",
        "internal_hartree": r"Internal unco\+D energy in solvent=\s*([-+\d.]+)",
    }
    values = {}
    for name, pattern in patterns.items():
        matches = re.findall(pattern, text)
        if not matches:
            raise ValueError(f"missing final PCM term: {name}")
        values[name] = float(matches[-1])
    if not all(math.isfinite(value) for value in values.values()):
        raise ValueError("nonfinite PCM energy")
    components = sum(values[name] for name in ("electrostatic", "cavitation", "dispersion", "repulsion"))
    if not math.isclose(components, values["total"], abs_tol=0.005):
        raise ValueError("printed PCM components do not close")
    energy_increment = (values["free_hartree"] - values["internal_hartree"]) * HARTREE_TO_KCAL_MOL
    if not math.isclose(energy_increment, values["total"], abs_tol=0.005):
        raise ValueError("PCM total does not match free-minus-internal energy")
    values["automatic_idisp_disable_messages"] = len(re.findall(r"Setting IDISP=0 for NCAV2=", text))
    return values


def require_consistent_pcm_cycle(decks: dict[str, str], logs: dict[str, str]) -> dict:
    states = {"complex", "receptor", "ligand"}
    if set(decks) != states or set(logs) != states:
        raise ValueError("exact complex/receptor/ligand state triple required")
    options = {}
    models = {}
    inventories = {}
    for state in sorted(states):
        if re.search(r"\bMODMOL\s*=", decks[state], re.I):
            raise ValueError("independent full-energy states required, not selective/subsystem decks")
        matches = re.findall(r"\$PCM\s+(.*?)\$END", decks[state], re.S | re.I)
        if len(matches) != 1:
            raise ValueError("one explicit PCM group required in every state")
        options[state] = dict(re.findall(r"([A-Z]+)\s*=\s*([^\s,]+)", matches[0].upper()))
        models[state] = {}
        for name in ("DFTB", "BASIS", "DFT", "PCMCAV", "TESCAV"):
            model_groups = re.findall(r"\$" + name + r"\s+(.*?)\$END", decks[state], re.S | re.I)
            if len(model_groups) > 1:
                raise ValueError("ambiguous electronic/solvent model group")
            models[state][name] = " ".join(model_groups[0].upper().split()) if model_groups else None
        if not all(name in options[state] for name in ("SOLVNT", "IEF", "ICAV", "IDISP", "IFMO")):
            raise ValueError("explicit solvent controls required")
        inventories[state] = pcm_energy_inventory(logs[state])
    if any(option != options["complex"] for option in options.values()):
        raise ValueError("different declared PCM options across states")
    if any(model != models["complex"] for model in models.values()):
        raise ValueError("different electronic, dispersion, cavity or tessellation model across states")
    if options["complex"]["IDISP"] != "0" and any(
            inventory["automatic_idisp_disable_messages"] for inventory in inventories.values()):
        raise ValueError("requested dispersion was automatically disabled; rerun all states with a consistent supported model")
    if options["complex"]["IDISP"] == "0" and any(
            abs(inventory[name]) > 0.005 for inventory in inventories.values()
            for name in ("dispersion", "repulsion")):
        raise ValueError("unexpected nonzero dispersion/repulsion in IDISP0 cycle")
    return {"states": inventories, "controls": options["complex"], "consistent": True}
