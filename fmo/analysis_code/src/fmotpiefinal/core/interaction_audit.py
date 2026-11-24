"""Label-free diagnostics for FMO interface energies and free-energy gaps.

The module deliberately does not fit coefficients or turn an endpoint energy
into a binding free energy.  Its purpose is to make the bookkeeping explicit:
PIEDA components decompose an interface score, whereas three-state and
solvation quantities are separate observables that require a defined
thermodynamic cycle before they can be combined.
"""
from __future__ import annotations

import csv
import hashlib
import math
from pathlib import Path
from statistics import mean, pstdev


PIEDA_COMPONENTS = ("ES", "EX", "CT", "DI", "SOLV")
OPTIONAL_OBSERVABLES = ("three_state", "pb_polar", "sa_nonpolar", "total_potential")

_ALIASES = {
    "case": ("case", "id", "ligand"),
    "target": ("target",),
    "tpie": ("tpie", "tpie_shared", "tpie_per_ligand", "interaction_kcal_mol"),
    "three_state": ("three_state", "qm_three_state_shared", "three_state_delta"),
    "pb_polar": ("pb_polar", "pb_polar_shared", "pbsa_polar"),
    "sa_nonpolar": ("sa_nonpolar", "sa_nonpolar_shared", "surface_area_nonpolar"),
    "total_potential": ("total_potential", "delta_u", "delta_u_kcal_mol"),
}
for _component in PIEDA_COMPONENTS:
    _ALIASES[_component] = (_component, _component.lower(), f"{_component}_kcal_mol")


def _number(value, field, row_number):
    if value is None or value == "":
        return None
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"row {row_number}: {field} must be numeric") from error
    if not math.isfinite(number):
        raise ValueError(f"row {row_number}: {field} must be finite")
    return number


def _first_value(row, names):
    for name in names:
        if name in row and row[name] not in (None, ""):
            return row[name]
    return None


def normalize_rows(rows):
    """Normalize CSV/JSON records without reading experimental labels."""
    normalized = []
    for row_number, source in enumerate(rows, start=1):
        if not isinstance(source, dict):
            raise ValueError(f"row {row_number}: expected an object")
        target = _first_value(source, _ALIASES["target"])
        case = _first_value(source, _ALIASES["case"])
        if case is None:
            case = f"row-{row_number}"
        if target:
            case = f"{target}/{case}"
        values = {field: _number(_first_value(source, aliases), field, row_number)
                  for field, aliases in _ALIASES.items()
                  if field not in ("case", "target")}
        if values["tpie"] is None and not all(values[name] is not None for name in PIEDA_COMPONENTS):
            raise ValueError(f"row {row_number}: tPIE or all five PIEDA components are required")
        normalized.append({"case": str(case), **values})
    if not normalized:
        raise ValueError("at least one energy row is required")
    return normalized


def _summary(values):
    if not values:
        return {"n": 0, "mean": None, "std": None, "min": None, "max": None, "range": None}
    return {"n": len(values), "mean": mean(values), "std": pstdev(values),
            "min": min(values), "max": max(values), "range": max(values) - min(values)}


def audit_rows(rows, *, closure_tolerance=0.05):
    """Return label-free bookkeeping and cancellation diagnostics.

    ``closure_tolerance`` is only for printed/serialized endpoint values.  It
    is not a physical uncertainty estimate and is never used to improve a
    ranking.
    """
    if not math.isfinite(closure_tolerance) or closure_tolerance <= 0:
        raise ValueError("closure tolerance must be positive and finite")
    normalized = normalize_rows(rows)
    closure_rows = []
    failures = []
    cancellation = []
    for row in normalized:
        components = [row[name] for name in PIEDA_COMPONENTS]
        if row["tpie"] is not None and all(value is not None for value in components):
            component_sum = sum(components)
            residual = row["tpie"] - component_sum
            closure_rows.append({"case": row["case"], "tpie": row["tpie"],
                                 "component_sum": component_sum, "residual": residual,
                                 "pass": abs(residual) <= closure_tolerance})
            if abs(residual) > closure_tolerance:
                failures.append(closure_rows[-1])
            absolute_sum = sum(abs(value) for value in components)
            cancellation.append({"case": row["case"],
                                 "ratio": abs(component_sum) / absolute_sum
                                 if absolute_sum else None,
                                 "absolute_component_sum": absolute_sum,
                                 "dominant_component": PIEDA_COMPONENTS[max(
                                     range(len(components)), key=lambda index: abs(components[index]))]})

    component_summaries = {
        name: _summary([row[name] for row in normalized if row[name] is not None])
        for name in ("tpie",) + PIEDA_COMPONENTS + OPTIONAL_OBSERVABLES
    }
    ratios = [item["ratio"] for item in cancellation if item["ratio"] is not None]
    dominant_counts = {name: sum(item["dominant_component"] == name for item in cancellation)
                       for name in PIEDA_COMPONENTS}
    optional_present = {name: sum(row[name] is not None for row in normalized)
                        for name in OPTIONAL_OBSERVABLES}
    warnings = [
        "PIEDA components are a decomposition of tPIE; adding tPIE to its components double-counts the same interface energy.",
        "Three-state, PB polar, and nonpolar surface terms are separate observables; they are not additive corrections without a declared thermodynamic cycle.",
        "A single endpoint interaction score is not a total Hamiltonian difference and cannot be passed to a perturbation estimator as delta-U.",
        "Experimental columns are not accessed, used, or fitted by this audit.",
    ]
    if optional_present["total_potential"]:
        warnings.append("total_potential values are inventoried only; ensemble provenance, overlap, and independent-repeat checks are still required.")
    status = "PASS_BOOKKEEPING"
    if failures:
        status = "REVIEW_COMPONENT_CLOSURE"
    elif len(closure_rows) != len(normalized):
        status = "REVIEW_INCOMPLETE_COMPONENT_COVERAGE"
    return {
        "schema": "fmotpiefinal-interaction-energy-audit-v1",
        "status": status,
        "ranking_eligible": False,
        "n_rows": len(normalized),
        "component_closure": {
            "n_complete": len(closure_rows),
            "n_unchecked": len(normalized) - len(closure_rows),
            "n_pass": len(closure_rows) - len(failures),
            "n_fail": len(failures),
            "tolerance_kcal_mol": closure_tolerance,
            "max_abs_residual_kcal_mol": max((abs(item["residual"]) for item in closure_rows), default=None),
            "failures": failures[:20],
        },
        "label_free_diagnostics": {
            "summaries_kcal_mol": component_summaries,
            "net_to_gross_component_ratio": _summary(ratios),
            "dominant_component_counts": dominant_counts,
            "optional_observable_coverage": optional_present,
        },
        "fep_gap": {
            "missing_from_endpoint_table": [
                "bound_state_ensemble",
                "solution_state_ensemble",
                "matched_total_potential_delta_u",
                "alchemical_or_free_energy_estimator",
                "overlap_and_effective_sample_size",
                "independent_repeat_convergence",
                "chemical_state_population_model",
            ],
            "ranking_eligible": False,
        },
        "overlap_warnings": warnings,
    }


def audit_csv(path: Path, *, closure_tolerance=0.05):
    with path.open(newline="", encoding="utf-8") as handle:
        report = audit_rows(csv.DictReader(handle), closure_tolerance=closure_tolerance)
    report["source"] = {"path": str(path),
                        "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    return report
