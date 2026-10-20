"""Audit and aggregate FMO interfacial PIEDA for any two fragment sets."""
from __future__ import annotations

import math
from collections.abc import Iterable

from .parser import PairInteraction


COMPONENTS = {
    "ES": "electrostatic_kcal_mol",
    "EX": "exchange_repulsion_kcal_mol",
    "CT": "charge_transfer_kcal_mol",
    "DI": "dispersion_kcal_mol",
    "SOLV": "solvation_kcal_mol",
}


def cross_interface(pairs: Iterable[PairInteraction], group_a: Iterable[int],
                    group_b: Iterable[int], *, max_component_error: float = 0.02) -> dict:
    """Reject missing/duplicate cross pairs and incomplete PIEDA components."""
    first, second = set(group_a), set(group_b)
    if (not first or not second or first & second or
            any(not isinstance(fragment, int) or isinstance(fragment, bool) or fragment < 1
                for fragment in first | second)):
        raise ValueError("two nonempty disjoint positive integer fragment groups are required")
    if max_component_error <= 0 or not math.isfinite(max_component_error):
        raise ValueError("component rounding tolerance must be positive and finite")
    selected = {}
    for pair in pairs:
        fragment_i, fragment_j = pair.fragment_i, pair.fragment_j
        if not ((fragment_i in first and fragment_j in second) or
                (fragment_i in second and fragment_j in first)):
            continue
        key = (fragment_i, fragment_j) if fragment_i in first else (fragment_j, fragment_i)
        if key in selected:
            raise ValueError(f"duplicate interface pair {key}")
        values = [pair.value_kcal_mol] + [getattr(pair, field) for field in COMPONENTS.values()]
        if any(value is None or not math.isfinite(value) for value in values):
            raise ValueError(f"incomplete or non-finite PIEDA row {key}")
        if abs(pair.value_kcal_mol - sum(values[1:])) > max_component_error:
            raise ValueError(f"PIEDA components do not close for pair {key}")
        selected[key] = pair
    expected = {(left, right) for left in first for right in second}
    if set(selected) != expected:
        missing = sorted(expected - set(selected))
        raise ValueError(f"missing {len(missing)} interface pair(s): {missing[:5]}")
    total = float(sum(pair.value_kcal_mol for pair in selected.values()))
    components = {name: float(sum(getattr(pair, field) for pair in selected.values()))
                  for name, field in COMPONENTS.items()}
    return {"n_fragments_a": len(first), "n_fragments_b": len(second),
            "n_cross_pairs": len(expected), "interaction_kcal_mol": total,
            "favorable_score": -total, "components_kcal_mol": components,
            "component_roundoff_kcal_mol": total - sum(components.values()),
            "pair_keys": [list(key) for key in sorted(selected)]}
