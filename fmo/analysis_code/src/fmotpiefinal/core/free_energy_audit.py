"""Numerically stable diagnostics, not a certificate of ensemble convergence."""
from __future__ import annotations

import math

import numpy as np
from scipy.special import logsumexp


GAS_CONSTANT_KCAL_MOL_K = 0.00198720425864083


def perturbation_summary(delta_u_kcal_mol, *, temperature_kelvin,
                         statistical_inefficiency=None, energy_kind="total_potential"):
    if energy_kind != "total_potential":
        raise ValueError("a full Hamiltonian correction requires total_potential differences")
    values = np.asarray(delta_u_kcal_mol, dtype=float)
    if values.ndim != 1 or len(values) < 2 or not np.isfinite(values).all():
        raise ValueError("at least two finite, ordered energy differences required")
    if not math.isfinite(temperature_kelvin) or temperature_kelvin <= 0:
        raise ValueError("positive finite temperature required")
    if statistical_inefficiency is not None and (
            not math.isfinite(statistical_inefficiency) or statistical_inefficiency < 1):
        raise ValueError("statistical inefficiency must be finite and at least one")
    thermal_energy = GAS_CONSTANT_KCAL_MOL_K * temperature_kelvin
    log_weights = -values / thermal_energy
    normalizer = logsumexp(log_weights)
    weights = np.exp(log_weights - normalizer)
    effective_samples = float(1.0 / np.sum(weights ** 2))
    return {
        "schema": "one-sided-total-potential-perturbation-diagnostic-v1",
        "n_frames": len(values),
        "temperature_kelvin": temperature_kelvin,
        "delta_free_energy_kcal_mol": float(-thermal_energy * (normalizer - math.log(len(values)))),
        "mean_delta_u_kcal_mol": float(values.mean()),
        "std_delta_u_kcal_mol": float(values.std(ddof=1)),
        "weight_ess": effective_samples,
        "weight_ess_fraction": effective_samples / len(values),
        "largest_normalized_weight": float(weights.max()),
        "statistical_inefficiency": statistical_inefficiency,
        "heuristic_correlation_adjusted_weight_ess": (
            effective_samples / statistical_inefficiency
            if statistical_inefficiency is not None else None),
        "convergence_proven": False,
        "warning": "Weight ESS alone cannot detect unvisited states; require equilibrium provenance, independent repeats and overlap diagnostics.",
    }


def relative_binding_difference(bound_kcal_mol, solution_kcal_mol):
    if not all(math.isfinite(value) for value in (bound_kcal_mol, solution_kcal_mol)):
        raise ValueError("finite estimates for both thermodynamic legs required")
    return bound_kcal_mol - solution_kcal_mol


def relative_qmmm_correction(bound_a_delta_u_kcal_mol, bound_b_delta_u_kcal_mol,
                             solution_a_delta_u_kcal_mol, solution_b_delta_u_kcal_mol, *,
                             temperature_kelvin, energy_kind,
                             statistical_inefficiencies=None):
    """Diagnose the four-leg QM/MM correction to relative binding free energy.

    Each array contains same-frame ``U_high - U_low`` values for one endpoint
    state (bound/solution, ligand A/B), sampled under the corresponding low
    Hamiltonian.  The returned cycle is
    ``(corr_bound_B - corr_bound_A) - (corr_solution_B - corr_solution_A)``.
    One-sided exponential averaging is diagnostic-only: it cannot prove
    overlap, equilibrium, or convergence.
    """
    if energy_kind != "total_potential":
        raise ValueError("a relative QM/MM free-energy correction requires total_potential differences")
    arrays = {
        "bound_A": bound_a_delta_u_kcal_mol,
        "bound_B": bound_b_delta_u_kcal_mol,
        "solution_A": solution_a_delta_u_kcal_mol,
        "solution_B": solution_b_delta_u_kcal_mol,
    }
    if statistical_inefficiencies is None:
        statistical_inefficiencies = {}
    if not isinstance(statistical_inefficiencies, dict):
        raise ValueError("statistical_inefficiencies must map endpoint state names to values")
    estimates = {
        name: perturbation_summary(
            values,
            temperature_kelvin=temperature_kelvin,
            statistical_inefficiency=statistical_inefficiencies.get(name),
            energy_kind=energy_kind,
        )
        for name, values in arrays.items()
    }
    bound = estimates["bound_B"]["delta_free_energy_kcal_mol"] - estimates["bound_A"]["delta_free_energy_kcal_mol"]
    solution = estimates["solution_B"]["delta_free_energy_kcal_mol"] - estimates["solution_A"]["delta_free_energy_kcal_mol"]
    return {
        "schema": "relative-qmmm-total-potential-correction-diagnostic-v1",
        "endpoint_corrections": estimates,
        "bound_relative_correction_kcal_mol": bound,
        "solution_relative_correction_kcal_mol": solution,
        "delta_binding_correction_kcal_mol": bound - solution,
        "convergence_proven": False,
        "ranking_eligible": False,
        "warning": (
            "This is a Hamiltonian correction diagnostic, not a binding free-energy "
            "estimate. Require bidirectional/bridged overlap, independent repeats, "
            "chemical-state coverage and endpoint provenance before ranking."
        ),
    }
