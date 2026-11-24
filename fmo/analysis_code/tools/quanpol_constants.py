"""Shared defaults and safety checks for GAMESS QuanPol inputs."""
from __future__ import annotations

import math

DEFAULT_SWRB_ANGSTROM = 6000.0
QM_CENTER_RESTRAINT_THRESHOLD_ANGSTROM = 10000.0 * 0.529177210903


def validate_swrb(swrb, allow_qm_center_restraint=False):
    if not math.isfinite(swrb) or swrb <= 0:
        raise ValueError("SWRB must be a finite positive value")
    if swrb <= QM_CENTER_RESTRAINT_THRESHOLD_ANGSTROM and not allow_qm_center_restraint:
        raise ValueError("SWRB activates GAMESS' QM-center restraint branch")
