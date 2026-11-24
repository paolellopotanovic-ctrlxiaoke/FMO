"""Explicit optional dependency gates for the lightweight test environment."""
from importlib.util import find_spec

import pytest


COMPUTE_MODULES = {"test_special_chemistry.py", "test_component_inventory.py",
                   "test_counterion_preparation.py", "test_chemical_preflight.py"}
COMPUTE_CASES = {"test_component_hash_change_blocks_reuse",
                 "test_water_only_mm_fallback_fixes_clash_without_moving_oxygen"}


def pytest_collection_modifyitems(items):
    if find_spec("openmmforcefields") is None or find_spec("openff") is None:
        marker = pytest.mark.skip(reason="full compute integration test requires OpenFF/OpenMMForceFields")
        for item in items:
            if item.path.name in COMPUTE_MODULES or item.name in COMPUTE_CASES:
                item.add_marker(marker)
