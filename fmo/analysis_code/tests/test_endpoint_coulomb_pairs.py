from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / "tools"))
from analyze_endpoint_coulomb_pairs import pair_parameters


def test_pair_parameters_preserve_exclusions_and_scaled_exceptions():
    import openmm

    force = openmm.NonbondedForce()
    for charge in (.5, -.5, .2):
        force.addParticle(charge, .3, .1)
    force.addException(0, 1, 0, .3, 0)
    force.addException(2, 0, .03, .3, .1)
    first, second, products = pair_parameters(force)
    assert list(zip(first, second)) == [(0, 1), (0, 2), (1, 2)]
    assert products.tolist() == pytest.approx([0, .03, -.1])
