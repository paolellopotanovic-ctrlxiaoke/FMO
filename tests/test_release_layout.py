"""The standalone distribution must not depend on the research workspace."""
import ast
from pathlib import Path

from fmotpiefinal import __version__
from fmotpiefinal.core import bound_geometry


def test_release_version_and_live_geometry_helpers():
    assert __version__ == "0.1.0rc1"
    assert callable(bound_geometry.shared_pocket)
    assert callable(bound_geometry.load_bound_ligand)
    assert not hasattr(bound_geometry, "calculate_relative")


def test_no_research_imports_in_runtime():
    root = Path(bound_geometry.__file__).resolve().parents[1]
    for path in root.rglob("*.py"):
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.ImportFrom) and not node.level:
                assert (node.module or "").split(".")[0] not in {
                    "experiments", "validation", "sophosqm_baseline"}
            elif isinstance(node, ast.Import):
                assert all(alias.name.split(".")[0] not in {
                    "experiments", "validation", "sophosqm_baseline"} for alias in node.names)
