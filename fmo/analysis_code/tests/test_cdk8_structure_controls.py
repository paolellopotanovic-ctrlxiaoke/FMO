import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest


def test_local_sidechains_exclude_backbone_waters_remote_and_ligand():
    pytest.importorskip("openmm")
    pytest.importorskip("openff")
    from fmotpiefinal.core.minimization import _local_sidechain_indices
    residue = SimpleNamespace(chain=SimpleNamespace(id="A"), id="173", name="ASP")
    water = SimpleNamespace(chain=SimpleNamespace(id="A"), id="400", name="HOH")
    remote = SimpleNamespace(chain=SimpleNamespace(id="A"), id="174", name="ASP")
    atoms = [SimpleNamespace(name=name, residue=residue) for name in ("N", "CA", "C", "O", "H", "HA", "CB", "HB2", "OD1")]
    atoms += [SimpleNamespace(name="H1", residue=water), SimpleNamespace(name="CB", residue=remote)]
    atoms += [SimpleNamespace(name="CB", residue=residue)]
    assert _local_sidechain_indices(atoms, 11, (("A", "173", "ASP"),)) == [6, 7, 8]


def test_frozen_arm_names_and_topology_signature():
    tools = Path(__file__).parents[1] / "tools"
    sys.path.insert(0, str(tools))
    spec = importlib.util.spec_from_file_location("structure_controls", tools / "run_cdk8_structure_controls.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.ARMS == ("baseline_reproduction", "local_sidechains", "water_hydrogens")
    deck = "$FMO NFRAG=1 ICHARG(1)=0 FRGNAM(1)=LIG INDAT(1)=0 1 -3 0 $END"
    assert module.topology_signature(deck)["membership"] == ["0", "1", "-3", "0"]


def test_complete_structure_control_evidence():
    import json
    evidence = Path(__file__).parents[1] / "docs/CDK8_STRUCTURE_CONTROL_RESULTS.zh-CN.json"
    if not evidence.exists():
        pytest.skip("structure control evidence not generated")
    report = json.loads(evidence.read_text())
    assert report["baseline_matched"]
    assert not report["frozen_protocol"]["labels_used_for_generation"]
    assert not report["frozen_protocol"]["production_default_changed"]
    for arm, values in report["arms"].items():
        assert values["n_succeeded"] == values["n_expected"] == 10
        assert len(report["coordinate_audit"][arm]) == 10
        assert len(report["pairwise"][arm]["correct"]) + len(report["pairwise"][arm]["wrong"]) + len(report["pairwise"][arm]["tied"]) == 45
        assert all(row["bonds"]["n_flagged"] == 0 for row in report["ligand_geometry"][arm].values())
    assert report["arms"]["local_sidechains"]["metrics"]["spearman"] == pytest.approx(0.6606060606)
