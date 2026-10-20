from copy import deepcopy
import json

import pytest
from rdkit import Chem

from fmotpiefinal.core.tartrate import source_tartrate


def tartrate_record():
    rows = [
        ("O1", "O", (33.408, 8.321, -8.242)), ("O11", "O", (34.880, 8.942, -6.680)),
        ("C1", "C", (33.668, 8.691, -7.108)), ("C2", "C", (32.588, 8.809, -6.036)),
        ("O2", "O", (31.305, 8.554, -6.580)), ("C3", "C", (32.622, 10.179, -5.349)),
        ("O3", "O", (32.448, 11.219, -6.298)), ("C4", "C", (31.569, 10.256, -4.245)),
        ("O4", "O", (30.725, 11.146, -4.218)), ("O41", "O", (31.676, 9.314, -3.334)),
        ("H2", "H", (32.819, 8.136, -5.211)), ("HO", "H", (30.643, 8.634, -5.889)),
        ("H3", "H", (33.556, 10.285, -4.796)), ("HO", "H", (31.563, 11.171, -6.667))]
    return {"id": "A:TLA:502", "model": "1", "kind": "nonstandard_component",
            "atoms": [{"name": name, "element": element, "xyz": list(point), "altloc": ""} for name, element, point in rows]}


def test_source_tartrate_preserves_state_coordinates_and_resolves_duplicate_hydrogen_names():
    record = tartrate_record()
    before = deepcopy(record)
    molecule, names, aliases = source_tartrate(record)
    assert Chem.GetFormalCharge(molecule) == -2
    assert molecule.GetNumAtoms() == 14
    assert len(set(names)) == len(names)
    assert {row["parent"] for row in aliases} == {"C2", "C3", "O2", "O3"}
    assert record == before
    assert molecule.GetConformer().GetPositions().tolist() == [atom["xyz"] for atom in record["atoms"]]


def test_tartrate_missing_heavy_atom_does_not_silently_rebuild(tmp_path):
    record = tartrate_record()
    record["atoms"].pop(0)
    with pytest.raises(ValueError, match="complete"):
        source_tartrate(record)


def test_tartrate_missing_alcohol_hydrogen_is_not_deprotonation_guess():
    record = tartrate_record()
    record["atoms"].pop()
    with pytest.raises(ValueError, match="explicit complete"):
        source_tartrate(record)


def test_tartrate_bad_geometry_is_not_fixed_by_renaming():
    record = tartrate_record()
    record["atoms"][0]["xyz"] = [0., 0., 0.]
    with pytest.raises(ValueError, match="geometry"):
        source_tartrate(record)


def test_tartrate_alternate_location_is_not_silently_selected():
    record = tartrate_record()
    record["atoms"][0]["altloc"] = "A"
    with pytest.raises(ValueError, match="resolved"):
        source_tartrate(record)


def test_component_hash_change_blocks_reuse(tmp_path):
    from fmotpiefinal.core.tartrate import component_molecules
    source = tmp_path / "source.pdb"
    source.write_text("changed receptor\n")
    source.with_suffix(".components.json").write_text(json.dumps({"prepared_sha256": "0" * 64}))
    with pytest.raises(ValueError, match="changed"):
        component_molecules(source)
