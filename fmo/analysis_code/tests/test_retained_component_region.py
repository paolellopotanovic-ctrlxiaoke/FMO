import json
from types import SimpleNamespace

import pytest

from fmotpiefinal.core.component_inventory import sha
from fmotpiefinal.core.structure import Atom, Protein
from fmotpiefinal.pipeline.target_stages import validate_retained_component_region


def inputs(tmp_path, distance):
    source = tmp_path / "receptor.pdb"
    source.write_text("registered receptor\n")
    source.with_suffix(".components.json").write_text(json.dumps({"prepared_sha256": sha(source), "models": [{"id": "A:TLA:502"}]}))
    atom = Atom(0, "C", (distance, 0., 0.), "A", "TLA", "502", "HETATM", "C1")
    protein = Protein((atom,))
    ligands = [SimpleNamespace(atoms=[SimpleNamespace(xyz=(0., 0., 0.))])]
    return source, protein, ligands


def test_component_inside_qm_region_is_rejected_even_when_pocket_filter_would_drop_hetatm(tmp_path):
    source, protein, ligands = inputs(tmp_path, 2.)
    with pytest.raises(ValueError, match="special QM model required"):
        validate_retained_component_region(source, protein, ligands, 5.)


def test_remote_retained_component_has_explicit_region_evidence(tmp_path):
    source, protein, ligands = inputs(tmp_path, 17.)
    report = validate_retained_component_region(source, protein, ligands, 5.)
    assert report["component_ids"] == ["A:TLA:502"]
    assert report["minimum_contact_angstrom"] == 17.


def test_registered_component_cannot_disappear_before_qm(tmp_path):
    source, protein, ligands = inputs(tmp_path, 17.)
    with pytest.raises(ValueError, match="missing"):
        validate_retained_component_region(source, Protein(()), ligands, 5.)


@pytest.mark.parametrize('distance', [2., 17.])
def test_retained_terminal_cap_is_not_silently_dropped_by_qm_filter(tmp_path, distance):
    source, protein, ligands = inputs(tmp_path, distance)
    source.with_suffix('.components.json').unlink()
    source.with_suffix('.model.json').write_text(json.dumps({
        'schema': 'source-terminal-cap-retention-v1', 'prepared_sha256': sha(source),
        'retained': [{'raw_id': '_:UNK:890', 'prepared_id': 'A:ACE:2'}]}))
    protein = Protein((Atom(0, 'C', (distance, 0., 0.), 'A', 'ACE', '2', 'ATOM', 'C'),))
    if distance < 5.:
        with pytest.raises(ValueError, match='special QM model required'):
            validate_retained_component_region(source, protein, ligands, 5.)
    else:
        assert validate_retained_component_region(source, protein, ligands, 5.)['component_ids'] == ['A:ACE:2']


def test_blank_chain_ion_identity_survives_structure_reader(tmp_path):
    source, protein, ligands = inputs(tmp_path, 17.)
    source.with_suffix('.components.json').unlink()
    source.with_suffix('.model.json').write_text(json.dumps({
        'schema': 'retained-halide-model-v1', 'prepared_sha256': sha(source),
        'matched_ion_templates': {'_:CL:2': {'template': 'CL', 'charge': -1.}}}))
    protein = Protein((Atom(0, 'CL', (17., 0., 0.), '', 'CL', '2', 'HETATM', 'CL'),))
    assert validate_retained_component_region(source, protein, ligands, 5.)['component_ids'] == ['_:CL:2']
