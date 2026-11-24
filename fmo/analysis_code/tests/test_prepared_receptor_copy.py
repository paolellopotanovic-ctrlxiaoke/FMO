import pytest

from fmotpiefinal.pipeline.target_stages import copy_prepared_receptor


def test_prepared_override_carries_registered_component_model(tmp_path):
    source = tmp_path / "source.pdb"
    source.write_text("receptor\n")
    source.with_suffix(".components.json").write_text('{"model": "pinned"}')
    source.with_suffix(".components.sdf").write_text("explicit component coordinates")
    destination = tmp_path / "copy.pdb"
    copy_prepared_receptor(source, destination)
    assert destination.with_suffix(".components.json").read_bytes() == source.with_suffix(".components.json").read_bytes()
    assert destination.with_suffix(".components.sdf").read_bytes() == source.with_suffix(".components.sdf").read_bytes()
    copy_prepared_receptor(source, destination)


def test_conflicting_component_registration_not_silently_overwritten(tmp_path):
    source, destination = tmp_path / "source.pdb", tmp_path / "copy.pdb"
    source.write_text("receptor\n")
    source.with_suffix(".model.json").write_text("original")
    copy_prepared_receptor(source, destination)
    destination.with_suffix(".model.json").write_text("different model")
    with pytest.raises(ValueError, match="sidecar changed"):
        copy_prepared_receptor(source, destination)


def test_unregistered_destination_sidecar_is_rejected(tmp_path):
    source, destination = tmp_path / "source.pdb", tmp_path / "copy.pdb"
    source.write_text("receptor\n")
    destination.with_suffix(".model.json").write_text("unregistered")
    with pytest.raises(ValueError, match="unregistered"):
        copy_prepared_receptor(source, destination)
