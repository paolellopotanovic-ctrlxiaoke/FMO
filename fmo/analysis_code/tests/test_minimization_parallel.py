from concurrent.futures import ThreadPoolExecutor

import pytest


def test_same_named_receptors_use_distinct_openmm_temp_files(tmp_path):
    pytest.importorskip("openmm")
    pytest.importorskip("openff")
    pytest.importorskip("openmmforcefields")
    from fmotpiefinal.core.minimization import _water_compatible_pdb

    first, second = tmp_path / "first/protein.pdb", tmp_path / "second/protein.pdb"
    first.parent.mkdir()
    second.parent.mkdir()
    first.write_text("REMARK first receptor\nEND\n")
    second.write_text("REMARK second receptor\nEND\n")
    with ThreadPoolExecutor(max_workers=2) as pool:
        outputs = list(pool.map(_water_compatible_pdb, [first, second]))
    assert outputs[0] != outputs[1]
    assert "first receptor" in outputs[0].read_text()
    assert "second receptor" in outputs[1].read_text()
    assert _water_compatible_pdb(first) == outputs[0]


def test_shared_receptor_conversion_is_atomic_for_parallel_readers(tmp_path):
    pytest.importorskip("openmm")
    pytest.importorskip("openff")
    pytest.importorskip("openmmforcefields")
    from fmotpiefinal.core.minimization import _water_compatible_pdb
    source = tmp_path / "shared.pdb"
    content = "REMARK complete receptor\n" * 10000 + "END\n"
    source.write_text(content)

    def convert_and_read(index):
        return _water_compatible_pdb(source).read_text()

    with ThreadPoolExecutor(max_workers=8) as pool:
        assert all(text == content for text in pool.map(convert_and_read, range(32)))
