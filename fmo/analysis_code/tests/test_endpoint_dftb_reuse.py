import hashlib
import json
from pathlib import Path
import sys
from types import ModuleType

sys.path.insert(0, str(Path(__file__).parents[1] / "tools"))
from run_endpoint_dftb_internal import verified_cached_result


def test_reuse_checks_coordinates_deck_binary_parameters_and_log(tmp_path, monkeypatch):
    module = ModuleType("sophosqm_baseline.relative_total_energy")
    module.parse_full_energy = lambda text, input_text: float(text)
    monkeypatch.setitem(sys.modules, "sophosqm_baseline.relative_total_energy", module)
    source = tmp_path / "source"
    source.mkdir()
    parameters = tmp_path / "sk"
    parameters.mkdir()
    parameter = parameters / "C-C.skf"
    parameter.write_text("frozen SKF")
    installation = tmp_path / "gamess"
    installation.mkdir()
    binary = installation / "gamess.00.x"
    binary.write_text("frozen binary")
    old_deck = f"   C C {parameter}\n"
    new_deck = old_deck.replace(str(parameters), "/new/alias")
    checksum = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    parameter_hash = hashlib.sha256()
    parameter_hash.update(str(parameter).encode())
    parameter_hash.update(bytes.fromhex(checksum(parameter)))
    record = {"status": "SUCCEEDED", "energy_hartree": -1.,
              "input_sha256": hashlib.sha256(old_deck.encode()).hexdigest(),
              "parameters_sha256": parameter_hash.hexdigest(), "executable_sha256": checksum(binary)}
    (source / "result.json").write_text(json.dumps(record))
    log = source / "gamess.log"
    log.write_text("-1.0")
    task = {"short_parameter_path": "/new/alias", "runtime": {"gamess_dir": str(installation)},
            "reuse_candidate": {"directory": str(source), "short_parameter_path": str(parameters)}}
    assert verified_cached_result(task, new_deck)["cache_reused"]
    assert verified_cached_result(task, new_deck + "changed geometry") is None
    log.write_text("-2.0")
    assert verified_cached_result(task, new_deck) is None
    log.write_text("-1.0")
    parameter.write_text("modified SKF")
    assert verified_cached_result(task, new_deck) is None
