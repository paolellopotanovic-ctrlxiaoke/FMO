import csv
import hashlib
from pathlib import Path

from sophosqm_baseline.cohort import build_cohort


def test_reconstructs_exact_87_compound_four_target_cohort():
    structure_root = Path("vendor/public_binding_free_energy_benchmark/fep_benchmark_inputs/structure_inputs/jacs_set")
    if structure_root.is_dir():
        import pytest
        pytest.importorskip("rdkit")
        records = build_cohort(structure_root)
        counts = {target: sum(record.target == target for record in records)
                  for target in ["CDK2", "JNK1", "P38", "TYK2"]}
        assert len(records) == 87
        assert counts == {"CDK2": 16, "JNK1": 21, "P38": 34, "TYK2": 16}
        assert all(record.formal_charge == 0 for record in records)
        jnk_duplicates = [record for record in records if record.target == "JNK1" and record.source_record_count > 1]
        assert len(jnk_duplicates) == 17
        assert all("flip" not in record.ligand_name.lower() for record in jnk_duplicates)
        return

    frozen = Path("data/fmo87_cohort.csv")
    with frozen.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    counts = {target: sum(row["target"] == target for row in rows)
              for target in ["CDK2", "JNK1", "P38", "TYK2"]}
    assert len(rows) == 87
    assert counts == {"CDK2": 16, "JNK1": 21, "P38": 34, "TYK2": 16}
    assert all(int(row["formal_charge"]) == 0 for row in rows)
    assert hashlib.sha256(frozen.read_bytes()).hexdigest() == "19edfb8e2ccf3ab1029d848e4062d70cb35532d2b4fee63674f160b2a292997b"
