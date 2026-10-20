#!/usr/bin/env python
"""Score cross-partition PIEDA; never interpret as a binding free energy."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path


from fmotpiefinal.core.interface_pieda import cross_interface
from fmotpiefinal.core.parser import find_selected_tie, parse_pie


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def score(log_path: Path, partition_path: Path, result_path: Path | None = None) -> dict:
    partition = json.loads(partition_path.read_text(encoding="utf-8"))
    if set(partition) != {"fragment_ids_a", "fragment_ids_b", "interface_type"}:
        raise ValueError("partition requires exactly fragment_ids_a, fragment_ids_b, interface_type")
    if partition["interface_type"] not in ("protein_ligand", "protein_protein",
                                            "protein_rna", "rna_ligand", "custom"):
        raise ValueError("unsupported interface_type")
    first, second = partition["fragment_ids_a"], partition["fragment_ids_b"]
    if (not isinstance(first, list) or not isinstance(second, list) or
            len(set(first)) != len(first) or len(set(second)) != len(second)):
        raise ValueError("fragment groups must be nonduplicate lists")
    with log_path.open("rb") as handle:
        handle.seek(max(0, log_path.stat().st_size - 65536))
        if b"EXECUTION OF GAMESS TERMINATED NORMALLY" not in handle.read().upper():
            raise ValueError("GAMESS job did not terminate normally")
    interface = cross_interface(parse_pie(log_path), first, second)
    if len(second) == 1:
        selected_tie = find_selected_tie(log_path, second[0])
        if selected_tie is not None:
            if not math.isclose(interface["interaction_kcal_mol"], selected_tie,
                                rel_tol=0, abs_tol=max(0.02, .002 * interface["n_cross_pairs"])):
                if result_path is not None or partition["interface_type"] == "protein_ligand":
                    raise ValueError("cross-interface sum differs from GAMESS selected TIE")
            else:
                interface["gamess_selected_tie_kcal_mol"] = selected_tie
    result_hash = None
    if result_path is not None:
        result = json.loads(result_path.read_text(encoding="utf-8"))
        if result.get("status") != "SUCCEEDED" or result.get("ligand_fragment") != second[0] or (
                len(second) != 1 or not math.isclose(float(result["tpie_kcal_mol"]),
                interface["interaction_kcal_mol"], rel_tol=0,
                abs_tol=max(0.02, .002 * interface["n_cross_pairs"]))):
            raise ValueError("archived result disagrees with cross-interface score")
        result_hash = digest(result_path)
    return {"schema": "fmo-cross-interface-pieda-v1", "interface_type": partition["interface_type"],
            "log_sha256": digest(log_path), "partition_sha256": digest(partition_path),
            "result_sha256": result_hash, "rank_only_not_binding_free_energy": True,
            **interface}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--partition", type=Path, required=True)
    parser.add_argument("--result", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = score(args.log, args.partition, args.result)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: report[key] for key in ("interface_type", "n_cross_pairs", "interaction_kcal_mol")}))


if __name__ == "__main__":
    main()
