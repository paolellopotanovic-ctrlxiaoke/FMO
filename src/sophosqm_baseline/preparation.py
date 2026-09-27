"""Open, reproducible protein preparation used before FMO fragmentation."""
from __future__ import annotations

from pathlib import Path
import json
import shutil
import subprocess


def _filter_protein_input(source: Path, target: Path, *, keep_water: bool) -> int:
    """Write PDB2PQR input records and return retained water-oxygen count."""
    water_oxygens = 0
    with source.open(encoding="utf-8") as handle, target.open("w", encoding="utf-8") as output:
        for line in handle:
            record = line[:6].strip()
            residue = line[17:20].strip().upper()
            atom_name = line[12:16].strip().upper()
            is_water_oxygen = (
                record == "HETATM"
                and residue in {"HOH", "WAT", "TIP"}
                and atom_name in {"O", "OW", "OX"}
            )
            water_oxygens += int(is_water_oxygen)
            if record == "ATOM" or (
                record == "HETATM"
                and residue in {"HOH", "WAT", "TIP"}
                and keep_water
            ):
                output.write(line)
        output.write("END\n")
    return water_oxygens


def prepare_protein_pdb2pqr(input_pdb: str | Path, output_pdb: str | Path,
                            *, keep_water: bool = True) -> dict:
    """Run the installed PDB2PQR backend and persist auditable provenance.

    Non-protein ligand HETATM records are excluded before preparation; solvent
    waters are retained and the ligand remains sourced from SDF.
    """
    source = Path(input_pdb)
    output = Path(output_pdb)
    output.parent.mkdir(parents=True, exist_ok=True)
    filtered = output.with_suffix(".protein_input.pdb")
    input_water_oxygens = _filter_protein_input(source, filtered, keep_water=keep_water)
    executable = shutil.which("pdb2pqr")
    if executable is None:
        raise RuntimeError("pdb2pqr is required for reproducible protein preparation")
    pqr = output.with_suffix(".pqr")
    command = [
        executable,
        "--ff=AMBER",
        "--titration-state-method=propka",
        "--with-ph=7.4",
        "--pdb-output",
        str(output),
    ]
    if not keep_water:
        command.append("--drop-water")
    command.extend([str(filtered), str(pqr)])
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    log_path = output.with_suffix(".pdb2pqr.log")
    log_path.write_text(completed.stdout + "\n" + completed.stderr, encoding="utf-8")
    if completed.returncode != 0 or not output.exists():
        raise RuntimeError(f"PDB2PQR failed; inspect {log_path}")
    metadata = {
        "backend": "pdb2pqr",
        "executable": executable,
        "forcefield": "AMBER",
        "titration_state_method": "propka",
        "ph": 7.4,
        "keep_water": keep_water,
        "explicit_water_policy": "retained" if keep_water else "removed_before_pdb2pqr_propka",
        "input_water_oxygens": input_water_oxygens,
        "retained_water_oxygens": input_water_oxygens if keep_water else 0,
        "filtered_input": str(filtered),
        "prepared_pdb": str(output),
        "pqr": str(pqr),
        "log": str(log_path),
        "returncode": completed.returncode,
    }
    output.with_suffix(".preparation.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata
