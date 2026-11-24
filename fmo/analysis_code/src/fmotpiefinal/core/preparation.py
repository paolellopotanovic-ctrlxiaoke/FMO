"""Open, reproducible protein preparation used before FMO fragmentation."""
from __future__ import annotations

from pathlib import Path
import json
import shutil
import subprocess


def prepare_protein_pdb2pqr(input_pdb: str | Path, output_pdb: str | Path,
                            *, ph: float = 7.4, forcefield: str = "AMBER",
                            keep_water: bool = True) -> dict:
    """Run the installed PDB2PQR backend and persist auditable provenance.

    Non-protein ligand HETATM records are excluded before preparation; solvent
    waters are retained by default and the ligand remains sourced from SDF.
    """
    source = Path(input_pdb)
    output = Path(output_pdb)
    output.parent.mkdir(parents=True, exist_ok=True)
    filtered = output.with_suffix(".protein_input.pdb")
    with source.open(encoding="utf-8") as handle, filtered.open("w", encoding="utf-8") as target:
        for line in handle:
            record = line[:6].strip()
            residue = line[17:20].strip().upper()
            if record == "ATOM" or (keep_water and record == "HETATM" and residue in {"HOH", "WAT", "TIP"}):
                target.write(line)
        target.write("END\n")
    executable = shutil.which("pdb2pqr")
    if executable is None:
        raise RuntimeError("pdb2pqr is required for reproducible protein preparation")
    pqr = output.with_suffix(".pqr")
    command = [
        executable,
        f"--ff={forcefield}",
        "--titration-state-method=propka",
        f"--with-ph={ph}",
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
        "forcefield": forcefield,
        "titration_state_method": "propka",
        "ph": ph,
        "keep_water": keep_water,
        "filtered_input": str(filtered),
        "prepared_pdb": str(output),
        "pqr": str(pqr),
        "log": str(log_path),
        "returncode": completed.returncode,
    }
    output.with_suffix(".preparation.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata
