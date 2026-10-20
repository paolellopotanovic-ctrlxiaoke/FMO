"""Prepare auditable GAMESS jobs from Ross structure inputs."""
from __future__ import annotations

from pathlib import Path
from dataclasses import replace
import json

from .fragmentation import make_facio_like_fragments, write_fragment_manifest
from .gamess import GamessConfig, write_job_manifest
from .fmo_input import build_deck_atoms, write_fmo2_pieda_input
from .structure import bridge_one_residue_sequence_gaps, read_pdb, read_sdf, select_pocket


def prepare_job(protein_path: str | Path, ligand_sdf: str | Path, ligand_name: str,
                output_dir: str | Path, dftb_parameter_dir: str | Path,
                minimized_complex: str | Path | None = None,
                *, include_waters: bool = True) -> Path:
    protein = read_pdb(protein_path)
    ligands = read_sdf(ligand_sdf)
    ligand = next((item for item in ligands if item.name == ligand_name), None)
    if ligand is None:
        available = ", ".join(item.name for item in ligands[:10])
        raise KeyError(f"Ligand {ligand_name!r} not found; first available names: {available}")
    ligand_coordinates_source = "input SDF"
    if minimized_complex is not None:
        complex_protein = read_pdb(minimized_complex)
        complex_ligand_atoms = [atom for atom in complex_protein.atoms
                                if atom.record == "HETATM" and atom.residue_name == "UNK"]
        if len(complex_ligand_atoms) != len(ligand.atoms) or any(
            left.element.upper() != right.element.upper()
            for left, right in zip(complex_ligand_atoms, ligand.atoms, strict=True)
        ):
            raise ValueError("minimized complex ligand does not match the selected SDF atom order")
        ligand = replace(ligand, atoms=tuple(
            replace(atom, xyz=complex_atom.xyz)
            for atom, complex_atom in zip(ligand.atoms, complex_ligand_atoms, strict=True)
        ))
        ligand_coordinates_source = str(minimized_complex)
    cutoff = 5.0
    pocket = select_pocket(protein, ligand, cutoff)
    pocket, bridged_residues = bridge_one_residue_sequence_gaps(pocket, protein)
    preparation_metadata = {}
    pocket, fragments, preparation_metadata = make_facio_like_fragments(
        pocket, ligand, include_waters=include_waters, reference_protein=protein
    )
    if not dftb_parameter_dir:
        raise ValueError("FMO87 job preparation requires a 3OB-3-1 parameter directory")
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    structure_manifest = output / "structure_manifest.json"
    structure_manifest.write_text(json.dumps({
        "protein_path": str(protein_path), "ligand_sdf": str(ligand_sdf),
        "ligand_name": ligand_name, "cutoff_angstrom": cutoff,
        "protein_atoms_in_pocket": len(pocket.atoms), "ligand_atoms": len(ligand.atoms),
        "ligand_coordinates_source": ligand_coordinates_source,
        "include_explicit_waters": include_waters,
        "bridged_one_residue_sequence_gaps": bridged_residues,
    }, indent=2), encoding="utf-8")
    fragment_manifest = output / "fragment_manifest.json"
    write_fragment_manifest(fragment_manifest, fragments, {
        "protein_path": str(protein_path), "ligand_sdf": str(ligand_sdf),
        "ligand_name": ligand_name, "cutoff_angstrom": cutoff,
        "ligand_fragment_id": fragments[-1].fragment_id,
        "readiness": "GAMESS_DECK_GENERATED_FULL_ACE_NME_SHIFTED_C_ALPHA_C",
        "qm_protocol": {
            "engine_method": "FMO2-DFTB3/PIEDA",
            "solvation": "water PCM (IFMO=-1) for DFTB3",
        },
        "explicit_water_policy": "retained" if include_waters else "removed_from_fmo_pocket",
        "preparation": preparation_metadata,
        "limitations": [
            "FACIO is approximated by the documented shifted C-alpha/C-prime partition",
            "full neutral ACE/NME residues saturate pocket segment termini",
            "input PDB hydrogens are retained; PDB2PQR/PROPKA status must be recorded by the caller",
            "atom coordinates and fragment membership are serialized to a GAMESS FMO deck",
            "fragment charge/electron parity is validated at GAMESS deck generation",
        ],
    })
    input_path = output / "job.inp"
    deck_atoms = build_deck_atoms(pocket, ligand, fragments)
    fragment_bonds = [tuple(pair) for pair in preparation_metadata.get("hop_pairs_gamess_indices", [])]
    write_fmo2_pieda_input(input_path, pocket, ligand, fragments,
                           title=f"FMO87 {ligand_name}",
                           dftb_parameter_dir=dftb_parameter_dir,
                           fragment_bonds=fragment_bonds)
    return write_job_manifest(output,
                              GamessConfig("GAMESS_CONFIGURE_ME", str(input_path),
                                           method="DFTB3",
                                           basis="DFTB3/3OB-3-1"),
                              str(structure_manifest), str(fragment_manifest), {
                                  "status": "GAMESS_DECK_GENERATED_NOT_RUN",
                                  "input_path": str(input_path),
                                  "method": "FMO2-DFTB3/PIEDA",
                                  "fragment_count": len(fragments),
                                  "ligand_fragment_id": fragments[-1].fragment_id,
                              })
