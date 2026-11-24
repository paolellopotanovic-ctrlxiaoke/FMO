"""Bounded, target-wide repair of preparation-induced amide deformation."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from rdkit import Chem

from fmotpiefinal.core.bound_geometry import load_bound_ligand
from fmotpiefinal.core.sidechain_motion import add_displaced_clash_seeds
from fmotpiefinal.core.sidechain_network import extend_one_contact_shell
from fmotpiefinal.core.sidechain_response import decide_target, diagnose_case
from fmotpiefinal.core.structure import Protein, residue_key
from fmotpiefinal.sidechains import audit_sidechains


VERSION = "induced-amide-response-v1"


def write_frozen(path, record):
    path = Path(path)
    text = json.dumps(record, indent=2, allow_nan=False) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.read_text() != text:
        raise ValueError(f"structure-response provenance changed: {path}")
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text)
    temporary.replace(path)


def apply_response(protein_path, sdf, complexes, protein, pocket, ligands, args, root, rows, signature):
    from fmotpiefinal.core.fragmentation import make_facio_like_fragments
    from fmotpiefinal.core.minimization import minimize_complex
    from .target_stages import (directory_key, receptor_at_minimized_geometry, score_one,
                                validate_retained_component_region)

    names = [row["ligand"] for row in rows]
    molecules = {molecule.GetProp("_Name").strip(): molecule for molecule in
                 Chem.SDMolSupplier(str(sdf), removeHs=False) if molecule is not None}
    if set(molecules) != set(names) or len(names) != len(set(names)):
        raise ValueError("complete unique response inventory required")
    bound = {name: np.asarray([atom.xyz for atom in ligand.atoms])
             for name, ligand in zip(names, ligands, strict=True)}
    receptors = {name: receptor_at_minimized_geometry(protein, Path(complexes[name])) for name in names}
    cases = {name: diagnose_case(molecules[name], molecules[name].GetConformer().GetPositions(),
             bound[name], bound[name], receptors[name], receptors[name]) for name in names}
    folder = Path(root) / "sidechain_response"
    source_paths = [Path(protein_path), Path(sdf), *(Path(complexes[name]) for name in names)]
    source_paths += [Path(protein_path).with_suffix(suffix) for suffix in
                     (".components.json", ".components.sdf", ".model.json")
                     if Path(protein_path).with_suffix(suffix).exists()]
    provenance = {"schema": VERSION, "names": names, "labels_used": False,
                  "source_sha256": {str(path.resolve()): hashlib.sha256(path.read_bytes()).hexdigest()
                                    for path in source_paths},
                  "baseline_predictions": rows}
    write_frozen(folder / "baseline.json", provenance)
    cached_path = folder / "decision.json"
    if cached_path.exists():
        cached = json.loads(cached_path.read_text())
        for path, expected in cached.get("candidate_sha256", {}).items():
            if hashlib.sha256(Path(path).read_bytes()).hexdigest() != expected:
                raise ValueError("cached response evidence changed")
        return cached["candidate_predictions"] if cached.get("accepted") else rows
    if not any(case["baseline_new_amide_flags"] for case in cases.values()):
        decision = decide_target(cases, names, geometry_verified=True, qm_verified=True)
        decision.update(production_default_changed=True, schema=VERSION, candidate_executed=False)
        write_frozen(folder / "decision.json", decision)
        return rows
    seed = audit_sidechains(protein_path, sdf, complexes)
    graph = extend_one_contact_shell(protein, add_displaced_clash_seeds(protein, molecules, bound, seed), n_hops=2)
    write_frozen(folder / "graph.json", graph)
    selected = tuple(map(tuple, graph["selected_residues"]))
    candidates = {}
    candidate_rows = []
    geometry_verified = False
    qm_verified = False
    error = None
    try:
        if not selected:
            raise ValueError("induced amide deformation has no safe graph-selected sidechains")
        keys = {residue_key(atom) for atom in pocket.atoms}
        for name in names:
            minimized = minimize_complex(protein_path, sdf, name, folder / "minimized" / directory_key(name),
                restraint_k_kj_mol_a2=10.0, move_protein=False,
                local_flexible_residues=selected, receptor_restraint_k_kj_mol_a2=2.5)
            receptor = receptor_at_minimized_geometry(protein, minimized.complex_pdb)
            ligand = load_bound_ligand(sdf, name, minimized.complex_pdb)
            validate_retained_component_region(protein_path, receptor, [ligand], getattr(args, "cutoff", 5.0))
            points = np.asarray([atom.xyz for atom in ligand.atoms])
            cases[name] = diagnose_case(molecules[name], molecules[name].GetConformer().GetPositions(),
                bound[name], points, receptors[name], receptor)
            candidate_pocket = Protein(tuple(atom for atom in receptor.atoms if residue_key(atom) in keys))
            prepared, fragments, metadata = make_facio_like_fragments(candidate_pocket, ligand,
                include_waters=True, reference_protein=receptor, use_hop=True,
                merge_disulfides=args.merge_disulfides, repair_peptide_numbering=args.repair_peptide_numbering)
            current = (prepared, fragments[:-1], metadata["hop_pairs_gamess_indices"])
            identities = lambda atoms: [(atom.element, atom.chain, atom.residue_number,
                                        atom.residue_name, atom.atom_name) for atom in atoms]
            if (identities(prepared.atoms) != identities(signature[0].atoms) or
                    current[1:] != signature[1:]):
                raise ValueError("response changed frozen fragment/charge/cap topology")
            candidates[name] = (receptor, candidate_pocket, ligand, current)
        geometry_verified = True
        preliminary = decide_target(cases, names, geometry_verified=True, qm_verified=True)
        if preliminary["action"] == "ACCEPT_LOCAL_REPAIR_RESEARCH_CONTROL":
            for name in names:
                receptor, candidate_pocket, ligand, current = candidates[name]
                candidate_rows.append(score_one(name, receptor, candidate_pocket, ligand, args, folder, current))
            qm_verified = len(candidate_rows) == len(names)
    except Exception as failure:
        error = f"{type(failure).__name__}: {failure}"
    decision = decide_target(cases, names, geometry_verified=geometry_verified, qm_verified=qm_verified)
    accepted = decision["action"] == "ACCEPT_LOCAL_REPAIR_RESEARCH_CONTROL"
    decision.update(schema=VERSION, production_default_changed=True, candidate_executed=True,
                    accepted=accepted, error=error, cases=cases, candidate_predictions=candidate_rows)
    evidence_paths = [path for path in folder.rglob("*") if path.is_file() and
                      path.name in {"minimized_complex.pdb", "minimization.json", "input.inp", "gamess.log", "result.json"}]
    decision["candidate_sha256"] = {str(path.resolve()): hashlib.sha256(path.read_bytes()).hexdigest()
                                    for path in evidence_paths}
    if accepted:
        decision["action"] = "ACCEPT_LOCAL_REPAIR"
    write_frozen(folder / "decision.json", decision)
    return candidate_rows if accepted else rows
