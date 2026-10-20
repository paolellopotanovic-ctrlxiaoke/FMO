#!/usr/bin/env python
"""Run a complete protein–ligand target only after label-free intake and geometry gates."""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
import time
from pathlib import Path

from .geometry_audit import audit as audit_geometry
from .target_intake import digest


ROOT = Path.cwd()
RUNNER = ("-m", "fmotpiefinal.pipeline.target_stages")


def prerequisites(geometry_root: Path) -> tuple[dict, Path, Path, list[str], str]:
    gate_path = geometry_root / "geometry_gate.json"
    gate = json.loads(gate_path.read_text())
    if gate.get("schema") != "fmo-new-target-geometry-gate-v1" or (
            gate.get("status") not in ("STRUCTURE_TREATMENT_COMPLETE", "K10_REPAIR_REQUIRED") or
            gate.get("ranking_eligible") is not False or
            (gate.get("status") == "STRUCTURE_TREATMENT_COMPLETE" and
             gate.get("n_residual_severe_contacts") != 0)):
        raise ValueError("whole-target structure geometry not ready for production")
    cohort = Path(gate["treated_cohort"])
    with cohort.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != gate["n_ligands"] or not rows or (
            len({row["protein_pdb"] for row in rows}) != 1 or
            len({row["ligand_sdf"] for row in rows}) != 1 or
            len({row["protein_chain"] for row in rows}) != 1 or
            len({row["target"] for row in rows}) != 1 or rows[0]["target"] != gate["target"]):
        raise ValueError("treated cohort no longer matches the target-level gate")
    protein = Path(rows[0]["protein_pdb"])
    ligands = Path(rows[0]["ligand_sdf"])
    if digest(protein) != gate["treated_protein_sha256"] or (
            digest(ligands) != gate["ligand_sdf_sha256"]):
        raise ValueError("treated receptor or ligand poses changed after structure audit")
    names = [row["ligand_name"] for row in rows]
    if len(set(names)) != len(names):
        raise ValueError("duplicate new-target ligand names")
    preflight = json.loads((geometry_root / "treated_structure_preflight.json").read_text())
    if preflight.get("ready_for_full_target_run") is not True or preflight.get("n") != len(rows):
        raise ValueError("post-water structure preflight not complete")
    return gate, protein, ligands, names, rows[0]["protein_chain"]


def enforce_prep_geometry(report: dict, stage: str, allow_nonwater_pre_k10: bool) -> None:
    severe = [contact for case in report["cases"] for contact in case["severe_contacts"]]
    if any(contact["receptor_water"] for contact in severe):
        raise ValueError(f"{stage} introduced unexpected severe receptor-water clash")
    if severe and not allow_nonwater_pre_k10:
        raise ValueError(f"{stage} has severe nonwater clash outside registered k10 repair")


def copy_if_unchanged(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        if digest(source) != digest(destination):
            raise ValueError(f"frozen source file changed: {destination}")
    else:
        shutil.copyfile(source, destination)


def select_receptor_chains(source: Path, destination: Path, chain_spec: str) -> None:
    selected = {chain.strip() for chain in chain_spec.split(",") if chain.strip()}
    if not selected:
        raise ValueError("no receptor chains registered in geometry cohort")
    atoms = []
    for line in source.read_text().splitlines(keepends=True):
        if not line.startswith(("ATOM  ", "HETATM")) or line[21].strip() not in selected:
            continue
        if line.startswith("ATOM  ") or line[17:20].strip() in ("HOH", "WAT", "H2O"):
            atoms.append(line)
    if not atoms or not any(line.startswith("ATOM  ") for line in atoms):
        raise ValueError("selected receptor chains contain no protein atoms")
    contents = "".join(atoms) + "END\n"
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        if destination.read_text() != contents:
            raise ValueError("selected receptor atom set changed")
    else:
        destination.write_text(contents)


def reset_prepared_hydrogens(source: Path, destination: Path) -> dict:
    lines = source.read_text().splitlines(keepends=True)
    protein = [line for line in lines if line.startswith(("ATOM  ", "HETATM"))]
    heavy = [line for line in protein if line[76:78].strip().upper() not in ("H", "D")]
    if not heavy or len(heavy) > len(protein) or not all(line[76:78].strip() for line in protein):
        raise ValueError("PDBFixer output has no valid heavy atoms or missing element fields")
    contents = "".join(heavy) + "END\n"
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        if destination.read_text() != contents:
            raise ValueError("PDB2PQR heavy-only receptor fingerprint changed")
    else:
        destination.write_text(contents)
    return {"schema": "fmo-pdbfixer-pdb2pqr-hydrogen-reset-v1",
            "pdbfixer_input_sha256": digest(source), "pdb2pqr_heavy_input_sha256": digest(destination),
            "n_heavy_atoms_preserved": len(heavy), "n_hydrogen_atoms_removed": len(protein) - len(heavy),
            "purpose": "let PDB2PQR/PROPKA place all hydrogens; do not retain PDBFixer HXT"}


def prepare_commands(target: str, count: int, source: Path, output: Path,
                     gamess: Path, library: Path, params: Path,
                     workers: int, ncores: int, timeout: int,
                     repair_peptide_numbering: bool = False,
                     pocket_heavy_only: bool = False,
                     preserve_input_water_hydrogens: bool = False,
                     water_source: Path | None = None) -> list[list[str]]:
    base = [sys.executable, *RUNNER, "--source", str(source), "--set-name", "newtarget",
            "--target-name", target, "--expected", str(count), "--output", str(output),
            "--workers", str(workers),
            "--ncores", str(ncores), "--timeout", str(timeout), "--merge-disulfides"]
    if repair_peptide_numbering:
        base.append("--repair-peptide-numbering")
    if pocket_heavy_only:
        base.append("--pocket-heavy-only")
    if preserve_input_water_hydrogens:
        if water_source is None:
            raise ValueError("source-water preservation requires the frozen treated receptor")
        base += ["--reference-water-hydrogens", str(water_source)]
    return [[*base, "--stage", stage] + (["--gamess-dir", str(gamess), "--library-dir",
               str(library), "--dftb-parameter-dir", str(params)] if stage == "score" else [])
            for stage in ("prepare", "minimize", "score")]


def execute(geometry_root: Path, output: Path, gamess: Path, library: Path, params: Path,
            workers: int, ncores: int, timeout: int, *, dry_run: bool = False,
            repair_peptide_numbering: bool = False,
            pocket_heavy_only: bool = False,
            preserve_input_water_hydrogens: bool = False) -> dict:
    geometry_root = geometry_root.resolve()
    output = output.resolve()
    gamess = (ROOT / gamess).resolve() if not gamess.is_absolute() else gamess.resolve()
    library = (ROOT / library).resolve() if not library.is_absolute() else library.resolve()
    params = (ROOT / params).resolve() if not params.is_absolute() else params.resolve()
    gate, protein, poses, names, chains = prerequisites(geometry_root)
    if min(workers, ncores, timeout) < 1:
        raise ValueError("workers, cores and timeout must be positive")
    if not dry_run and any(not folder.is_dir() for folder in (gamess, library, params)):
        raise ValueError("GAMESS, library and parameter directories must exist")
    source_dir = output / "source"
    source_protein = source_dir / "newtarget_protein.pdb"
    source_ligands = source_dir / "newtarget_ligands.sdf"
    raw_water_treated = output / "input_target_water_treated.pdb"
    repaired_protein = output / "pdbfixer_receptor.pdb"
    run_dir = output / "full_target_fmo"
    commands = prepare_commands(gate["target"], len(names), source_dir, run_dir,
                                gamess, library, params,
                                workers, ncores, timeout, repair_peptide_numbering,
                                pocket_heavy_only, preserve_input_water_hydrogens,
                                raw_water_treated)
    plan = {"schema": "fmo-new-target-full-execution-v1", "target": gate["target"],
            "n_ligands": len(names), "names": names, "geometry_gate_sha256": digest(
                geometry_root / "geometry_gate.json"), "input_sha256": {
                "protein": digest(protein), "ligands": digest(poses)},
            "ranking_eligible": False, "status": "PLAN_ONLY" if dry_run else "RUNNING",
            "repair_peptide_numbering": repair_peptide_numbering,
            "pocket_heavy_only": pocket_heavy_only,
            "preserve_input_water_hydrogens": preserve_input_water_hydrogens,
            "runner_stages": commands, "stage_completed": []}
    if dry_run:
        return plan
    output.mkdir(parents=True, exist_ok=True)
    report_file = output / "run_status.json"

    def checkpoint() -> None:
        report_file.write_text(json.dumps(plan, indent=2) + "\n")

    select_receptor_chains(protein, raw_water_treated, chains)
    copy_if_unchanged(poses, source_ligands)
    checkpoint()
    try:
        from .structure_repair import prepare

        if repaired_protein.exists():
            repaired_audit = output / "pdbfixer_audit.json"
            if not repaired_audit.exists():
                raise ValueError("PDBFixer receptor exists without its audit")
            prior = json.loads(repaired_audit.read_text())
            if prior["source"] != str(raw_water_treated) or prior["output"] != str(repaired_protein):
                raise ValueError("PDBFixer source/target changed")
        else:
            prepare(raw_water_treated, repaired_protein, output / "pdbfixer_audit.json", 7.4)
        reviewed = audit_geometry(repaired_protein, source_ligands)
        (output / "pdbfixer_geometry.json").write_text(json.dumps(reviewed, indent=2) + "\n")
        if len(reviewed["cases"]) != len(names):
            raise ValueError("PDBFixer geometry coverage incomplete")
        enforce_prep_geometry(reviewed, "PDBFixer", True)
        hydrogen_reset = reset_prepared_hydrogens(repaired_protein, source_protein)
        (output / "hydrogen_reset.json").write_text(json.dumps(hydrogen_reset, indent=2) + "\n")
        plan["stage_completed"].append("pdbfixer+geometry")
        checkpoint()
        for stage, command in zip(("prepare", "minimize", "score"), commands, strict=True):
            from fmotpiefinal.core.stage_execution import execute_stage

            stage_started = time.perf_counter()
            execution_mode = execute_stage(command, ROOT)
            plan.setdefault("stage_execution", []).append({
                "stage": stage, "mode": execution_mode,
                "wall_seconds": time.perf_counter() - stage_started,
            })
            if stage == "prepare":
                prepared = run_dir / f"{gate['target'].lower()}_prepared.pdb"
                from fmotpiefinal.core.component_inventory import inventory
                special = [row for row in inventory(protein)["residues"] if row["kind"] not in {"standard_residue", "water"}]
                if special and not (prepared.with_suffix(".components.json").exists() or prepared.with_suffix(".model.json").exists()):
                    from fmotpiefinal.core.special_dispatch import prepare_components
                    retained = output / "retained_components" / prepared.name
                    prepare_components(protein, prepared, retained)
                    shutil.copyfile(retained, prepared)
                    for suffix in (".components.json", ".components.sdf", ".model.json", ".dispatch.json"):
                        if retained.with_suffix(suffix).exists():
                            shutil.copyfile(retained.with_suffix(suffix), prepared.with_suffix(suffix))
                    from scipy.spatial import cKDTree
                    from rdkit import Chem
                    import numpy as np
                    poses = list(Chem.SDMolSupplier(str(run_dir / "ligands_no_properties.sdf"), removeHs=False))
                    tree = cKDTree(np.concatenate([molecule.GetConformer().GetPositions() for molecule in poses]))
                    points = [atom["xyz"] for row in special for atom in row["atoms"]]
                    if float(min(tree.query(points)[0])) <= 5.0:
                        raise ValueError("retained special component intersects QM pocket; registered in-pocket QM model required")
                from fmotpiefinal.core.graph_preparation import prepare_graph_receptor
                graph_root = output / "graph_preparation"
                if graph_root.exists():
                    graph_report = json.loads((graph_root / "graph_preparation.json").read_text())
                    if (graph_report.get("status") not in {"REPAIRED_AND_VERIFIED", "UNCHANGED_VERIFIED"}
                            or digest(prepared) != graph_report.get("output_sha256")):
                        raise ValueError("graph-prepared receptor registration changed; use a fresh run directory")
                else:
                    candidate, graph_report = prepare_graph_receptor(prepared, graph_root,
                                                                     source_water=raw_water_treated)
                    if graph_report["coordinates_modified"]:
                        copy_if_unchanged(prepared, graph_root / "original_prepared.pdb")
                        shutil.copyfile(candidate, prepared)
                        for suffix in (".components.json", ".components.sdf", ".model.json"):
                            sidecar = candidate.with_suffix(suffix)
                            if sidecar.exists():
                                shutil.copyfile(sidecar, prepared.with_suffix(suffix))
                plan["graph_preparation_sha256"] = digest(graph_root / "graph_preparation.json")
                plan["stage_completed"].append("chemical_graph+repair+reaudit")
                prepared_geometry = audit_geometry(prepared, run_dir / "ligands_no_properties.sdf")
                (output / "pdb2pqr_geometry.json").write_text(json.dumps(prepared_geometry, indent=2) + "\n")
                enforce_prep_geometry(prepared_geometry, "PDB2PQR", True)
            if stage == "minimize":
                relaxed = audit_geometry(run_dir / f"{gate['target'].lower()}_prepared.pdb",
                                         run_dir / "ligands_no_properties.sdf",
                                         run_dir / "minimized")
                (output / "post_k10_geometry.json").write_text(json.dumps(relaxed, indent=2) + "\n")
                if len(relaxed["cases"]) != len(names) or any(
                        case["minimized_pose"] is None or
                        case["minimized_pose"]["n_severe_contacts"] for case in relaxed["cases"]):
                    raise ValueError("missing or severe post-minimization pose; do not score")
                from .pocket_audit import audit_run

                pocket_report = output / "pocket_selection_audit.json"
                audit_run(run_dir, pocket_report, names)
                plan["pocket_selection_audit_sha256"] = digest(pocket_report)
            plan["stage_completed"].append(stage)
            checkpoint()
        predictions = json.loads((run_dir / "blind_predictions.json").read_text())
        if (predictions.get("n_requested") != len(names) or
                predictions.get("n_succeeded") != len(names) or
                [row["ligand"] for row in predictions["predictions"]] != names):
            raise ValueError("GAMESS FMO output does not cover the full new-target cohort")
        plan["status"] = "BLIND_FMO_COMPLETE_REQUIRES_EXTERNAL_VALIDATION"
        plan["ranking_eligible"] = True
        plan["blind_predictions_sha256"] = digest(run_dir / "blind_predictions.json")
        checkpoint()
    except Exception as error:
        plan["status"] = "REVIEW_INCOMPLETE"
        plan["ranking_eligible"] = False
        plan["failure"] = f"{type(error).__name__}: {error}"
        checkpoint()
        raise
    return plan


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--gamess-dir", type=Path, required=True)
    parser.add_argument("--library-dir", type=Path, required=True)
    parser.add_argument("--dftb-parameter-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--ncores", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=3600)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--repair-peptide-numbering", action="store_true",
                        help="use one audited target-wide HOP bond reconstruction policy")
    parser.add_argument("--pocket-heavy-only", action="store_true",
                        help="use target-wide heavy-atom pocket selection, not H contacts")
    args = parser.parse_args()
    print(json.dumps(execute(args.geometry_dir, args.output, args.gamess_dir,
                             args.library_dir, args.dftb_parameter_dir,
                             args.workers, args.ncores, args.timeout, dry_run=args.dry_run,
                             repair_peptide_numbering=args.repair_peptide_numbering,
                             pocket_heavy_only=args.pocket_heavy_only), indent=2))


if __name__ == "__main__":
    main()
