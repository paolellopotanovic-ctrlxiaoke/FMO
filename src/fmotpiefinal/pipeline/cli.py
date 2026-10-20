#!/usr/bin/env python
"""Unified fail-closed entry points for FMO intake, structure, and interface audits."""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

from fmotpiefinal.scoring import score_run


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="stage", required=True)
    normalize = commands.add_parser("normalize", help="loss-audited PDB/mmCIF and residue-consistent altloc selection")
    normalize.add_argument("--protein", type=Path, required=True)
    normalize.add_argument("--structure-model")
    normalize.add_argument("--output", type=Path, required=True)
    inspect = commands.add_parser("inspect", help="parse receptor and report chemistry routes without launching or approving QM")
    inspect.add_argument("--protein", type=Path, required=True)
    inspect.add_argument("--structure-model")
    inspect.add_argument("--output", type=Path, required=True)
    chemistry = commands.add_parser("chemistry-audit", help="label-free component/chemical graph diagnostics; never edits coordinates")
    chemistry.add_argument("--raw-protein", type=Path, required=True)
    chemistry.add_argument("--prepared-protein", type=Path)
    chemistry.add_argument("--dftb-parameter-dir", type=Path)
    chemistry.add_argument("--ccd-cache", type=Path)
    chemistry.add_argument("--output", type=Path, required=True)
    special = commands.add_parser("special-prepare", help="chemistry-driven retained-ion/phosphoresidue MM models; not a QM model")
    special.add_argument("--raw-protein", type=Path, required=True)
    special.add_argument("--prepared-protein", type=Path, required=True)
    special.add_argument("--output", type=Path, required=True)
    new = commands.add_parser("intake", help="inspect SMILES and optional receptor-frame poses")
    new.add_argument("--target", required=True)
    new.add_argument("--ligands", type=Path, required=True)
    new.add_argument("--protein", type=Path)
    new.add_argument("--poses", type=Path)
    new.add_argument("--component-decisions", type=Path)
    new.add_argument("--chain", default="auto",
                     help="auto: audit the unique receptor chain contacting every bound ligand")
    new.add_argument("--output", type=Path, required=True)
    structural = commands.add_parser("geometry", help="audit full-target water and clashes")
    structural.add_argument("--intake-dir", type=Path, required=True)
    structural.add_argument("--output", type=Path, required=True)
    interaction = commands.add_parser("interface", help="score explicit FMO fragment partitions")
    interaction.add_argument("--log", type=Path, required=True)
    interaction.add_argument("--partition", type=Path, required=True)
    interaction.add_argument("--result", type=Path)
    interaction.add_argument("--output", type=Path, required=True)
    general = commands.add_parser("partner-intake", help="audit pre-aligned protein/RNA complexes")
    general.add_argument("--complex", type=Path, required=True)
    general.add_argument("--interface-type", choices=("protein_protein", "protein_rna",
                                                     "rna_ligand"), required=True)
    general.add_argument("--chain-a", required=True)
    general.add_argument("--chain-b", required=True)
    general.add_argument("--ligand-sdf", type=Path)
    general.add_argument("--output", type=Path, required=True)
    for stage in ("run", "adaptive-run"):
        production = commands.add_parser(stage, help="PDBFixer→PDB2PQR→k10→FMO"
                                         + ("; label-blind registered rescue on failure" if stage == "adaptive-run" else ""))
        production.add_argument("--geometry-dir", type=Path, required=True)
        production.add_argument("--output", type=Path, required=True)
        production.add_argument("--gamess-dir", type=Path, required=True)
        production.add_argument("--library-dir", type=Path, required=True)
        production.add_argument("--dftb-parameter-dir", type=Path, required=True)
        production.add_argument("--workers", type=int, default=1)
        production.add_argument("--ncores", type=int, default=1)
        production.add_argument("--timeout", type=int, default=3600)
        production.add_argument("--dry-run", action="store_true")
        production.add_argument("--no-repair-peptide-numbering", dest="repair_peptide_numbering",
                                action="store_false", help="legacy-only: disable geometry-checked peptide connectivity")
        production.set_defaults(repair_peptide_numbering=True)
        production.add_argument("--pocket-heavy-only", action="store_true")
        production.add_argument("--preserve-input-water-hydrogens", action="store_true",
                                help="retain validated source water orientations after PDB2PQR")
    blind = commands.add_parser("blind-run", help="label-free intake→geometry→adaptive QM for any cohort")
    blind.add_argument("--target", required=True)
    blind.add_argument("--ligands", type=Path, required=True)
    blind.add_argument("--protein", type=Path, required=True)
    blind.add_argument("--poses", type=Path, required=True)
    blind.add_argument("--component-decisions", type=Path)
    blind.add_argument("--chain", default="auto")
    blind.add_argument("--structure-model", help="explicit source MODEL number; required for multi-model input")
    blind.add_argument("--output", type=Path, required=True)
    blind.add_argument("--gamess-dir", type=Path, required=True)
    blind.add_argument("--library-dir", type=Path, required=True)
    blind.add_argument("--dftb-parameter-dir", type=Path, required=True)
    blind.add_argument("--workers", type=int, default=1)
    blind.add_argument("--ncores", type=int, default=1)
    blind.add_argument("--timeout", type=int, default=3600)
    blind.add_argument("--preserve-input-water-hydrogens", action="store_true",
                       help="experimental source-water arm; default uses PDB2PQR orientations")
    blind.add_argument("--score-reference", type=Path,
                       help="previously frozen target-local reference; otherwise report local diagnostic")
    ranking = commands.add_parser("score", help="rank a validated complete blind-run output")
    ranking.add_argument("--run-dir", type=Path, required=True)
    ranking.add_argument("--score-reference", type=Path)
    args = parser.parse_args()
    if args.stage == "inspect":
        from fmotpiefinal.core.structure_input import inspect_structure
        print(json.dumps(inspect_structure(args.protein, args.output, args.structure_model)))
        return
    if args.stage == "normalize":
        from fmotpiefinal.core.structure_input import normalize_structure
        _, report = normalize_structure(args.protein, args.output, args.structure_model)
        print(json.dumps(report))
        return
    if args.stage == "chemistry-audit":
        from fmotpiefinal.core.chemical_preflight import audit_chemistry
        report = audit_chemistry(args.raw_protein, args.output, args.prepared_protein,
                                 args.dftb_parameter_dir, args.ccd_cache)
        print(json.dumps(report))
        if report["status"] != "PASS_CURRENT_CHEMISTRY_CHECKS":
            raise SystemExit(2)
        return
    if args.stage == "score":
        scored = score_run(args.run_dir, args.score_reference)
        print(json.dumps({"stage": "score", "ranking": str(args.run_dir / "ranking/ranking.json"),
                          "n": scored["n"], "softcap_defined": scored["softcap_defined"]}))
        return
    from .target_geometry import audit as audit_target_structure
    from .interface_intake import audit as audit_molecular_interface
    from .target_intake import audit as audit_intake
    from .target_workflow import execute as execute_target
    from .single_case_rescue import rescue as rescue_target, strict_failure
    from .cohort_rescue import eligible_failure, rescue_cohort
    from .interface_scoring import score as score_interface
    if args.stage == "blind-run":
        if min(args.workers, args.ncores, args.timeout) < 1:
            parser.error("workers, cores and timeout must be positive")
        root = args.output.resolve()
        from fmotpiefinal.core.structure_input import normalize_structure
        protein, normalization = normalize_structure(args.protein, root / "structure_input", args.structure_model)
        chain = normalization["chain_mapping"].get(args.chain, args.chain)
        intake = audit_intake(args.target, args.ligands, root / "intake", protein,
                              chain, args.poses, args.component_decisions)
        (root / "intake/intake.json").write_text(json.dumps(intake, indent=2) + "\n")
        geometry = audit_target_structure(root / "intake", root / "geometry")
        (root / "geometry/geometry_gate.json").write_text(json.dumps(geometry, indent=2) + "\n")
        command = [sys.executable, "-m", "fmotpiefinal.pipeline.cli", "adaptive-run",
                   "--geometry-dir", str(root / "geometry"), "--output", str(root / "qm"),
                   "--gamess-dir", str(args.gamess_dir), "--library-dir", str(args.library_dir),
                   "--dftb-parameter-dir", str(args.dftb_parameter_dir), "--workers",
                   str(args.workers), "--ncores", str(args.ncores), "--timeout", str(args.timeout)]
        if args.preserve_input_water_hydrogens:
            command.append("--preserve-input-water-hydrogens")
        with (root / "adaptive.stdout").open("w") as handle:
            completed = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT,
                                       check=False)
        if completed.returncode:
            raise SystemExit(f"REVIEW_REQUIRED: {root / 'adaptive.stdout'}")
        report = json.loads((root / "adaptive_status.json").read_text())
        if report.get("ranking_eligible") is not True:
            raise SystemExit("REVIEW_REQUIRED: no validated blind ranking")
        scored = score_run(root, args.score_reference)
        print(json.dumps({"stage": args.stage, "output": str(root / "adaptive_status.json"),
                          "status": report["status"], "ranking_eligible": True,
                          "ranking": str(root / "ranking/ranking.json"),
                          "n": scored["n"], "softcap_defined": scored["softcap_defined"]}))
        return
    if args.stage == "special-prepare":
        from fmotpiefinal.core.special_dispatch import prepare_components
        report = prepare_components(args.raw_protein, args.prepared_protein, args.output)
        print(json.dumps({"stage": args.stage, "output": str(args.output),
                          "n_retained_components": len(report["retained"]),
                          "chemistry_route": report["dispatch"]["route"],
                          "forcefield_system_created": report["forcefield_system_created"],
                          "ranking_eligible": False}))
        return
    if args.stage == "intake":
        report = audit_intake(args.target, args.ligands, args.output, args.protein,
                              args.chain, args.poses, args.component_decisions)
        destination = args.output / "intake.json"
    elif args.stage == "geometry":
        report = audit_target_structure(args.intake_dir, args.output)
        destination = args.output / "geometry_gate.json"
    elif args.stage == "interface":
        report = score_interface(args.log, args.partition, args.result)
        destination = args.output
    elif args.stage == "partner-intake":
        report = audit_molecular_interface(args.complex, args.interface_type, args.chain_a,
                                           args.chain_b, args.ligand_sdf)
        destination = args.output
    elif args.stage == "adaptive-run":
        if args.dry_run or args.pocket_heavy_only or not args.repair_peptide_numbering:
            parser.error("adaptive-run requires registered numbered default parameters")
        gate = json.loads((args.geometry_dir / "geometry_gate.json").read_text())
        name = gate["target"]
        case = args.output.parent
        if gate.get("n_ligands") == 1 and args.workers != 1:
            parser.error("single-crystal rescue requires one QM worker")
        from .target_workflow import prerequisites
        ligand_names = prerequisites(args.geometry_dir)[3]
        cohort_names = ligand_names if gate.get("n_ligands") != 1 else None
        ligand_name = ligand_names[0] if cohort_names is None else None
        prior = (strict_failure(case, ligand_name, args.output) if cohort_names is None else
                 eligible_failure(args.output, cohort_names))
        if prior is None:
            try:
                report = execute_target(args.geometry_dir, args.output, args.gamess_dir,
                                        args.library_dir, args.dftb_parameter_dir,
                                        args.workers, args.ncores, args.timeout,
                                        repair_peptide_numbering=True,
                                        preserve_input_water_hydrogens=args.preserve_input_water_hydrogens)
                (args.output / "run_status.json").write_text(json.dumps(report, indent=2) + "\n")
            except Exception:
                prior = (strict_failure(case, ligand_name, args.output) if cohort_names is None else
                         eligible_failure(args.output, cohort_names))
                if prior is None:
                    raise
        if prior is not None:
            if cohort_names is None:
                report = rescue_target(ligand_name, prior, case, args.library_dir, args.ncores,
                                       args.timeout, args.gamess_dir, args.dftb_parameter_dir,
                                       qm_path=args.output, repair_peptide_numbering=True)
            else:
                report = rescue_cohort(name, args.output, cohort_names, args.library_dir,
                                       args.ncores, args.timeout, args.workers,
                                       args.gamess_dir, args.dftb_parameter_dir)
        destination = case / "adaptive_status.json"
    else:
        report = execute_target(args.geometry_dir, args.output, args.gamess_dir,
                                args.library_dir, args.dftb_parameter_dir, args.workers,
                                args.ncores, args.timeout, dry_run=args.dry_run,
                                repair_peptide_numbering=args.repair_peptide_numbering,
                                pocket_heavy_only=args.pocket_heavy_only,
                                preserve_input_water_hydrogens=args.preserve_input_water_hydrogens)
        destination = args.output / ("plan.json" if args.dry_run else "run_status.json")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage": args.stage, "output": str(destination),
                      "status": report.get("status", "SCORED_INTERFACE_NOT_AFFINITY"),
                      "ranking_eligible": report.get("ranking_eligible", False)}))
    if args.stage == "adaptive-run" and report.get("status") == "REVIEW_REQUIRED":
        raise SystemExit("REVIEW_REQUIRED: no registered, quality-checked rescue succeeded")


def main() -> None:
    try:
        _main()
    except (ValueError, FileNotFoundError, FileExistsError, ImportError, subprocess.CalledProcessError) as error:
        report = {"schema": "fmo-tool-error-v1", "status": "REVIEW_REQUIRED", "ranking_eligible": False,
                  "error_type": type(error).__name__, "message": str(error),
                  "next_steps": ["Inspect component_model_plan.json and preparation/audit logs in the output directory.",
                                 "Retain original inputs; do not delete unsupported chemistry to force a score."]}
        if "--output" in sys.argv:
            index = sys.argv.index("--output")
            if index + 1 < len(sys.argv):
                output = Path(sys.argv[index + 1])
                if output.is_dir():
                    report["output_directory"] = str(output.resolve())
                    if "--protein" in sys.argv:
                        source_index = sys.argv.index("--protein")
                        diagnostic = output / "failure_chemistry"
                        if source_index + 1 < len(sys.argv) and not diagnostic.exists():
                            try:
                                from fmotpiefinal.core.chemical_preflight import audit_chemistry
                                diagnosis = audit_chemistry(Path(sys.argv[source_index + 1]), diagnostic)
                                report["chemistry_diagnosis"] = {"status": diagnosis["status"], "directory": str(diagnostic.resolve())}
                            except (ValueError, OSError, ImportError, KeyError) as diagnostic_error:
                                report["diagnostic_error"] = str(diagnostic_error)
                    destination = output / "tool_error.json"
                    if not destination.exists():
                        destination.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report), file=sys.stderr)
        raise SystemExit(2) from None


if __name__ == "__main__":
    main()
