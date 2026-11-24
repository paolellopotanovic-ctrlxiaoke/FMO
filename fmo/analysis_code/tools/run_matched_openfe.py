"""Freeze complete networks and execute explicit-solvent, paired OpenFE legs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import time
import traceback


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def read_edges(path, names):
    with Path(path).open() as handle:
        edges = [(row["Ligand 1"].strip(), row["Ligand 2"].strip())
                 for row in csv.DictReader(handle)]
    if not edges or any(first == second or first not in names or second not in names for first, second in edges):
        raise ValueError("network contains invalid or unknown ligand identities")
    if len({tuple(sorted(edge)) for edge in edges}) != len(edges):
        raise ValueError("duplicate network edges")
    reached = {names[0]}
    while True:
        expanded = reached | {name for edge in edges if reached.intersection(edge) for name in edge}
        if reached == expanded:
            break
        reached = expanded
    if reached != set(names):
        raise ValueError("complete connected cohort required")
    return edges


def make_protocol(profile, platform):
    from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol
    from openff.units import unit

    settings = RelativeHybridTopologyProtocol.default_settings()
    settings.forcefield_settings.small_molecule_forcefield = "openff-2.2.0"
    settings.partial_charge_settings.partial_charge_method = "nagl"
    settings.partial_charge_settings.nagl_model = "openff-gnn-am1bcc-1.0.0.pt"
    settings.engine_settings.compute_platform = platform
    settings.protocol_repeats = 1 if profile == "smoke" else 3
    if profile == "smoke":
        settings.simulation_settings.equilibration_length = 2.5 * unit.picosecond
        settings.simulation_settings.production_length = 5 * unit.picosecond
    settings.output_settings.output_indices = "all"
    settings.output_settings.positions_write_frequency = 2.5 * unit.picosecond if profile == "smoke" else 100 * unit.picosecond
    settings.output_settings.checkpoint_interval = 2.5 * unit.picosecond if profile == "smoke" else 100 * unit.picosecond
    return RelativeHybridTopologyProtocol(settings)


def freeze(workspace, output):
    root = workspace / "FMOtpiefinal_HIV1lab/runs"
    registration = root / "graph_sidechain_transfer_20261021/frozen_protocol.json"
    cohorts = json.loads(registration.read_text())["cohorts"]
    source_names = {"BACE": "jacs_set/bace", "CDK8": "misc/cdk8_koehler"}
    targets = {}
    for target, suffix in source_names.items():
        cohort = cohorts[target]
        for path, expected in cohort["input_sha256"].items():
            if digest(path) != expected:
                raise ValueError(f"changed registered input: {path}")
        poses = dict(cohort["complexes"])
        if target == "CDK8":
            candidate = root / "main_response_regression_20261021/CDK8/sidechain_response"
            decision = json.loads((candidate / "decision.json").read_text())
            if not decision["accepted"] or decision["affinity_labels_used"]:
                raise ValueError("verified, structure-only CDK8 repair required")
            poses = {name: str(candidate / "minimized" / name / "minimized_complex.pdb")
                     for name in cohort["names"]}
        edges_path = workspace / "vendor/public_binding_free_energy_benchmark/fep_benchmark_inputs/structure_inputs" / (suffix + "_edges.csv")
        edges = read_edges(edges_path, cohort["names"])
        paths = [Path(cohort["protein"]), Path(cohort["sdf"]), edges_path, *map(Path, poses.values())]
        paths += list(Path(cohort["protein"]).parent.glob(Path(cohort["protein"]).stem + ".components.*"))
        targets[target] = {"names": cohort["names"], "edges": edges,
            "protein": cohort["protein"], "sdf": cohort["sdf"], "poses": poses,
            "reference_receptor_pose": poses[sorted(poses)[0]] if target == "CDK8" else None,
            "input_sha256": {str(path): digest(path) for path in paths},
            "assay_ph": None, "assay_ph_verified": False,
            "chemical_state_policy": "conditional on current explicit source states; not a pH population prediction"}
    plan = {"schema": "matched-openfe-research-v1", "targets": targets,
        "driver_sha256": digest(__file__), "registration_sha256": digest(registration),
        "experimental_labels_used_for_generation": False, "weights_fitted": False,
        "production_main_pipeline_changed": False,
        "network_policy": "published edges, first two identity columns only; no experimental ddG used",
        "geometry_policy": "existing main-pipeline poses; one common receptor per target, CDK8 lexical-first repaired pose; receptor subsequently mobile",
        "forcefield": "AMBER ff14SB, OpenFF 2.2.0, NAGL am1bcc-1.0.0, TIP3P",
        "production_schedule": {"equilibration_ns": 1, "production_ns": 5, "repeats": 3, "lambda_windows": 11},
        "smoke_schedule": {"equilibration_ps": 2.5, "production_ps": 5, "repeats": 1, "lambda_windows": 11},
        "acceptance": "smoke cannot establish affinity accuracy; require complete paired legs, independent repeat and overlap/time convergence, then full-cohort evaluation",
        "qm_correction_status": "not implemented; no interaction-only or PCM energy substituted for an explicit-solvent Hamiltonian correction"}
    path = output / "frozen_protocol.json"
    if path.exists() and json.loads(path.read_text()) != plan:
        raise ValueError("existing frozen protocol differs; use a new run directory")
    write_json(path, plan)
    return plan


def execute(workspace, output, target, edge_index, leg, profile, platform, threads):
    import openfe
    import openmm
    from openmm import app
    from openff.toolkit import Molecule
    from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
    from openff.units import unit
    from rdkit import Chem
    from gufe.protocols import execute_DAG
    from fmotpiefinal.core.bound_geometry import load_bound_ligand
    from fmotpiefinal.core.minimization import _water_compatible_pdb
    from fmotpiefinal.core.tartrate import component_molecules

    plan = json.loads((output / "frozen_protocol.json").read_text())
    if digest(__file__) != plan["driver_sha256"]:
        raise ValueError("driver changed after freeze")
    cohort = plan["targets"][target]
    for path, expected in cohort["input_sha256"].items():
        if digest(path) != expected:
            raise ValueError(f"changed frozen input: {path}")
    first, second = cohort["edges"][edge_index]
    directory = output / profile / target / f"edge_{edge_index:03d}" / leg
    directory.mkdir(parents=True, exist_ok=True)
    if (directory / "status.json").exists():
        raise ValueError("existing execution status; use preserved checkpoint/result rather than overwriting")
    started = time.monotonic()
    write_json(directory / "status.json", {"status": "RUNNING", "pid": os.getpid(), "first": first, "second": second})
    try:
        molecules = {molecule.GetProp("_Name").strip(): molecule
                     for molecule in Chem.SDMolSupplier(cohort["sdf"], removeHs=False)}
        components = {}
        for name in (first, second):
            molecule = Chem.Mol(molecules[name])
            pose = load_bound_ligand(cohort["sdf"], name, cohort["poses"][name])
            for index, atom in enumerate(pose.atoms):
                molecule.GetConformer().SetAtomPosition(index, atom.xyz)
            offmol = Molecule.from_rdkit(molecule, allow_undefined_stereo=False, hydrogens_are_explicit=True)
            offmol.assign_partial_charges(partial_charge_method="openff-gnn-am1bcc-1.0.0.pt", toolkit_registry=NAGLToolkitWrapper())
            components[name] = openfe.SmallMoleculeComponent.from_openff(offmol, name=name)
        if Chem.GetFormalCharge(molecules[first]) != Chem.GetFormalCharge(molecules[second]):
            raise ValueError("charge-changing edge requires separately validated correction protocol")
        mapper = openfe.LomapAtomMapper(time=20, element_change=False)
        mapping = next(mapper.suggest_mappings(components[first], components[second]), None)
        if mapping is None:
            raise ValueError("no validated atom mapping; retain failed edge, do not omit ligand")
        mapping.to_json(directory / "mapping.json")
        common = {"solvent": openfe.SolventComponent(ion_concentration=0.15 * unit.molar)}
        component_report = []
        if leg == "bound":
            receptor_molecules, models = component_molecules(cohort["protein"])
            source = app.PDBFile(str(_water_compatible_pdb(cohort["protein"])))
            if cohort["reference_receptor_pose"]:
                reference = app.PDBFile(cohort["reference_receptor_pose"])
                receptor_positions = reference.positions[:source.topology.getNumAtoms()]
                if len(reference.positions) != source.topology.getNumAtoms() + molecules[sorted(cohort["names"])[0]].GetNumAtoms():
                    raise ValueError("common receptor reference atom inventory mismatch")
            else:
                receptor_positions = source.positions
            modeller = app.Modeller(source.topology, receptor_positions)
            for model, offmol in zip(models, receptor_molecules, strict=True):
                chain, residue_name, residue_number = model["id"].split(":")
                selected = [residue for residue in modeller.topology.residues()
                            if (residue.chain.id, residue.name, residue.id) == (chain, residue_name, residue_number)]
                if len(selected) != 1:
                    raise ValueError("registered component identity does not match receptor")
                modeller.delete(selected)
                component = openfe.SmallMoleculeComponent.from_openff(offmol, name=model["id"])
                common[f"component_{len(component_report)}"] = component
                component_report.append(model)
            receptor_path = directory / "common_receptor.pdb"
            with receptor_path.open("w") as handle:
                app.PDBFile.writeFile(modeller.topology, modeller.positions, handle, keepIds=True)
            common["protein"] = openfe.ProteinComponent.from_pdb_file(receptor_path)
        write_json(directory / "components.json", {"retained_registered_components": component_report,
            "source_receptor_sha256": digest(cohort["protein"]), "silent_component_deletion": False})
        protocol = make_protocol(profile, platform)
        settings = protocol.settings
        state_a = openfe.ChemicalSystem({**common, "ligand": components[first]}, name=f"{target}_{first}_{leg}")
        state_b = openfe.ChemicalSystem({**common, "ligand": components[second]}, name=f"{target}_{second}_{leg}")
        transformation = openfe.Transformation(stateA=state_a, stateB=state_b, mapping=mapping, protocol=protocol)
        transformation.to_json(directory / "transformation.json")
        write_json(directory / "execution.json", {"openfe": openfe.__version__, "openmm": openmm.__version__,
            "profile": profile, "platform": platform, "cpu_threads": threads, "full_solvent_coordinates_saved": True,
            "equilibrium_affinity_claim_allowed": False, "assay_ph_verified": False,
            "input_sha256": cohort["input_sha256"], "mapping_sha256": digest(directory / "mapping.json"),
            "transformation_sha256": digest(directory / "transformation.json")})
        shared, scratch = directory / "shared", directory / "scratch"
        shared.mkdir(exist_ok=True)
        scratch.mkdir(exist_ok=True)
        dag = transformation.create()
        result = execute_DAG(dag, shared_basedir=shared, scratch_basedir=scratch,
                             keep_shared=True, keep_scratch=True, raise_error=False)
        result.to_json(directory / "dag_result.json")
        record = {"status": "SUCCEEDED" if result.ok() else "FAILED", "profile": profile,
            "first": first, "second": second, "leg": leg, "wall_seconds": time.monotonic() - started,
            "equilibrium_affinity_claim_allowed": False}
        if result.ok():
            gathered = protocol.gather([result])
            gathered.to_json(directory / "protocol_result.json")
            estimate = gathered.get_estimate()
            uncertainty = gathered.get_uncertainty() if settings.protocol_repeats > 1 else None
            record.update(estimate_kcal_mol=float(estimate.to(unit.kilocalorie_per_mole).m),
                          repeat_std_kcal_mol=(float(uncertainty.to(unit.kilocalorie_per_mole).m)
                                               if uncertainty is not None else None))
        write_json(directory / "status.json", record)
        print(json.dumps(record), flush=True)
        return result.ok()
    except Exception as error:
        write_json(directory / "status.json", {"status": "FAILED", "first": first, "second": second,
            "leg": leg, "error": f"{type(error).__name__}: {error}", "traceback": traceback.format_exc(),
            "wall_seconds": time.monotonic() - started})
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--freeze", action="store_true")
    parser.add_argument("--target", choices=("BACE", "CDK8"))
    parser.add_argument("--edge", type=int, default=0)
    parser.add_argument("--leg", choices=("bound", "solution"))
    parser.add_argument("--profile", choices=("smoke", "production"), default="smoke")
    parser.add_argument("--platform", choices=("CPU", "CUDA"), default="CPU")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    if args.threads < 1 or args.edge < 0:
        parser.error("positive threads and nonnegative edge required")
    os.environ["OPENMM_CPU_THREADS"] = str(args.threads)
    os.environ["OMP_NUM_THREADS"] = str(args.threads)
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    if args.freeze:
        plan = freeze(args.workspace.resolve(), args.output.resolve())
        print(json.dumps({target: {"n": len(cohort["names"]), "edges": len(cohort["edges"])}
                          for target, cohort in plan["targets"].items()}))
    else:
        if not args.target or not args.leg:
            parser.error("execution requires target and leg")
        if not execute(args.workspace.resolve(), args.output.resolve(), args.target,
                       args.edge, args.leg, args.profile, args.platform, args.threads):
            raise SystemExit(1)


if __name__ == "__main__":
    main()
