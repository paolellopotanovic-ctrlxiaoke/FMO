"""Choose a supported component preparation by chemistry, never target or labels."""
import json
import shutil
from pathlib import Path
from tempfile import TemporaryDirectory

from .component_inventory import inventory, sha


def choose_route(graph):
    if graph["n_models"] != 1:
        raise ValueError("register a single source model before component preparation")
    special = [row for row in graph["residues"] if row["kind"] not in {"standard_residue", "water"}]
    if not special:
        return "NO_SPECIAL_COMPONENTS"
    if all(row["id"].split(":")[1] == "TLA" for row in special):
        from .tartrate import source_tartrate
        for row in special:
            source_tartrate(row)
        return "SOURCE_TARTRATE_MM"
    if all(row["kind"] == "ion" and len(row["atoms"]) == 1 for row in special):
        return "MONATOMIC_HALIDE_MM"
    from .counterion_preparation import supported_monatomic_components
    if supported_monatomic_components(special):
        return "MONATOMIC_ION_MM"
    if all(row["id"].split(":")[1] in {"ACE", "NMA", "NME", "UNK"} for row in special):
        from .terminal_caps import cap_identity
        for row in special:
            cap_identity(row)
        return "SOURCE_TERMINAL_CAP_MM"
    supported = all((row["kind"] == "metal" and len(row["atoms"]) == 1 and set(row["elements"]) <= {"NA", "CA"})
                    or row["id"].split(":")[1] in {"TPO", "SEP", "PTR"} for row in special)
    if supported:
        return "PHOSPHORESIDUE_ALKALI_ALKALINE_EARTH_MM"
    unsupported = [row["id"] for row in special]
    raise ValueError(f"no validated complete component route for {unsupported}; preserve raw inventory, do not delete or mutate")


def prepare_components(raw, prepared, output):
    from .special_chemistry import assemble_retained
    from .counterion_preparation import assemble_counterions
    raw, prepared, output = Path(raw), Path(prepared), Path(output)
    if output.exists():
        raise FileExistsError("component preparation requires a fresh receptor output")
    graph = inventory(raw)
    route = choose_route(graph)
    if any(output.with_suffix(suffix).exists() for suffix in (".model.json", ".dispatch.json")):
        raise FileExistsError("component provenance already exists; use a fresh output")
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".special-model-", dir=output.parent) as temporary:
        staged = Path(temporary) / output.name
        if route == "NO_SPECIAL_COMPONENTS":
            staged.write_bytes(prepared.read_bytes())
            result = {"retained": [], "forcefield_system_created": False,
                "status": "REGISTERED_STANDARD_PREPARATION_UNCHANGED", "prepared_sha256": sha(staged)}
        elif route in {"MONATOMIC_HALIDE_MM", "MONATOMIC_ION_MM"}:
            result = assemble_counterions(raw, prepared, staged)
        elif route == "SOURCE_TERMINAL_CAP_MM":
            from .terminal_caps import assemble_terminal_caps
            result = assemble_terminal_caps(raw, prepared, staged)
        elif route == "SOURCE_TARTRATE_MM":
            from .tartrate import restore_tartrate, component_molecules
            shutil.copyfile(prepared, staged)
            report = restore_tartrate(raw, staged)
            from openmm.app import PDBFile, NoCutoff
            from openmmforcefields.generators import SystemGenerator
            molecules, models = component_molecules(staged)
            pdb = PDBFile(str(staged))
            generator = SystemGenerator(forcefields=["amber14/protein.ff14SB.xml", "amber14/tip3p.xml"],
                small_molecule_forcefield="openff-2.2.0", molecules=molecules,
                nonperiodic_forcefield_kwargs={"nonbondedMethod": NoCutoff})
            system = generator.create_system(pdb.topology)
            result = {**report, "retained": models, "forcefield_system_created": True, "n_particles": system.getNumParticles()}
        else:
            result = assemble_retained(raw, prepared, staged)
        dispatch = {"schema": "chemistry-only-preparation-dispatch-v1", "route": route,
            "raw_sha256": graph["source_sha256"], "baseline_prepared_sha256": sha(prepared), "output_sha256": sha(staged),
            "experimental_labels_used": False, "ranking_eligible": False,
            "scope": "MM representation only; downstream pocket/fragment/QM validation remains required"}
        staged.with_suffix(".dispatch.json").write_text(json.dumps(dispatch, indent=2) + "\n")
        for artifact in sorted(Path(temporary).iterdir()):
            destination = output.parent / artifact.name
            if destination.exists():
                raise FileExistsError("component artifact already exists; use a fresh output")
        for artifact in sorted(Path(temporary).iterdir(), key=lambda path: (path == staged, path.name)):
            artifact.replace(output.parent / artifact.name)
    return {**result, "dispatch": dispatch}
