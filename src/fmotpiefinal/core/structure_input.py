"""Loss-audited receptor format and alternate-location normalization."""
from collections import defaultdict
import json
import math
from pathlib import Path
import shutil

from .component_inventory import inventory, sha


def normalize_structure(source, output, model=None):
    source, output = Path(source).resolve(), Path(output).resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    if output.exists():
        raise FileExistsError(f"use a new normalization directory: {output}")
    output.mkdir(parents=True)
    shutil.copyfile(source, output / ("original" + source.suffix.lower()))
    report = {"schema": "structure-input-v1", "source": str(source), "source_sha256": sha(source),
              "experimental_labels_used": False, "coordinates_optimized": False,
              "ranking_eligible": False, "status": "RUNNING", "chain_mapping": {}}
    try:
        if source.suffix.lower() in {".cif", ".mmcif"}:
            import gemmi
            structure = gemmi.read_structure(str(source))
            names = sorted({chain.name for structure_model in structure for chain in structure_model})
            alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
            used = {name for name in names if len(name) == 1 and name in alphabet}
            free = iter(character for character in alphabet if character not in used)
            mapping = {}
            for name in names:
                if len(name) == 1 and name in alphabet:
                    mapping[name] = name
                else:
                    replacement = next(free, None)
                    if replacement is None:
                        raise ValueError("more chains than lossless PDB downstream interface supports")
                    mapping[name] = replacement
            original_atoms = sum(structure_model.count_atom_sites() for structure_model in structure)
            if any(structure_model.count_atom_sites() > 99999 for structure_model in structure):
                raise ValueError("model exceeds downstream decimal PDB atom-serial capacity")
            for structure_model in structure:
                for chain in structure_model:
                    chain.name = mapping[chain.name]
                    for residue in chain:
                        if len(residue.name) > 3 or not -999 <= residue.seqid.num <= 9999:
                            raise ValueError("residue identity cannot be represented in downstream PDB")
                        if any(len(atom.name) > 4 for atom in residue):
                            raise ValueError("atom identity cannot be represented in downstream PDB")
            for connection in structure.connections:
                for partner in (connection.partner1, connection.partner2):
                    if partner.chain_name not in mapping:
                        raise ValueError("connection refers to an absent source chain")
                    partner.chain_name = mapping[partner.chain_name]
            text = structure.make_pdb_string()
            if len(gemmi.read_pdb_string(text).connections) != len(structure.connections):
                raise ValueError("mmCIF conversion lost explicit structural connections")
            if sum(line.startswith(("ATOM  ", "HETATM")) for line in text.splitlines()) != original_atoms:
                raise ValueError("mmCIF conversion changed atom count")
            report["chain_mapping"] = mapping
        elif source.suffix.lower() in {".pdb", ".ent"}:
            text = source.read_text()
        else:
            raise ValueError("receptor formats: .pdb, .ent, .cif, .mmcif; supply ligand chemistry separately as SDF")
        decoded = output / "decoded.pdb"
        decoded.write_text(text)
        raw_inventory = inventory(decoded)
        (output / "raw_component_inventory.json").write_text(json.dumps(raw_inventory, indent=2) + "\n")
        lines = text.splitlines(keepends=True)
        models = defaultdict(list)
        current = "1"
        for line in lines:
            if line.startswith("MODEL "):
                current = line[10:14].strip()
            if line.startswith(("ATOM  ", "HETATM")):
                models[current].append(line)
        if not models:
            raise ValueError("receptor contains no coordinate atoms")
        if model is None and len(models) != 1:
            raise ValueError(f"multiple source models {list(models)}; choose --structure-model explicitly")
        selected = str(model) if model is not None else next(iter(models))
        if selected not in models:
            raise ValueError(f"requested model {selected} not present: {list(models)}")
        groups = defaultdict(list)
        for line in models[selected]:
            groups[line[21:27]].append(line)
        kept = set()
        alternates = []
        for residue_id, atoms in groups.items():
            alternatives = sorted({line[16:17] for line in atoms} - {" "})
            chosen = None
            if alternatives:
                occupancies = {alternative: sum(float(line[54:60].strip() or "0") for line in atoms
                    if line[16:17] == alternative) / sum(line[16:17] == alternative for line in atoms)
                    for alternative in alternatives}
                if not all(math.isfinite(value) and 0 <= value <= 1 for value in occupancies.values()):
                    raise ValueError(f"invalid altloc occupancy: {residue_id}")
                chosen = min(alternatives, key=lambda alternative: (-occupancies[alternative], alternative))
                alternates.append({"residue": residue_id.strip(), "selected": chosen,
                                   "mean_occupancy": occupancies, "policy": "highest_mean_occupancy_then_lexical"})
            selected_atoms = [line for line in atoms if line[16:17] == " " or line[16:17] == chosen]
            identities = [line[12:16] for line in selected_atoms]
            if len({line[17:20] for line in selected_atoms}) != 1:
                raise ValueError(f"ambiguous residue identity: {residue_id}")
            duplicated = {name for name in identities if identities.count(name) > 1}
            permitted_tartrate_h = not alternatives and all(line[17:20] == "TLA" and
                line[76:78].strip().upper() == "H" for line in selected_atoms if line[12:16] in duplicated)
            if duplicated and (alternatives or not permitted_tartrate_h):
                raise ValueError(f"ambiguous atom/residue identity after altloc selection: {residue_id}")
            if set(identities) != {line[12:16] for line in atoms}:
                raise ValueError(f"selected altloc lacks atoms present only in another state: {residue_id}")
            kept.update(line[6:11] for line in selected_atoms)
        serials = [line[6:11] for line in models[selected]]
        if len(serials) != len(set(serials)):
            raise ValueError("duplicate source atom serials cannot safely resolve connectivity")
        candidate = output / "receptor.pdb"
        unchanged = source.suffix.lower() in {".pdb", ".ent"} and len(models) == 1 and not alternates
        if unchanged:
            shutil.copyfile(source, candidate)
        else:
            normalized = []
            current = "1"
            for line in lines:
                if line.startswith("MODEL "):
                    current = line[10:14].strip()
                    continue
                if line.startswith(("ATOM  ", "HETATM")):
                    if current == selected and line[6:11] in kept:
                        normalized.append(line[:16] + " " + line[17:])
                elif line.startswith("TER"):
                    if current == selected:
                        normalized.append(line)
                elif line.startswith("CONECT"):
                    fields = [line[start:start + 5] for start in range(6, len(line.rstrip()), 5)]
                    if fields and fields[0] in kept:
                        normalized.append("CONECT" + "".join(field for field in fields if field in kept) + "\n")
                elif not line.startswith(("MODEL ", "ENDMDL", "END", "TER", "ANISOU", "MASTER")):
                    normalized.append(line)
            candidate.write_text("".join(normalized) + "END\n")
        normalized_inventory = inventory(candidate)
        if sum(len(row["atoms"]) for row in normalized_inventory["residues"]) != len(kept):
            raise ValueError("normalization atom accounting failed")
        for suffix in (".components.json", ".components.sdf", ".model.json"):
            sidecar = source.with_suffix(suffix)
            if sidecar.exists():
                if not unchanged:
                    raise ValueError("registered chemical sidecars require revalidation after format/altloc changes")
                shutil.copyfile(sidecar, candidate.with_suffix(suffix))
        (output / "component_inventory.json").write_text(json.dumps(normalized_inventory, indent=2) + "\n")
        report.update(status="NORMALIZED", source_models=list(models), selected_model=selected,
                      selected_model_atoms=len(serials), output_atoms=len(kept), altloc_choices=alternates,
                      output_sha256=sha(candidate), byte_identical=unchanged,
                      output=str(candidate), excluded_models=[name for name in models if name != selected])
    except Exception as error:
        report.update(status="REVIEW_REQUIRED", failure=f"{type(error).__name__}: {error}")
        raise
    finally:
        (output / "normalization.json").write_text(json.dumps(report, indent=2) + "\n")
    return candidate, report


def inspect_structure(source, output, model=None):
    from .chemical_preflight import component_plan
    output = Path(output)
    if output.exists():
        raise FileExistsError("inspection requires a new output directory")
    candidate, normalization = normalize_structure(source, output / "structure_input", model)
    graph = inventory(candidate)
    plan = component_plan(graph)
    from .chemical_model_requirements import model_requirements
    requirements = model_requirements(graph)
    (output / "chemical_model_requirements.json").write_text(json.dumps(requirements, indent=2) + "\n")
    (output / "component_model_plan.json").write_text(json.dumps(plan, indent=2) + "\n")
    report = {"schema": "structure-inspection-v1", "status": "PARSED", "normalization": normalization,
        "n_residues": len(graph["residues"]), "component_plan": plan,
        "experimental_labels_used": False, "ranking_eligible": False,
        "scope": "format, identity and preparation-route inspection; no MM/QM or affinity approval"}
    (output / "inspection.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
