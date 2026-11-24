"""Apply bounded graph-selected repairs before the common k10/QM pipeline."""
import json
from pathlib import Path
import shutil

from .component_inventory import sha
from .chemical_graph_policy import plan_actions
from .receptor_chemical_graph import audit
from .local_bond_repair import identify_amide_nitrogens, repair_amide_nitrogens
from .water_hydrogen_clash import hydrogen_clashes, repair_water_hydrogen_clashes
from .carbonyl_repair import selected_carbonyl_oxygens, repair_carbonyl_oxygens


def prepare_graph_receptor(source, output, source_water=None):
    source, output = Path(source), Path(output)
    output.mkdir(parents=True, exist_ok=False)
    before = audit(source)
    actions = plan_actions(before)
    (output / "before_graph.json").write_text(json.dumps(before, indent=2) + "\n")
    (output / "before_plan.json").write_text(json.dumps(actions, indent=2) + "\n")
    candidate = output / "receptor.pdb"
    selected = identify_amide_nitrogens(actions["tasks"])
    report = {"schema": "default-graph-preparation-v1", "source_sha256": sha(source),
        "experimental_labels_used": False, "selected_nitrogens": selected,
        "status": "RUNNING", "ranking_improvement_claimed": False}
    try:
        if selected and source.with_suffix(".components.json").exists():
            raise ValueError("combined local repair and special-component forcefield requires joint validation")
        if selected:
            repair_amide_nitrogens(source, actions["tasks"], candidate)
        else:
            shutil.copyfile(source, candidate)
            for suffix in (".components.json", ".components.sdf", ".model.json"):
                sidecar = source.with_suffix(suffix)
                if sidecar.exists():
                    shutil.copyfile(sidecar, candidate.with_suffix(suffix))
        carbonyl_repair = None
        carbonyl_tasks = plan_actions(audit(candidate))['tasks'] if selected else actions['tasks']
        if selected_carbonyl_oxygens(carbonyl_tasks):
            repaired_carbonyl = output / 'receptor_carbonyl_repaired.pdb'
            carbonyl_repair = repair_carbonyl_oxygens(candidate, carbonyl_tasks, repaired_carbonyl)
            shutil.copyfile(repaired_carbonyl, candidate)
        report['carbonyl_repair'] = carbonyl_repair
        water_repair = None
        if hydrogen_clashes(candidate):
            repaired_water = output / "receptor_water_H_repaired.pdb"
            water_repair = repair_water_hydrogen_clashes(candidate, source_water, repaired_water)
            shutil.copyfile(repaired_water, candidate)
            for suffix in (".components.json", ".components.sdf", ".model.json"):
                sidecar = repaired_water.with_suffix(suffix)
                if sidecar.exists():
                    shutil.copyfile(sidecar, candidate.with_suffix(suffix))
        report["water_hydrogen_repair"] = water_repair
        after = audit(candidate)
        after_plan = plan_actions(after)
        (output / "after_graph.json").write_text(json.dumps(after, indent=2) + "\n")
        (output / "after_plan.json").write_text(json.dumps(after_plan, indent=2) + "\n")
        if after_plan["quality_status"] != "PASS_CURRENT_CHEMISTRY_CHECKS":
            raise ValueError("unresolved graph quality tasks; inspect after_plan.json before k10/QM")
        report.update(status="REPAIRED_AND_VERIFIED" if selected or water_repair or carbonyl_repair else "UNCHANGED_VERIFIED",
                      output_sha256=sha(candidate), coordinates_modified=sha(candidate) != sha(source))
    except Exception as error:
        report.update(status="REVIEW_REQUIRED", failure=f"{type(error).__name__}: {error}")
        raise
    finally:
        (output / "graph_preparation.json").write_text(json.dumps(report, indent=2) + "\n")
    return candidate, report
