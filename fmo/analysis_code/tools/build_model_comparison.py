"""Audit saved FMO, Boltz2 and FEP+ scores on the same complete cohorts."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from rdkit import Chem
from scipy.stats import kendalltau, pearsonr, spearmanr


COHORTS = {
    "CDK2": (16, "jacs_set/cdk2", "jacs_set/cdk2_out.csv"),
    "JNK1": (21, "jacs_set/jnk1", "jacs_set/jnk1_manual_flips_symbmcorr_out.csv"),
    "P38": (34, "jacs_set/p38", "jacs_set/p38_out.csv"),
    "TYK2": (16, "jacs_set/tyk2", "jacs_set/tyk2_out.csv"),
    "Galectin3": (26, "misc/galectin3_extra", "misc/galectin3_extra_out.csv"),
    "FXa": (27, "opls_stress/fxa_yoshikawa_set", "opls_stress/fxa_yoshikawa_set_out.csv"),
    "CHK1": (13, "waterset/chk1", "waterset/chk1_out.csv"),
    "CDK8": (10, "misc/cdk8_koehler", "misc/cdk8_koehler_213_out.csv"),
    "BACE": (36, "jacs_set/bace", "jacs_set/bace_out.csv"),
    "HSP90": (11, "waterset/hsp90_kung", "waterset/hsp90_kung_out.csv"),
    "PTP1B": (23, "jacs_set/ptp1b", "jacs_set/ptp1b_out.csv"),
    "FAAH": (24, "misc/hfaah", "misc/hfaah_out.csv"),
    "HIV1PROT": (13, "misc/hiv_prot_ekegren", "misc/hiv_prot_ekegren_out.csv"),
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def unique_rows(rows: list[dict], field: str) -> dict:
    result = {str(row[field]).strip(): row for row in rows}
    if len(result) != len(rows):
        raise ValueError(f"duplicate identities in {field}")
    return result


def metrics(observed: list[float], predicted: list[float]) -> dict:
    if len(observed) != len(predicted) or len(observed) < 3:
        raise ValueError("matched complete observations required")
    if not all(math.isfinite(value) for value in observed + predicted):
        raise ValueError("nonfinite model score or label")
    values = {"pearson": float(pearsonr(observed, predicted).statistic),
              "spearman": float(spearmanr(observed, predicted).statistic),
              "kendall": float(kendalltau(observed, predicted).statistic)}
    if not all(math.isfinite(value) for value in values.values()):
        raise ValueError("undefined complete-cohort metric")
    return values


def check_metrics(actual: dict, expected: dict) -> None:
    for metric, value in actual.items():
        if not math.isclose(value, expected[metric], abs_tol=1e-10, rel_tol=0):
            raise ValueError(f"saved metric differs: {metric}: {value} vs {expected[metric]}")


def format_metrics(values: dict | None) -> str:
    if values is None:
        return "未找到完整结果"
    return " / ".join(f"{values[key]:.3f}".replace("-", "−") for key in ("pearson", "spearman", "kendall"))


def render_report(result: dict) -> str:
    complete = not result["missing_complete_Boltz2"]
    description = ("缺失六组117配体已补跑并完整验收，本表已自动更新。" if complete else
                   "**重要更正：不能说14个队列的Boltz2都已完成。** 还有六组117配体未获完整验收结果。")
    coverage = ("Boltz2：13队列本地完整推理，加1组OpenBind发布对照；六组补跑全部独立验收。" if complete else
                "Boltz2：7队列本地完整推理，加1组OpenBind发布对照；缺CHK1/CDK8/BACE/HSP90/PTP1B/FAAH共117配体的完整记录。")
    weak_statement = ("CHK1/CDK8/BACE：FEP+明显优于当前tPIE；新增Boltz2结果见完整表，不事后改分数或删例。" if complete else
                      "CHK1/CDK8/BACE：FEP+明显优于当前tPIE；Boltz2尚无完整核验结果，不能填数。")
    lines = ["# 全队列FMO / Boltz2 / FEP+对照", "", "核查日期：2026-10-20。FMO/FEP列为已有分数复核；新增Boltz2推理如有完成，单独记录在来源中。", "",
        description, "OpenBind是benchmark发布的490化合物Boltz2分数，不是本地重跑。", "",
        "## 完整指标表", "", "每格为 **Pearson / Spearman / Kendall**；FMO为最新2004版完整tPIE，另列软压缩Pearson。", "",
        "| 集合 | n化合物 | FMO原始 P/S/K | FMO软压缩P | Boltz2 P/S/K | FEP+21-4 P/S/K |",
        "| --- | ---: | --- | ---: | --- | --- |"]
    for target, group in result["groups"].items():
        fep = format_metrics(group["FEPplus_21_4"]) if target != "OpenBind" else "未找到匹配参考"
        lines.append(f"| {target} | {group['n']} | {format_metrics(group['FMO_raw'])} | "
                     f"{group['FMO_softcap']['pearson']:.3f} | {format_metrics(group['Boltz2'])} | {fep} |")
    weighted = result["original87_weighted"]
    lines.extend([f"| **原87按例数加权，非新增队列** | **87** | **{format_metrics(weighted['FMO_raw'])}** | "
                  f"**{weighted['FMO_softcap']['pearson']:.3f}** | **{format_metrics(weighted['Boltz2'])}** | "
                  f"**{format_metrics(weighted['FEPplus_21_4'])}** |", "",
        "## 覆盖与比较口径", "",
        "- FMO：14队列、909结构全部独立验收；OpenBind为639晶体聚合到490化合物。",
        "- " + coverage,
        "- FEP+：13个Ross队列，共270配体；逐一匹配配体身份与实验标签。OpenBind未找到对应490化合物的FEP参考，不能拿别的系列替代。",
        "- FEP+21-4是benchmark发布的校正后节点预测，不是本地新跑，也不是OpenFE。包含构象/旋转异构、对称性及pKa/互变异构修正；均值平移不影响相关。",
        "- JNK1的18628-1与18628-1_flip经原始SDF规范异构SMILES验证为同一化学身份后匹配；不是随意模糊匹配，也没有增加样本数。",
        "- 所有模型都用最新FMO评估口径的实验标签；来源标签允许最大0.001 kcal/mol的记录精度差。FXa早期Boltz指标可能在末位有微小变化。",
        "- FMO使用来源参考姿态加结构制备/k10；本地Boltz2从序列/SMILES生成自身复合物，不是同晶体姿态重打分。OpenBind直接采用其发布分数，不能套用本地推理参数说明。",
        "- FMO软压缩是单调尺度变换，S/K不变；不为Boltz2/FEP额外挑后处理。开发集合不能改称无训练重叠的前瞻盲测。", "",
        "## 直接结论", "",
        "- 原87：FEP+加权P/S/K约0.767/0.779/0.615，Boltz2约0.652/0.692/0.510，FMO原始约0.635/0.631/0.469。",
        "- Galectin3/FXa：当前FMO低于对应Boltz2和FEP+，不能把差距归为所有物理方法都差。",
        "- " + weak_statement,
        "- HIV1PROT：本系列FEP+接近无相关，Boltz2负相关，FMO也弱；这支持多方法都困难，但不能证明同一个根因。",
        "- OpenBind：当前软压缩FMO Pearson约0.416，发布Boltz2约0.413；FMO的S/K较高。仅是点估计，未在本轮重算配对显著性。", "",
        "## 证据与复现", "",
        "`MODEL_COMPARISON.json`保存精确指标、逐队列覆盖状态、名称映射和原始材料SHA256；不复制分子结构、模型权重或QM日志。",
        "`MODEL_COMPARISON.csv`便于制作表格。所有FEP原始文件来自本地benchmark的`21_4_results/ligand_predictions/`，版本说明见该目录README。", "",
        "```bash", "python tools/build_model_comparison.py --workspace /path/to/sophosqm-style-fmo-baseline", "```", "",
        "本脚本不启动QM、Boltz2或FEP。缺失项代表当前没有可报告的完整证据，不表示模型预测为零，也不保证这些运行从未存在。", ""])
    return "\n".join(lines)


def build(workspace: Path, root: Path, boltz_supplement: Path | None = None) -> dict:
    workspace = workspace.resolve()
    lab = workspace / "FMOtpiefinal_HIV1lab"
    run = lab / "runs/all_accelerated_robust_20261020"
    vendor = workspace / "vendor/public_binding_free_energy_benchmark"
    sources = {}

    def register(path: Path) -> Path:
        sources[str(path.relative_to(workspace))] = digest(path)
        return path

    def load(path: Path) -> dict:
        return json.loads(register(path).read_text())

    evaluation = load(run / "evaluation.json")
    status = load(run / "final_status.json")
    if (status["status"] != "COMPLETE_INDEPENDENTLY_VERIFIED" or
            status["n_verified"] != 909 or status["evaluation_sha256"] != digest(run / "evaluation.json")):
        raise ValueError("full FMO regression has not passed")
    register(vendor / "21_4_results/README.md")
    original = read_csv(register(workspace / "data/fmo87_cohort.csv"))
    boltz_original = read_csv(register(workspace / "runs/fmo87/final/boltz_predictions.csv"))
    hiv_boltz = load(lab / "runs/HIV1PROT_boltz2_20261018/evaluation.json")
    supplement = load(boltz_supplement.resolve()) if boltz_supplement is not None else None
    if supplement is not None and (supplement["status"] != "COMPLETE_INDEPENDENTLY_VERIFIED" or
            supplement["n_verified"] != 117 or set(supplement["groups"]) != {"CHK1", "CDK8", "BACE", "HSP90", "PTP1B", "FAAH"}):
        raise ValueError("full six-cohort supplemental Boltz2 verification required")
    result = {"schema": "complete-cohort-model-comparison-v1", "date_utc": "2026-10-20",
              "new_QM_or_Boltz_or_FEP_inference": False,
              "fep_source": "Published FEP+21-4 corrected ligand predictions, not local OpenFE or new FEP",
              "label_tolerance_kcal_mol": 0.001, "groups": {}, "sources_sha256": sources}
    if supplement is not None:
        result["new_Boltz_supplement_inference"] = True
        result["Boltz_supplement_timing"] = {
            "end_to_end_wall_seconds": supplement["timing"]["end_to_end_wall_seconds"],
            "per_target": {target: {stage: sum(attempt["wall_seconds"] for attempt in attempts)
                           for stage, attempts in timings.items()}
                           for target, timings in supplement["timing"]["timings"].items()},
            "scope": "wall time including imports/model loads; not GPU kernel or full CPU accounting"}
    for target, (expected, sdf_stem, fep_file) in COHORTS.items():
        if target in ("CDK2", "JNK1", "P38", "TYK2"):
            truth_rows = unique_rows([row for row in original if row["target"] == target], "ligand_name")
            truth = {name: -float(row["experimental_dg_kcal_mol"]) for name, row in truth_rows.items()}
        else:
            sdf = register(vendor / "fep_benchmark_inputs/structure_inputs" / f"{sdf_stem}_ligands.sdf")
            molecules = list(Chem.SDMolSupplier(str(sdf), removeHs=False))
            if any(molecule is None for molecule in molecules):
                raise ValueError(f"source molecule parse failure: {target}")
            truth = {molecule.GetProp("_Name").strip(): -float(molecule.GetProp("r_exp_dg"))
                     for molecule in molecules}
            if len(truth) != len(molecules):
                raise ValueError(f"duplicate source molecule identity: {target}")
        verification = load(run / target / "independent_verification.json")
        records = unique_rows(verification["records"], "ligand")
        fep_path = register(vendor / "21_4_results/ligand_predictions" / fep_file)
        fep = unique_rows(read_csv(fep_path), "Ligand name")
        aliases = {}
        if target == "JNK1":
            source = register(workspace / truth_rows["18628-1"]["ligand_sdf"])
            alias_molecules = {molecule.GetProp("_Name").strip(): molecule
                               for molecule in Chem.SDMolSupplier(str(source), removeHs=False)
                               if molecule is not None}
            identities = [Chem.MolToSmiles(Chem.RemoveHs(alias_molecules[name]), isomericSmiles=True)
                          for name in ("18628-1", "18628-1_flip")]
            if identities[0] != identities[1] or identities[0] != truth_rows["18628-1"]["canonical_isomeric_smiles"]:
                raise ValueError("JNK1 flipped-pose alias is not the same chemical identity")
            if "18628-1" in fep:
                raise ValueError("ambiguous JNK1 FEP pose alias")
            fep["18628-1"] = fep.pop("18628-1_flip")
            aliases = {"18628-1_flip": "18628-1"}
        if len(truth) != expected or set(truth) != set(records) or set(truth) != set(fep):
            raise ValueError(f"complete FMO/label/FEP cohort identities differ: {target}")
        if target == "Galectin3":
            label_rows = unique_rows(read_csv(register(lab / "datasets/evaluation_only/Galectin3.csv")), "ligand")
            latest_truth = {name: float(row["experimental_score_minus_dG"]) for name, row in label_rows.items()}
        elif target == "FXa":
            latest_truth = {name: -float(row["Exp. dG (kcal/mol)"]) for name, row in fep.items()}
        else:
            latest_truth = truth
        if set(latest_truth) != set(truth) or any(abs(latest_truth[name] - truth[name]) > 0.001 for name in truth):
            raise ValueError(f"latest labels differ from structural source: {target}")
        truth = latest_truth
        names = sorted(truth)
        observed = [truth[name] for name in names]
        label_difference = max(abs(truth[name] + float(fep[name]["Exp. dG (kcal/mol)"])) for name in names)
        if label_difference > 0.001:
            raise ValueError(f"FEP experimental labels differ: {target}: {label_difference}")
        raw = metrics(observed, [records[name]["raw_score"] for name in names])
        soft = metrics(observed, [records[name]["softcap_score"] for name in names])
        check_metrics(raw, evaluation["groups"][target]["raw"])
        check_metrics(soft, evaluation["groups"][target]["softcap"])
        group = {"n": expected, "complete_identity_match": True, "FMO_raw": raw, "FMO_softcap": soft,
                 "FEPplus_21_4": metrics(observed, [-float(fep[name]["Pred. dG (kcal/mol)"]) for name in names]),
                 "FEP_max_label_difference_kcal_mol": label_difference}
        group["FEP_verified_pose_aliases"] = aliases
        boltz = None
        if target in ("CDK2", "JNK1", "P38", "TYK2"):
            rows = [row for row in boltz_original if row["target"] == target]
            if any(not row["status"].startswith("SUCCEEDED") for row in rows):
                raise ValueError("incomplete original87 Boltz inference")
            boltz_rows = unique_rows(rows, "ligand_id")
            if set(boltz_rows) != set(truth):
                raise ValueError(f"Boltz identities differ: {target}")
            boltz = {name: -float(boltz_rows[name]["affinity_pred_value"]) for name in names}
            if any(abs(float(boltz_rows[name]["experimental_dg_kcal_mol"]) + truth[name]) > 1e-10 for name in names):
                raise ValueError("original87 Boltz labels differ")
            group["Boltz2_origin"] = "local complete inference; sequence/SMILES, not fixed-pose rescoring"
        elif target in ("Galectin3", "FXa"):
            path = workspace / "runs/fmo87/boltz_external_analysis" / target
            blind = load(path / "blind_predictions.json")
            unblinded = load(path / "unblinded_comparison.json")
            if blind["complete"] is not True or digest(path / "blind_predictions.json") != unblinded["blind_predictions_sha256"]:
                raise ValueError("external Boltz coverage/fingerprint mismatch")
            boltz_rows = unique_rows(unblinded["cases"], "ligand")
            if set(boltz_rows) != set(truth):
                raise ValueError("external Boltz identities differ")
            if any(abs(float(boltz_rows[name]["experimental_dg_kcal_mol"]) + truth[name]) > 0.001 for name in names):
                raise ValueError("external Boltz labels differ")
            boltz = {name: -float(boltz_rows[name]["Boltz2_affinity_pred_value"]) for name in names}
            group["Boltz2_origin"] = "local complete inference; sequence/SMILES, not fixed-pose rescoring"
        elif target == "HIV1PROT":
            boltz_rows = unique_rows(hiv_boltz["cases"], "ligand")
            if hiv_boltz["complete"] is not True or set(boltz_rows) != set(truth):
                raise ValueError("HIV Boltz identities differ")
            if any(abs(boltz_rows[name]["experimental_score_minus_dG"] - truth[name]) > 1e-10 for name in names):
                raise ValueError("HIV Boltz labels differ")
            boltz = {name: boltz_rows[name]["boltz2_rank_score_minus_affinity"] for name in names}
            group["Boltz2_origin"] = "local complete native-homodimer inference; sequence/SMILES"
        elif supplement is not None:
            completed = supplement["groups"][target]
            if completed["status"] != "COMPLETE_INDEPENDENTLY_VERIFIED" or completed["n"] != expected:
                raise ValueError("supplemental Boltz2 complete coverage not verified")
            boltz_rows = unique_rows(completed["cases"], "ligand")
            if set(boltz_rows) != set(truth):
                raise ValueError("supplemental Boltz2 identities differ")
            if any(not math.isclose(row["score"], -row["affinity_pred_value"], abs_tol=1e-12)
                   for row in boltz_rows.values()):
                raise ValueError("supplemental affinity score direction differs")
            boltz = {name: boltz_rows[name]["score"] for name in names}
            group["Boltz2_origin"] = "new complete missing-six comparator; all source protein chains, sequence/SMILES"
        group["Boltz2"] = metrics(observed, [boltz[name] for name in names]) if boltz is not None else None
        group["Boltz2_status"] = "COMPLETE_SAVED_RESULTS" if boltz is not None else "NO_COMPLETE_RESULT_LOCATED"
        result["groups"][target] = group
    paired = read_csv(register(workspace / "runs/openbind/score_distribution_audit/paired_compounds.csv"))
    if len(paired) != 490 or len(unique_rows(paired, "compound_index")) != 490:
        raise ValueError("OpenBind complete unique compound cohort required")
    published = load(workspace / "runs/openbind/score_distribution_audit/audit.json")
    openbind_boltz = metrics([float(row["experimental_pKD"]) for row in paired],
                            [float(row["boltz_2"]) for row in paired])
    check_metrics(openbind_boltz, published["metrics"]["boltz_2"])
    group = evaluation["groups"]["OpenBind"]
    latest_compounds = unique_rows(group["compound_predictions"], "compound_index")
    paired_compounds = unique_rows(paired, "compound_index")
    benchmark_path = register(workspace / "vendor/openbind/EV-A71_2A_benchmark/affinity/outputs/compound_level_prediction_analysis.csv")
    benchmark = unique_rows([row for row in read_csv(benchmark_path) if row["method"] == "boltz_2"], "smiles")
    if set(latest_compounds) != set(paired_compounds) or len(benchmark) != 490:
        raise ValueError("OpenBind published/latest compound identities differ")
    for identifier, row in paired_compounds.items():
        published_row = benchmark[row["compound_smiles"]]
        if (abs(float(row["experimental_pKD"]) - latest_compounds[identifier]["experimental_pKD"]) > 1e-10 or
                abs(float(row["experimental_pKD"]) - float(published_row["experimental_pKD"])) > 1e-10 or
                abs(float(row["boltz_2"]) - float(published_row["predicted_affinity"])) > 1e-10):
            raise ValueError("OpenBind published score/experimental label mismatch")
    result["groups"]["OpenBind"] = {"n": 490, "n_crystals": 639, "FMO_raw": group["raw"],
        "FMO_softcap": group["softcap"], "Boltz2": openbind_boltz,
        "Boltz2_status": "PUBLISHED_BENCHMARK_SCORES_NOT_LOCAL_RERUN",
        "Boltz2_origin": "benchmark boltz_2 column, matched 490 compounds; not local inference",
        "FEPplus_21_4": None, "FEP_status": "NO_MATCHED_REFERENCE_LOCATED"}
    for model in ("FMO_raw", "FMO_softcap", "Boltz2", "FEPplus_21_4"):
        result.setdefault("original87_weighted", {})[model] = {
            metric: sum(result["groups"][name]["n"] * result["groups"][name][model][metric]
                        for name in ("CDK2", "JNK1", "P38", "TYK2")) / 87
            for metric in ("pearson", "spearman", "kendall")}
    result["missing_complete_Boltz2"] = [name for name, group in result["groups"].items() if group["Boltz2"] is None]
    destination = root / "docs/MODEL_COMPARISON.json"
    destination.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    (root / "docs/MODEL_COMPARISON.zh-CN.md").write_text(render_report(result), encoding="utf-8")
    fields = ["cohort", "n"] + [f"{model}_{metric}" for model in
              ("FMO_raw", "FMO_softcap", "Boltz2", "FEPplus_21_4") for metric in ("pearson", "spearman", "kendall")]
    with (root / "docs/MODEL_COMPARISON.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for target, group in result["groups"].items():
            row = {"cohort": target, "n": group["n"]}
            for model in ("FMO_raw", "FMO_softcap", "Boltz2", "FEPplus_21_4"):
                row.update({f"{model}_{metric}": value for metric, value in (group[model] or {}).items()})
            writer.writerow(row)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", required=True, type=Path)
    parser.add_argument("--boltz-supplement", type=Path)
    arguments = parser.parse_args()
    result = build(arguments.workspace, Path(__file__).resolve().parents[1], arguments.boltz_supplement)
    for name, group in result["groups"].items():
        formatted = {model: "/".join(f"{value:.3f}" for value in group[model].values())
                     if group[model] is not None else "MISSING" for model in ("Boltz2", "FEPplus_21_4")}
        print(name, formatted)
    print("Original87 weighted:", result["original87_weighted"])
