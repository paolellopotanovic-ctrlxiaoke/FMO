#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from sophosqm_baseline.exact_evaluation import (
    PAPER_BOLTZ2,
    PAPER_FMO,
    evaluate_exact_fmo87,
    finite_json,
)


def markdown_table(frame: pd.DataFrame, *, index: bool = False) -> str:
    shown = frame.reset_index() if index else frame
    headers = [str(column) for column in shown.columns]
    rows = []
    for row in shown.itertuples(index=False, name=None):
        values = []
        for value in row:
            if isinstance(value, float):
                values.append(f"{value:.6g}")
            else:
                values.append(str(value).replace("|", "\\|"))
        rows.append(values)
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    lines.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", default="data/fmo87_cohort.csv")
    parser.add_argument("--fmo", default="runs/fmo87/final/fmo_features.csv")
    parser.add_argument("--boltz", default="runs/fmo87/final/boltz_predictions.csv")
    parser.add_argument("--boltz-input-manifest", default="runs/fmo87/boltz_inputs/manifest.json")
    parser.add_argument("--output-dir", default="runs/fmo87/final")
    parser.add_argument("--bootstrap-repetitions", type=int, default=10000)
    args = parser.parse_args()

    evaluation = evaluate_exact_fmo87(
        args.cohort, args.fmo, args.boltz,
        bootstrap_repetitions=args.bootstrap_repetitions,
        boltz_input_manifest=args.boltz_input_manifest,
    )
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    evaluation.compounds.to_csv(output / "compounds.csv", index_label="cohort_key")
    evaluation.models.to_json(output / "models.json", orient="records", indent=2)
    (output / "metrics.json").write_text(
        json.dumps(finite_json({"metrics": evaluation.metrics, "audit": evaluation.audit}), indent=2) + "\n",
        encoding="utf-8",
    )

    weighted = evaluation.models.set_index("model")["sample_weighted"].apply(pd.Series)
    weighted.insert(0, "n_compounds", evaluation.models.set_index("model")["n_compounds"])
    weighted_columns = [
        "n_compounds", "pearson_r", "kendall_tau", "pairwise_mae_kcal_mol",
        "mae_noncentered_kcal_mol", "mae_centered_kcal_mol",
        "percent_within_1_noncentered_kcal_mol", "percent_within_1_centered_kcal_mol",
        "percent_within_2_noncentered_kcal_mol", "percent_within_2_centered_kcal_mol",
    ]
    comparison = pd.DataFrame(evaluation.metrics["paper_comparison"])
    comparison.to_csv(output / "paper_comparison.csv", index=False)
    per_target_rows = []
    for model, targets in evaluation.metrics["per_target"].items():
        for target, metrics in targets.items():
            per_target_rows.append({"model": model, "target": target, **metrics})
    per_target = pd.DataFrame(per_target_rows)
    per_target.to_csv(output / "per_target_metrics.csv", index=False)
    bootstrap_rows = []
    for model, metrics in evaluation.metrics["bootstrap_95"].items():
        for metric, interval in metrics.items():
            bootstrap_rows.append({
                "model": model, "metric": metric,
                "low_95": interval["low_95"], "high_95": interval["high_95"],
                "repetitions": interval["repetitions"],
            })
    bootstrap = pd.DataFrame(bootstrap_rows)
    bootstrap.to_csv(output / "bootstrap_95.csv", index=False)
    failure_text = "strict 87/87 complete" if evaluation.audit["strict_87_complete"] else "DIAGNOSTIC INCOMPLETE RUN"
    report = [
        "# FMO87 / Boltz-2 exact cohort evaluation",
        "",
        f"## Audit: {failure_text}",
        "",
        f"- Cohort: {evaluation.audit['cohort']['n_records']} compounds",
        f"- FMO succeeded: {evaluation.audit['fmo'].get('all_required_succeeded')}",
        f"- Boltz2 succeeded: {evaluation.audit['boltz2'].get('all_required_succeeded')}",
        "",
        "## Sample-weighted target-average ranking",
        "",
        markdown_table(weighted[weighted_columns], index=True),
        "",
        "## Paper Table 12 comparison",
        "",
        markdown_table(comparison, index=False),
        "",
        "## Per-target ranking",
        "",
        markdown_table(per_target[[
            "model", "target", "n", "pearson_r", "kendall_tau",
            "pairwise_mae_kcal_mol", "mae_centered_kcal_mol",
        ]]),
        "",
        "## Stratified bootstrap 95% intervals",
        "",
        markdown_table(bootstrap),
        "",
        "`fmo_target_fit` and `fmo_loo` are secondary SophosQM-style calibrated models and are not compared by delta to the paper's raw in-house FMO score.",
        "",
        "## Interpretation",
        "",
        "- `fmo_raw_tpie` is an unfitted enthalpic score and is never reported as a fitted affinity.",
        "- `fmo_target_fit` follows the SophosQM per-target α·tPIE + β·clogP + γ regression, but is in-sample and optimistic.",
        "- `fmo_loo` refits each target while withholding each ligand and is the main prospective self-test.",
        f"- Paper reference: Boltz-2 {PAPER_BOLTZ2['pearson_r']:.2f}/{PAPER_BOLTZ2['kendall_tau']:.2f}; FMO {PAPER_FMO['pearson_r']:.2f}/{PAPER_FMO['kendall_tau']:.2f} (Pearson/Kendall).",
        "- Local Boltz-2 uses seed 0, shared target MSAs, 200 affinity sampling steps, and 5 affinity diffusion samples; exact protocol differences are recorded in the final audit.",
        "",
    ]
    (output / "REPORT.md").write_text("\n".join(report), encoding="utf-8")
    print(json.dumps({
        "output_dir": str(output),
        "strict_87_complete": evaluation.audit["strict_87_complete"],
        "all_required_runs_succeeded": evaluation.audit["all_required_runs_succeeded"],
        "compounds": str(output / "compounds.csv"),
        "metrics": str(output / "metrics.json"),
        "report": str(output / "REPORT.md"),
    }, indent=2))


if __name__ == "__main__":
    main()
