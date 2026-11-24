# Artifact manifest

## Intended git-tracked project files

- `src/sophosqm_baseline/` — production Python package.
- `scripts/` — cohort construction, preparation, runners, audits, merge, evaluation, and ablation.
- `tests/` — production regression suite covering parsing, fragmentation, cohort, and evaluation.
- `configs/fmo87_production.yaml` — frozen production settings.
- `data/fmo87_cohort.csv` — exact 87-compound cohort.
- `data/fmo87_cohort.manifest.json` — source structure hashes.
- `data/references/` — cached paper/SI PDFs and extracted text used by audits.
- `docs/` — submission-ready project documentation.
- `README.md`, `pyproject.toml`, `.gitignore`.

## Generated locally, not committed by default

All content below is under ignored directories:

- `runs/fmo87/prepared_proteins/`
- `runs/fmo87/final_fmo_v4_standard_geometry/`
- `runs/fmo87/boltz_inputs/`
- `runs/fmo87/boltz_predictions_final/`
- `runs/fmo87/final/`
- `vendor/`
- `software/`

The final local report directory contains:

- `metrics.json`, `models.json`, `compounds.csv`
- `paper_comparison.csv`, `per_target_metrics.csv`, `bootstrap_95.csv`
- `production_audit.json`, `boltz_production_audit.json`
- `REPORT.md`, `PRODUCTION_AUDIT.md`, `BOLTZ_PRODUCTION_AUDIT.md`
- `CALIBRATION_ABLATION_CN.md` and calibration ablation CSV/JSON outputs

Regenerate these with `docs/REPRODUCTION_CN.md`. The only committed benchmark input is `data/fmo87_cohort.csv`.

## Suggested release archive

If raw runs must accompany a paper or supplement, archive them separately:

```text
fmo87-code/
fmo87-final-tables/
fmo87-raw-gamess/
fmo87-boltz-logs/
external-assets-manifest/
```

Keep vendor repositories, GAMESS binaries, and Boltz checkpoints out of the code repository and record their checksums separately.
