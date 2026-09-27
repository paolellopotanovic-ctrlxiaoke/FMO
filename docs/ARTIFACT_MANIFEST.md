# Artifact manifest

## Git-tracked

- Production Python package: `src/sophosqm_baseline/`
- Ordered workflow CLIs: `scripts/`
- Total controller: `run_pipeline.sh`
- Frozen protocol: `configs/fmo87_production.yaml`
- Exact cohort: `data/fmo87_cohort.csv` and manifest
- 87-case structures: `data/fmo87_structures/`
- Shared Boltz MSAs: `data/boltz_msa/`
- Complete GAMESS engine and native runtime: `engines/gamess/`
- Reference papers/SI for private internal use: `data/references/`
- Frozen final tables/audits: `reference_results/final/`
- Tests and documentation: `tests/`, `docs/`, root guides

## Generated locally

- Python environment: `envs/fmo87/`
- Boltz official checkpoints: `engines/boltz-cache/`
- Raw GAMESS logs, minimized complexes, Boltz JSON: `runs/`
- Build products and caches: `dist/`, `__pycache__/`, `.pytest_cache/`

Boltz checkpoints are intentionally downloadable rather than tracked because individual official weights exceed GitHub's ordinary file-size limit.
