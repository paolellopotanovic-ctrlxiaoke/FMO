# Independent GAMESS FMO87 / Boltz-2 pipeline

This is a **private, internal research repository**. It contains the exact Boltz-2 four-target FMO cohort, the complete production pipeline, and a bundled GAMESS runtime. Do not redistribute GAMESS or publisher PDFs outside the permitted internal environment.

The frozen benchmark is **87/87 neutral compounds** from CDK2, JNK1, P38, and TYK2—the same four-target cohort used for the Boltz-2 paper's FMO comparison.

## Pull And Run

```bash
./scripts/bootstrap_engines.sh
./run_pipeline.sh verify
./run_pipeline.sh check-inputs
./run_pipeline.sh download-boltz
./run_pipeline.sh all
```

No structure download or system GAMESS installation is required. The repository includes the 87-case PDB/SDF inputs, shared MSAs, GAMESS executable/source/parameters, MKL and OpenMPI libraries, frozen protocol, tests, audits, and final tables.

Geometry uses OpenMM/OpenFF directly; the separate OpenFE command-line environment is not required.

Boltz-2 Python is installed by bootstrap. Its official multi-GB checkpoints cannot be committed as ordinary GitHub blobs, so `download-boltz` fetches them from the official release URLs into the ignored `engines/boltz-cache/` directory and preserves their hashes for audit.

For parallel FMO:

```bash
./run_pipeline.sh prepare
FMO_JOBS=8 ./run_pipeline.sh fmo-parallel
```

## Documentation

- Chinese quick start: `START_HERE_CN.md`
- Pipeline map: `PIPELINE_MAP_CN.md`
- Engines and data: `docs/ENGINE_SETUP_CN.md`
- Full reproduction: `docs/REPRODUCTION_CN.md`
- Frozen protocol: `docs/EXPERIMENT_PROTOCOL_CN.md`
- Code index: `docs/CODE_REFERENCE_CN.md`
- Data contract: `docs/DATA_FLOW_CN.md`
- Final analysis: `docs/FINAL_REPORT_CN.md`
- Water ablation: `reference_results/final/WATER_ABLATION_CN.md`
- Provenance/artifacts: `docs/PROVENANCE.md`, `docs/ARTIFACT_MANIFEST.md`

## Frozen Protocol

- Protein preparation: PDB2PQR 3.6.1, AMBER, PROPKA 3.5.1, pH 7.4, hydrogen-bond optimization, OXT repair; P38 retains 20 prepared waters.
- Geometry: OpenMM 8.4, AMBER14 `protein.ff14SB`, OpenFF 2.2.0, GNN AM1-BCC ligand charges, flexible receptor/waters, and 100 kJ/mol/Å² ligand restraints.
- Fragmentation: shifted FACIO-like Cα–C′ partition, neutral ACE/NME terminal caps, internal `$FMOBND HOP_C`, and closed-shell/electron checks.
- FMO engine: GAMESS FMO2-DFTB3/3OB-3-1/PIEDA, all dimers, `MODORB=3`, water PCM `IFMO=-1`, and strict SCF/PIEDA/TIE gates.
- Boltz-2: version 2.2.1, seed 0, recycling 3, affinity sampling 200, diffusion samples 5, one shared MSA per target.

The bundled GAMESS runtime passed a real FMO smoke case (`P38_2ee`): `SUCCEEDED`, tPIE `-83.082 kcal/mol`, wall time `43.70 s`.

## Final Comparison

| Model | Pearson | Kendall | Note |
|---|---:|---:|---|
| Boltz-2 paper | 0.66 | 0.48 | Published Table 12 reference |
| Boltz-2 paper FMO | 0.55 | 0.38 | Private in-house protocol |
| Local Boltz-2 | 0.651763 | 0.509589 | Validates cohort/evaluator reconstruction |
| Our raw FMO tPIE | 0.479233 | 0.320055 | Label-free retained-water protocol |
| Our no-water raw tPIE | 0.413778 | 0.276553 | Full preparation-level ablation |

The fitted `tPIE+clogP` rows are `0.608842/0.499633` with water and `0.557731/0.482793` without explicit water. They are per-target in-sample sensitivity analyses and are not label-free comparisons. Full uncertainty and per-target tables are in `docs/FINAL_REPORT_CN.md`.

## Repository Layout

```text
src/sophosqm_baseline/       13 production modules
scripts/                     14 ordered workflow tools
configs/fmo87_production.yaml
data/fmo87_structures/       bundled 87-case inputs
data/boltz_msa/              four frozen shared MSAs
engines/boltz/               official Boltz source at the frozen commit
engines/gamess/              complete GAMESS and portable runtime
reference_results/final/     frozen 87/87 audits, metrics, reports
tests/                       production regression suite
```

Generated environments, caches, and raw logs are ignored. The final reproducible summaries are retained under `reference_results/final/`.
