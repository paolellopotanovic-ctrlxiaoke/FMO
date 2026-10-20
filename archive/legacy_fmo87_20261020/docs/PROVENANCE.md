# Provenance

## Bundled and external sources

The exact 87-case structures and shared MSAs are bundled under `data/`. GAMESS-US 2024.2.1, its executable, parameters, and native runtime are bundled under `engines/gamess/`. The official Boltz source is bundled at commit `b1ebfc46ecf57f5414e0d1a6f9027bbb122c53bc`; its Python dependencies are locked in `requirements-lock.txt`, while official weights are downloaded to an ignored cache because they exceed GitHub's file-size limit.

- Public FEP+ benchmark: `https://github.com/schrodinger/public_binding_free_energy_benchmark.git`.
- Boltz source repository: `https://github.com/jwohlwend/boltz.git`.

The original vendor checkout is no longer required. Its exact selected inputs, shared MSAs, and source hashes are bundled in this private repository.

## Frozen project data

- `data/fmo87_cohort.csv`, SHA-256 `19edfb8e2ccf3ab1029d848e4062d70cb35532d2b4fee63674f160b2a292997b`.
- `data/fmo87_cohort.manifest.json`, including ligand and protein source hashes.
- `data/references/README.md` records source URLs and SHA-256 hashes. The two reference PDFs are bundled for private internal use.

The benchmark inputs, frozen cohort, shared MSAs, reference PDFs, and complete GAMESS engine are committed for private internal use.

## Software versions used for production

| Component | Version |
|---|---:|
| GAMESS-US | 2024.2.1 |
| DFTB parameters | 3OB-3-1 |
| PDB2PQR | 3.6.1 |
| PROPKA | 3.5.1 |
| OpenMM | 8.4.0 |
| OpenFF Toolkit | 0.18.1 |
| OpenFF force field | 2.2.0 |
| RDKit | 2025.09.6 |
| Boltz | 2.2.1 |
| PyTorch | 2.10.0 |

## Scientific provenance and substitutions

- SophosQM method: Guareschi et al., *ACS Omega* 2023, DOI `10.1021/acsomega.2c08132`.
- Boltz-2 Table 12 and Appendix D.2.2 define the four-target FMO cohort, count-weighted metric aggregation, and published reference values.
- Boltz-2 engine provenance: the official repository releases the model, command-line inference, general YAML/MSA processing, tokenization, featurization, and output parser. As of 2026-09-26 (`jwohlwend/boltz` commit `b1ebfc46ecf57f5414e0d1a6f9027bbb122c53bc`), it does not release the paper-specific FEP+/four-target benchmark harness, cohort construction manifest, author prediction dump, or Table 12 metric scripts; `README.md` and `docs/evaluation.md` state that these evaluation scripts and affinity predictions are forthcoming.
- Boltz-2 internal-FMO provenance: the paper identifies the FMO row as an in-house DFTB3-based FMO code using Glide docked poses, but does not release its pose manifest, scoring post-processing, calibration equation, or coefficients. The project therefore reports raw tPIE as its main public reconstruction and fitted/LOO variants only as sensitivity analyses.
- Local Boltz-2 provenance: this repository independently reconstructs the frozen 87-compound cohort, generates protein-sequence/SMILES/MSA inputs, invokes official Boltz 2.2.1, converts affinity units, and implements Appendix D.2.2 metrics. It is therefore a public-engine reconstruction and positive control, not a bit-for-bit rerun of the authors' unavailable benchmark harness.
- Water provenance: the FMO reconstruction retains existing crystallographic waters and adds none; the production cohort contains explicit waters in 34/87 cases (all P38; 403 total). Local Boltz-2 YAML inputs contain no explicit water. The Boltz-2 paper does not disclose whether its private FMO baseline retained waters.
- Water-ablation provenance: `runs/fmo87/no_water_scorer` reuses retained-water minimized complexes and excludes waters only from FMO pockets; `runs/fmo87/no_water_preparation` removes water before PDB2PQR/PROPKA, then repeats OpenMM minimization and FMO. Both branches pass strict 87/87 audits, retain water PCM, and are summarized in `reference_results/final/WATER_ABLATION_CN.md`.
- Public Ross/Chen reference poses replace unavailable private Glide SP poses.
- OpenMM AMBER14/OpenFF replaces unavailable MacroModel/OPLS3e geometry.
- RDKit Crippen clogP replaces VolSurf only in secondary calibration experiments.
- GAMESS is the quantum/FMO numerical engine; this repository is an independent orchestration and benchmark implementation, not Boltz authors' private FMO code.
