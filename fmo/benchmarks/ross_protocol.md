# Frozen protocol for the public Ross et al. 2023 comparison

## Gate 0: data freeze

- Dataset repository: `bench/data/ross2023` at commit `a65d86b782a678c75388566ce00446e45283ee63`.
- The benchmark has 13 public categories and 92 matching protein/SDF/edge systems in `fep_benchmark_inputs/structure_inputs`; this observed count takes precedence over the repository README's stale statement that 103 protein structures are included.
- The repository is authoritative. Any external count, ligand identifier, pose, protonation state, tautomer state, protein structure, or experimental value that conflicts with it is rejected.
- Inputs are used exactly as supplied: no redocking, geometry minimization of supplied poses, protonation reassignment, tautomer reselection, stereoisomer regeneration, ligand renaming, or experimental-value replacement.
- Absolute experimental binding free energies are read only from each SDF record's `r_exp_dg` property, with `r_fepplus_exp_dg` accepted as the repository's documented naming variant in `charge_annhil/egfr`. Relative experimental values in `*_edges.csv` are retained for graph definitions and published-edge checks, but pairwise metrics are recomputed from the frozen absolute values.
- Each SDF record is treated as a distinct ligand state. The current freeze contains 1,343 records, of which 1,340 are evaluable; the three non-evaluable records are `macrocycles/2Q15_lig17to21_alpha05:ligand17_int` and `merck/syk_4puz_fullmap:intermediate` and `intermediate_1`. They remain in the inventory but are excluded from prediction CSVs and metrics because the repository supplies no absolute experimental value.

## Common prediction schema

Every method writes one UTF-8 CSV per category/system/method under `results/predictions/<method>/<category>/<system>.csv` for all and only that system's evaluable ligand records, with exactly:

```text
ligand,dG_pred,sigma_pred,dG_exp,sigma_exp
```

Units are kcal/mol. Missing values are forbidden. `sigma_exp` is recorded as zero unless an auditable per-measurement uncertainty is later added to the frozen reference layer.

## Frozen execution environment

- OpenFE runtime: `openfe=1.12.0`, `openmm=8.4.0`, `rdkit=2025.09.6`, with the environment recorded in `bench/metadata/environment/openfe.conda-explicit.txt` and `openfe.pip-freeze.txt`.
- Local analysis/preparation packages are installed editable and commit-pinned in `bench/metadata/software_repos.tsv`: `cinnabar`, `Lomap`, `kartograf`, `FMOkit`, and `abmptools`.
- Re-create the runtime with `envs/install.sh`; Boltz-2 dependencies and weights are not installed.
- Every benchmark file is covered by `bench/metadata/ross2023.sha256`; the official ABINIT-MP archive is covered by `bench/metadata/abinit_mp.sha256`.

## FMO scorer

- Pocket: supplied ligand pose plus all protein residues and waters with any heavy atom within 10 Å of any supplied ligand heavy atom; bridging waters remain explicit. Do not optimize supplied coordinates.
- Fragmentation: amino-acid residues, supplied ligand(s), explicit retained waters, and covalent/coordination boundaries handled by the engine's documented FMO fragmentation rules.
- Single points: FMO-DFTB3/3OB-3-1 with GAMESS-US, currently targeting the tested FMOkit reference version 15 JUL 2024 (R2 Patch 1). The public ABINIT-MP Open Ver.2 Rev.8 archive does not contain an executable or a public DFTB3 route and is retained only as an FMO2-MP2 ceiling reference if an institutional binary/build becomes available.
- Features: IFIE/PIEDA ligand-interaction fingerprints and target-level aggregates.
- Fitting: within each target only, leave-one-ligand-out linear regression from frozen features to experimental `r_exp_dg`. No cross-target fitting, feature selection on held-out ligands, deep learning, or prospective data.
- Accuracy ceiling: FMO2-MP2 on the JACS subset when compute permits; DFTB3 remains the scalable primary method.
- The public ABINIT-MP Open Ver.2 Rev.8 archive at `bench/refs/fmo/ABINIT-MP_Ver2Rev8.tar.gz` is frozen (SHA-256 `9eaed22ef43e2d41ad5b72e550d14c57400a5871ecba41940027e7c217abc578`). It is the official 3,521,809-byte sample/documentation/tool archive; its `binary-note.txt` documents the expected HPCI binaries, but no executable is included in the archive.
- GAMESS-US source registration was approved; the official GNU/x86_64 Linux archive is `gamess-current.tar.gz`, with expected MD5 `6210613bddd3930e67b7a7eba8c24e4c`. The archive is kept outside the project tree. The verified source is unpacked under `software/gamess` (GAMESS v2024.2.1 / July 15, 2024 R2 Patch 1). The current source-tree `install.info` generator emits `setenv` assignments although this release's Makefile and build scripts source it with Bash; `scripts/normalize_gamess_install_info.py` performs an explicit local conversion to Bash `export` syntax after generation. Do not store official download credentials in this workspace.
- Build machine: x86_64 Linux, GNU Fortran/GCC 13.3.0, Open MPI 4.1.7, OpenBLAS under `/usr/local`. Build and FMO execution remain pending until DDI/main compilation and version-matched smoke tests pass.
- Prepared GAMESS inputs cover only the four initial reproduction targets (`jacs_set/tyk2`, `jacs_set/cdk2`, `jacs_set/jnk1_manual_flips`, `jacs_set/p38`): 104 ligand pockets and 104 gas-phase FMO-DFTB3 single points. `work/fmo/gamess/run_manifest.csv` is the authoritative pending-job list until a real execution overwrites it.
- A GAMESS job enters parsing only if its output contains `EXECUTION OF GAMESS TERMINATED NORMALLY`; missing markers, nonzero launcher return codes, and timeouts are failures. Parsing uses FMOkit's `gamoutparser` and ligand-interacting `LIG0001L` rows to build fixed global PIEDA/IFIE aggregate features. Target-local leave-one-ligand-out ordinary least squares uses exactly `ES_sum`, `EX_sum`, `CT_sum`, `DI_sum`, `SOL_sum`, `TOTAL_sum`, and `n_interacting_fragments`; no held-out label influences its coefficient or uncertainty estimate.
- Reproduction command once an executable exists: `envs/openfe/bin/python scripts/run_fmo_gamess.py --rungms software/gamess/rungms --targets jacs_set/tyk2`, followed by `scripts/parse_fmo_pieda.py` and `scripts/fit_fmo_loo.py`.

## OpenFE physics baseline

- Start with TYK2 (`jacs_set/tyk2`), then expand only after cycle-closure residuals and TYK2 reproduction pass.
- Use the supplied protein and every supplied ligand state without structural/protonation edits.
- Plan minimal spanning-tree RBFE networks; retain Lomap atom mappings for audit.
- Use OpenFE's current stable released protocol as installed in the project environment; record protocol UUID, software version, force fields, λ schedule, solvent, salt, equilibration, production length, and random seeds in every run directory.
- Run three independent seeds where compute permits and report the mean prediction plus uncertainty propagated from network solving and seed variation.
- Use SepTop for scaffold hopping/core changes and the documented OpenFE co-alchemical-water treatment for charge-changing perturbations.
- The 2026-09-22 planning smoke test used the frozen TYK2 ligands and protein with Kartograf, a minimal spanning network, AM1-BCC charges, and three repeats. Its audited output is `work/openfe/tyk2_plan_smoke`; it proves loading, charge generation, mapping, connected-network generation, and serialization, but is not a simulation or predictive result.

## Excluded and deferred components

- Boltz-2 inference and weight download are intentionally deferred; no Boltz-2 weights are stored in this project.
- Zenodo metadata for all 38 records matching the exact OpenFE 2024 Industry Benchmarking Project query is frozen in `bench/refs/openfe_zenodo`. The declared payload totals about 220.74 GB and is not downloaded automatically; select public Ross records only when a specific reference comparison requires them.
- Private prospective layers from the original papers (for example Recursion assays and ByteDance de novo molecules) are out of scope.

## Evaluation

- Compute all pairwise ΔΔG values from ligand-level ΔG predictions and frozen experimental ΔG values; do not reuse published edge signs.
- Report pairwise RMSE and Pearson correlation on ΔG per category/system and per method.
- Aggregate the 13 categories compound-weighted, and mark methods as zero-fit or LOO-fitted.
- Bootstrap ligands within each target for uncertainty; preserve the target/category identity of every resampled compound.
