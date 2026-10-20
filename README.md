# FMOtpiefinal

[中文说明](README.zh-CN.md)

**Current tool:** the installable project at this repository root. The previous
FMO87 engineering snapshot is retained under `archive/`, not used by this package.

An auditable, label-free structure-preparation and GAMESS FMO full-tPIE ranking
toolkit. This standalone release candidate contains the maintained computation
path, not the historical research workspace.

**Status:** packaging/reorganization release candidate, based on validated
`0.1.0.dev2026102004`. It is not a universal chemistry engine or a binding-free-energy
estimator. Unsupported chemistry produces an explicit diagnostic rather than a
fabricated score. Publication still requires the owner's license/authorship review.

## Install

The complete calculation environment needs PDBFixer, PDB2PQR/PROPKA, OpenMM,
OpenFF, force-field/model assets, and a separately licensed GAMESS installation
with the DFTB parameter library. These external programs/assets are **not bundled**.
The conda environment captures versions used locally; a fresh remote environment
solve is not yet certified. Run from this project directory:

```bash
conda env create -f environment.yml
conda activate fmotpiefinal
fmotpiefinal --help
fmotpiefinal-doctor
pytest
```

To avoid downloading historical engines/data into a new working tree, clone the
maintained files with sparse checkout (the archive remains available on GitHub):

```bash
git clone --filter=blob:none --no-checkout https://github.com/paolellopotanovic-ctrlxiaoke/FMO.git
cd FMO
git sparse-checkout set --no-cone '/*' '!/archive/'
git checkout main
```

For parsing, scoring and dependency-light development only:

```bash
python -m pip install -e '.[test]'
pytest
```

The second installation is not enough to perform full structure preparation.
See [installation and troubleshooting](docs/INSTALLATION.md).

Before a QM batch, inspect the explicit runtime without launching calculations:

```bash
fmotpiefinal-doctor --gamess-dir /path/to/gamess \
  --library-dir /path/to/runtime-libraries \
  --dftb-parameter-dir /path/to/gamess/auxdata/DFTB/3OB-3-1
```

Exit code 2 means dependencies/configuration need attention. A passing preflight
does not replace a real QM smoke test or verify every model/parameter asset.

## One complete target

Inputs are a receptor PDB/mmCIF, a label-free CSV with `ligand_id,smiles`, and an
SDF containing the corresponding receptor-frame poses. SDF names must match the
CSV IDs exactly, including stereochemistry and formal charge. This tool does not
predict a binding pose from SMILES. Never include affinity labels in these inputs.

```bash
fmotpiefinal blind-run \
  --target MY_TARGET \
  --protein data/receptor.pdb \
  --ligands data/ligands.csv \
  --poses data/poses.sdf \
  --gamess-dir /path/to/gamess \
  --library-dir /path/to/runtime-libraries \
  --dftb-parameter-dir /path/to/gamess/auxdata/DFTB/3OB-3-1 \
  --workers 4 --ncores 1 --timeout 3600 \
  --output runs/MY_TARGET
```

`--workers` controls concurrent ligand work; `--ncores` is the per-GAMESS-job
request. Account separately for preparation/OpenMM threads and native-library
threads when allocating CPU resources. Use a new output directory for a new run.

The default command includes intake, chain/connectivity normalization, PDBFixer,
PDB2PQR/PROPKA, graph/chemistry and water handling, k10 minimization, full-cohort
pocket construction, FMO2-DFTB3/PCM, convergence checks, and registered label-blind
rescue. `run` and `adaptive-run` are lower-level staged commands, not alternative
scientific scoring methods. See [workflow](docs/WORKFLOW.md) for their boundaries.

Completed runs expose `adaptive_status.json`, preparation/audit artifacts, QM
logs and `ranking/ranking.json`. Do not read a partial directory as a complete
ranking. Raw score is `-tPIE`; the optional robust soft-cap is a dimensionless
ranking proxy, **not** a predicted pKd or kcal/mol free energy.

## Repository map

```text
src/fmotpiefinal/
  pipeline/       CLI and staged orchestration
  core/           parsing, graph audits, preparation, fragmentation and GAMESS
  references/     frozen retrospective soft-cap references (no molecular assets)
  scoring.py      complete-cohort checks and robust score transformation
  diagnostics.py  label-free diagnostics
  preparation_sensitivity.py  preparation comparison utilities
tests/            portable regression tests; no research workspace imports
docs/             methods, installation, validation and publication checklist
examples/         input-schema example and parameterized launch script
tools/            release integrity and build helpers
SOURCE_PROVENANCE.json  source hashes and module-rename mapping
```

In the maintained source package, no experiment scripts, old HF/three-state/PBSA branches, raw benchmark structures,
run directories, GAMESS binaries, private machine paths, or environment copies
are bundled. The historical repository contents remain isolated under `archive/`.
Retained production rescue code is intentional, not
an invitation to select results based on affinity labels.

## Evidence and limits

The source computation version completed 909/909 structures across 14 target
cohorts, including 639 OpenBind crystals aggregated to 490 compounds. Performance
is heterogeneous: some targets remain weak or anticorrelated. The results are
retrospective development evidence, not untouched external validation. The
renamed release's tests and smoke checks are recorded separately; those must not
be presented as a second 909-structure calculation.

- [Architecture](docs/ARCHITECTURE.md)
- [Scientific workflow](docs/WORKFLOW.md)
- [Validation and target metrics](docs/VALIDATION.md)
- [GitHub publication checklist](docs/PUBLISHING.md)
- [Repository history boundary](docs/REPOSITORY_LAYOUT.md)
- [Repository integration checks](docs/REPOSITORY_CHECKS.json)

## License

No project license has been selected by the owner. Do not imply an open-source
grant or publish third-party data/software/models without permission. The archived
engineering snapshot retains its existing restrictions. Resolve the checklist
before a public release or repository visibility change.
