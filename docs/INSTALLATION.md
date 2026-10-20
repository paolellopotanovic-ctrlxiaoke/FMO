# Installation and troubleshooting

## What is and is not installed

`pip install .` installs parsing/scoring dependencies and the CLI. The complete
preparation path also needs the conda/OpenFF stack in `environment.yml`.
PDBFixer is explicitly a runtime dependency of the compute environment, not an
optional step silently skipped when absent. `pdb2pqr` must be on `PATH`.
The ligand preparation uses OpenFF 2.2 and the NAGL
`openff-gnn-am1bcc-1.0.0.pt` model; verify those assets are available offline before
batch submission. Conda/model downloads can require network access.

Run `fmotpiefinal-doctor` after installation. Supply `--gamess-dir`, `--library-dir`
and `--dftb-parameter-dir` together for QM path checks. The JSON report identifies
missing imports and programs; exit status 2 means action is required. This offline
check does not start QM or download files, and does not prove that all force-field
assets, shared libraries or ligand-specific DFTB parameter pairs are valid.

GAMESS, its launch configuration and dynamic libraries must already work on the
execution host. Supply all three GAMESS/library/DFTB paths explicitly. The
validated host configuration does not establish compatibility with every MPI,
GAMESS build or scheduler. Run a small complete target before a large batch.

For tests, `.[test]` includes OpenMM and Gemmi. Tests requiring the complete
OpenFF/force-field stack skip explicitly when it is absent. CI is not a full QM
validation environment. For a build, use `python tools/build_release.py`; the
local setuptools/wheel backend must already be installed.

`python tools/archive_release.py` produces a separate, reproducible GitHub source
archive next to the project directory, excluding environments, runs and build
artifacts. It does not initialize a Git repository or contact GitHub.

## Diagnosing a stopped run

1. Read the structured CLI error/status and `ranking_eligible`, not just exit code.
2. Read component inventory, model plan, graph and preparation reports.
3. Check geometry/contact and protonation diagnostics before examining QM logs.
4. For QM failures, distinguish missing parameters from true nonconvergence.
5. Use registered recovery only; never delete atoms or choose a score using labels.

Unknown metal coordination, unsupported cofactors, covalent chemistry or ambiguous
chain/pose correspondence can legitimately stop a run. This is an explicit scope
limit, not permission to silently remove the component. Report the original input,
sanitized diagnostic and software versions when requesting a new model.

Do not run from inside a Python package directory, reuse another run's output,
or merge partial rankings. Keep hardware/thread settings stable when checking
numerical reproducibility; restrained minimization can reach nearby local states.
