# Publication checklist

## Prepared here

- Standalone `src/` package, a single named CLI and descriptive module names.
- Portable tests, CI template, environment recipe and build/integrity tools.
- Source provenance with exact source/export hashes and module-renaming map.
- Aggregate development results without raw third-party molecular structures.
- No run artifacts, experiment branches, credentials or local absolute paths.

## Owner actions required before public upload

- Choose an appropriate license and confirm permission from code contributors.
  There is intentionally no fabricated LICENSE or SPDX grant.
- Set repository URL and authors/maintainers, and add an accurate `CITATION.cff`
  with the permanent citation before a tagged release.
- Review rights for aggregate benchmark tables and frozen reference statistics,
  plus GAMESS/parameter/model installation and redistribution conditions.
- Solve `environment.yml` in a clean external machine and validate the local
  GAMESS launch configuration. Current CI does not provision proprietary QM.
- Review `docs/RELEASE_CHECKS.json` and distinguish release tests from the source
  computation's 909-structure validation. Do not claim an unseen blind benchmark.

The GitHub repository preserves its previous tree under `archive/`; those legacy
files include third-party assets and are not covered by the clean-package claims.
Review their rights and existing Git history before any visibility change. The
current upload does not change repository visibility or grant a source license.

Upload **this directory**, not its parent research workspace. The source archive
is a convenience copy. Publish a public wheel only after the license review.
Building or archiving alone never uploads or changes repository visibility.
