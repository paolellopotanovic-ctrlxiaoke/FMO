# Maintained package and history

The repository root is the supported `FMOtpiefinal` project. Install it with the
root `pyproject.toml`/`environment.yml`, and use `fmotpiefinal blind-run`.

`archive/legacy_fmo87_20261020/` is the exact previously tracked tree from commit
`7b99ca7bbc47a5d66ac4d40daf99cc13f33c390d`. Nothing in the default toolkit imports
it or implicitly selects its bundled engines. The archive is preserved for
historical interpretation, not advertised as a second maintained workflow.

The archive manifest records every original tracked file's SHA256 and symlink
status. The Git rename commit preserves existing objects and history; it does
not duplicate raw binary contents as new engine downloads. Future data, logs,
environments and binaries must not be added to the maintained package.

Wheel and sdist builds exclude `archive/`, as does the clean source archive helper.
GitHub CI checks out current files without the archive and runs portable tests;
full QM checks require a configured scientific host.

Historical third-party materials remain in `archive/` and Git history. Do not
change the repository to public without a separate rights review. No new source
license grant is introduced by this integration.
