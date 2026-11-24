# AI-vs-physics Ross 2023 reproduction workspace

This workspace is the frozen, public-only environment for the Ross et al. 2023
binding-free-energy benchmark comparison. It currently contains the complete
public benchmark, OpenFE/OpenMM execution environment, FMO preparation and
analysis tooling, cinnabar, pinned software repositories, Gate 0 inventories,
checksums, and a planning-only TYK2 OpenFE smoke test.

## Layout

- `bench/data/ross2023`: complete public benchmark repository at commit `a65d86b...`.
- `bench/refs/fmo`: official ABINIT-MP Open Ver.2 Rev.8 sample/tool archive and extraction.
- `bench/refs/openfe_zenodo`: metadata for all 38 matching OpenFE 2024 Zenodo records.
- `bench/metadata`: Gate 0 counts, repository provenance, environment locks, and SHA256 manifests.
- `software`: OpenFE, cinnabar, Lomap, kartograf, FMOkit, abmptools, benchmarks, and IndustryBenchmarks2024.
- `envs/openfe`: installed OpenFE 1.12.0 / OpenMM 8.4 runtime.
- `software/gamess`: locally obtained GAMESS-US v2024.2.1 source tree; source archive itself is kept outside this workspace.
- `work/fmo/gamess`: 104 prepared FMO-DFTB3 Ross/JACS inputs across TYK2, CDK2, JNK1, and p38; current run manifest is pending and no energies have been computed yet.
- `work/openfe/tyk2_plan_smoke`: audited TYK2 planning output; no simulation was run.
- `protocols/protocol.md`: locked methodological choices and current execution boundaries.

## Rebuild and validate

```bash
./envs/install.sh
./envs/openfe/bin/python ./scripts/gate0_inventory.py
./envs/openfe/bin/python ./scripts/validate_freeze.py --hashes
```

The data manifest is `bench/metadata/ross2023.sha256`. Verify it with:

```bash
sha256sum -c bench/metadata/ross2023.sha256
```

## Current boundaries

- Boltz-2 inference is intentionally deferred and no weights or runtime are stored.
- The official public ABINIT-MP archive documents expected binaries but does not contain an executable; an institutional `abinitmp_smp` binary or approved build is required before production FMO execution.
- GAMESS source archive checksum and gzip/tar integrity now pass. Local GNU/OpenMPI/OpenBLAS configuration is prepared; compilation and FMO runtime validation are in progress.
- Zenodo raw OpenFE payloads total 220,739,316,280 bytes and are not blindly downloaded. Complete per-record metadata, filenames, sizes, DOIs, and checksums are frozen locally.
- The TYK2 OpenFE smoke test stops after charge generation, mapping, and network serialization. No HREX/MD sampling, quickrun, network solve, or scoreboard number has been generated.
