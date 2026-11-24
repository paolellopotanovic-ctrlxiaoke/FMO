# Validation scope

## Source computation evidence

Source version `0.1.0.dev2026102004` independently verified **909/909 structures**
across 14 target cohorts, with no hidden failed-case exclusions. Its wheel SHA256:

```text
4781636cdab3984caf63e5f6cab93bc88eca3434e6774a9d15b7ab06ec007aa9
```

The underlying complete evaluation JSON SHA256:

```text
f1689f886bf13b3f7b24cefb3f3641b82901ab30bf7df523b846fa75dad0998c
```

These are retrospective development results. Some preparations were developed
using these targets. They do not establish prospective accuracy or uniform
performance across chemistry. Not every benchmark input is an unmodified crystal
pose. OpenBind has 639 evaluated crystals aggregated to 490 compounds; its table
row uses compound-level evaluation. Other rows use ligand-level evaluation.

| Target | n | Pearson raw | Pearson soft-cap | Spearman | Kendall |
| --- | ---: | ---: | ---: | ---: | ---: |
| CDK2 | 16 | 0.868 | 0.866 | 0.856 | 0.683 |
| JNK1 | 21 | 0.614 | 0.611 | 0.628 | 0.473 |
| P38 | 34 | 0.503 | 0.500 | 0.481 | 0.329 |
| TYK2 | 16 | 0.708 | 0.761 | 0.729 | 0.550 |
| Galectin3 | 26 | 0.562 | 0.568 | 0.544 | 0.387 |
| FXa | 27 | 0.613 | 0.623 | 0.569 | 0.394 |
| CHK1 | 13 | 0.464 | 0.442 | 0.443 | 0.348 |
| CDK8 | 10 | 0.274 | 0.267 | 0.236 | 0.156 |
| BACE | 36 | 0.212 | 0.225 | 0.168 | 0.112 |
| HSP90 | 11 | 0.749 | 0.750 | 0.725 | 0.598 |
| PTP1B | 23 | 0.745 | 0.741 | 0.758 | 0.615 |
| FAAH | 24 | 0.737 | 0.732 | 0.702 | 0.541 |
| HIV1PROT | 13 | -0.147 | -0.148 | -0.047 | 0.039 |
| OpenBind | 490 | 0.371 | 0.416 | 0.441 | 0.288 |

Full-precision aggregate values are in `development_results.json`. Soft-cap is
monotone within each fixed target/reference, so it is not a way to improve
Spearman or Kendall. It can improve or worsen Pearson. No uniformly improved
prediction claim is supported, especially for HIV1PROT, BACE or CDK8.

## Standalone release checks

The release reorganizes imports and output labels, removes unused experiments,
and extracts the two live bound-geometry helpers. Export verification compares
the retained source with the declared transformations; extracted helper ASTs
must match. `SOURCE_PROVENANCE.json` records the checked mapping and hashes.

Release checks comprise portable tests against an installed wheel from outside
the source tree, CLI/module import checks, source/wheel content inspection, and
fresh complete QM smoke runs for an ordinary crystal and a graph-triggered
carbonyl-repair crystal. Their actual outcomes and limitations are recorded in
`RELEASE_CHECKS.json`. Fresh smoke cases compare tPIE, QM input and minimized pose
against the validated source calculation. They are **not** a new full 909-case run.

GAMESS/QM is not included in GitHub-hosted CI. The conda recipe has not been solved
from scratch on an independent host. Larger scientific changes must be checked
against the full datasets, not just these packaging smoke tests. The raw datasets
and logs remain in the private research workspace, not in this source package.
