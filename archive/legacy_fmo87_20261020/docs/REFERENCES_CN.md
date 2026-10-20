# 参考文献

## Primary

- Guareschi et al., SophosQM / FMO affinity protocol, ACS Omega 2023; SI 保存在 `data/references/sophosqm_si.pdf`。
- Boltz-2 manuscript, Table 12 and Appendix D.2.2; PDF 保存在 `data/references/boltz2.pdf`。
- Ross/Chen public binding free-energy benchmark structures and labels; FMO87 输入已按 MIT 许可打包于 `data/fmo87_structures/`。

## Methods and engines

- GAMESS-US FMO/DFTB3/PIEDA official input manual and 3OB-3-1 parameter set.
- PDB2PQR/PROPKA preparation documentation.
- OpenMM, AMBER14, OpenFF 2.2.0 and GNN AM1-BCC charge model.
- Boltz-2 official code, CLI and model release.

## Interpretation

The Boltz-2 paper discloses its FMO reference metrics but not the private execution/scoring implementation. This project is therefore an independent official-like implementation, not a bit-for-bit reproduction.
