# 参考文献

## Primary methods and benchmark

- Guareschi et al. “SophosQM: Accurate Binding Affinity Prediction in Compound Optimization.” *ACS Omega* 2023, 8, 15083–15098. DOI: `10.1021/acsomega.2c08132`; PMID: `37151542`.
- Boltz-2 paper, Table 12 and Appendix D.2.2. Local cached extraction: `data/references/boltz2.txt`.
- Ross/Chen public binding free-energy benchmark structures and labels: `vendor/public_binding_free_energy_benchmark` when present.

## FMO / GAMESS

- Kitaura et al., Fragment Molecular Orbital method.
- Fedorov and Kitaura, FMO and PIEDA publications.
- GAMESS-US documentation and local registered source distribution, version 2024.2.1.
- GAMESS input manual sections for FMO, `FMOBND`/HOP, DFTB3, PIEDA, and FMO-PCM.

## Preparation and geometry

- PDB2PQR 3.6.1 and PROPKA 3.5.1 documentation.
- OpenMM 8.4 documentation.
- OpenFF 2.2.0 and OpenFF GNN AM1-BCC charge model documentation.
- RDKit Crippen logP contribution documentation.

## Interpretation notes

- SophosQM SI coefficients are transcribed in `scripts/ablate_sophosqm_calibration.py`.
- Official SI PDF and extracted text are cached under `data/references/`.
- The SI coefficients are source and series specific; they are not a universal calibration for the Boltz-2 four-target cohort.
