# Provenance

## External source repositories

- Public FEP+ benchmark: `https://github.com/schrodinger/public_binding_free_energy_benchmark.git`.
- Boltz source repository: `https://github.com/jwohlwend/boltz.git`.

Vendor checkouts are ignored in this repository and must be obtained under their own licenses. The frozen cohort and source-structure hashes are in `data/fmo87_cohort.manifest.json`.

## Frozen project data

- `data/fmo87_cohort.csv`, SHA-256 `19edfb8e2ccf3ab1029d848e4062d70cb35532d2b4fee63674f160b2a292997b`.
- `data/fmo87_cohort.manifest.json`, including ligand and protein source hashes.
- `data/references/boltz2.pdf` and extracted text for the Boltz-2 paper.
- `data/references/sophosqm_si.pdf` and extracted text for the SophosQM official SI.

The only benchmark input committed to the repository is the frozen FMO87 cohort and its source-hash manifest.

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
- Public Ross/Chen reference poses replace unavailable private Glide SP poses.
- OpenMM AMBER14/OpenFF replaces unavailable MacroModel/OPLS3e geometry.
- RDKit Crippen clogP replaces VolSurf only in secondary calibration experiments.
- GAMESS is the quantum/FMO numerical engine; this repository is an independent orchestration and benchmark implementation, not Boltz authors' private FMO code.
