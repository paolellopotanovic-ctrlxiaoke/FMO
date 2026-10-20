# Private reference artifacts

This private internal repository bundles the two verified reference PDFs for reproducibility:

```text
boltz2.pdf
sophosqm_si.pdf
```

Do not redistribute these files or the bundled GAMESS engine outside the permitted internal environment.

## Verified source artifacts

| Artifact | Source | SHA-256 |
|---|---|---|
| Boltz-2 paper PDF | `https://www.biorxiv.org/content/10.1101/2025.06.14.659707v1.full.pdf` | `9893d35ceeec5953cd1a56746676362b4eb50eadb099210e832656e646b81cbf` |
| SophosQM SI PDF | ACS Omega DOI `10.1021/acsomega.2c08132` | `12e0cbf2c604466ed1599ce8a031b8dc5e1a6c534eed4d0952aec78b57ee5f13` |

The exact Boltz-2 Table 12 values and Appendix D.2.2 metric rule are recorded in
`src/sophosqm_baseline/exact_evaluation.py` and `reference_results/final/`.
The SophosQM SI coefficients used by the optional ablation are transcribed with provenance in
`scripts/ablate_sophosqm_calibration.py`.
