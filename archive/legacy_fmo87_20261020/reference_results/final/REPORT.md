# FMO87 / Boltz-2 exact cohort evaluation

## Audit: strict 87/87 complete

- Cohort: 87 compounds
- FMO succeeded: True
- Boltz2 succeeded: True

## Sample-weighted target-average ranking

| model | n_compounds | pearson_r | kendall_tau | pairwise_mae_kcal_mol | mae_noncentered_kcal_mol | mae_centered_kcal_mol | percent_within_1_noncentered_kcal_mol | percent_within_1_centered_kcal_mol | percent_within_2_noncentered_kcal_mol | percent_within_2_centered_kcal_mol |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| boltz2_local | 87 | 0.651763 | 0.509589 | 0.88408 | 0.900475 | 0.625489 | 0.678161 | 0.827586 | 0.931034 | 0.977011 |
| fmo_raw_tpie | 87 | 0.479233 | 0.320055 | 10.2498 | 74.2409 | 7.311 | 0 | 0.103448 | 0 | 0.149425 |
| fmo_target_fit | 87 | 0.608842 | 0.499633 | 0.932494 | 0.630937 | 0.630937 | 0.747126 | 0.747126 | 1 | 1 |
| fmo_loo | 87 | 0.401952 | 0.330867 | 1.1565 | 0.771699 | 0.777047 | 0.735632 | 0.724138 | 0.965517 | 0.965517 |

## Paper Table 12 comparison

| model | local_pearson_r | local_kendall_tau | paper_pearson_r | paper_kendall_tau | pearson_delta | kendall_delta | local_pairwise_mae_kcal_mol | paper_pairwise_mae_kcal_mol | local_mae_noncentered_kcal_mol | paper_mae_noncentered_kcal_mol | local_mae_centered_kcal_mol | paper_mae_centered_kcal_mol | local_percent_within_1_noncentered_kcal_mol | paper_percent_within_1_noncentered_kcal_mol | local_percent_within_1_centered_kcal_mol | paper_percent_within_1_centered_kcal_mol | local_percent_within_2_noncentered_kcal_mol | paper_percent_within_2_noncentered_kcal_mol | local_percent_within_2_centered_kcal_mol | paper_percent_within_2_centered_kcal_mol |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| boltz2_local | 0.651763 | 0.509589 | 0.66 | 0.48 | -0.00823679 | 0.0295891 | 0.88408 | 0.85 | 0.900475 | 0.75 | 0.625489 | 0.59 | 0.678161 | 0.69 | 0.827586 | 0.83 | 0.931034 | 0.97 | 0.977011 | 0.98 |
| fmo_raw_tpie | 0.479233 | 0.320055 | 0.55 | 0.38 | -0.0707674 | -0.059945 | 10.2498 | nan | 74.2409 | nan | 7.311 | nan | 0 | nan | 0.103448 | nan | 0 | nan | 0.149425 | nan |

## Per-target ranking

| model | target | n | pearson_r | kendall_tau | pairwise_mae_kcal_mol | mae_centered_kcal_mol |
| --- | --- | --- | --- | --- | --- | --- |
| boltz2_local | CDK2 | 16 | 0.760145 | 0.6 | 0.919124 | 0.654052 |
| boltz2_local | JNK1 | 21 | 0.796728 | 0.6253 | 0.599322 | 0.407376 |
| boltz2_local | P38 | 34 | 0.437638 | 0.360871 | 1.0336 | 0.730353 |
| boltz2_local | TYK2 | 16 | 0.808132 | 0.583333 | 0.905048 | 0.660361 |
| fmo_raw_tpie | CDK2 | 16 | 0.796505 | 0.566667 | 11.7499 | 9.11815 |
| fmo_raw_tpie | JNK1 | 21 | 0.347779 | 0.186158 | 8.3012 | 5.93505 |
| fmo_raw_tpie | P38 | 34 | 0.568153 | 0.382416 | 13.6009 | 9.46885 |
| fmo_raw_tpie | TYK2 | 16 | 0.145538 | 0.116667 | 4.18597 | 2.72438 |
| fmo_target_fit | CDK2 | 16 | 0.858937 | 0.75 | 0.726341 | 0.494824 |
| fmo_target_fit | JNK1 | 21 | 0.579875 | 0.491648 | 0.801062 | 0.539548 |
| fmo_target_fit | P38 | 34 | 0.654035 | 0.465003 | 0.866203 | 0.585617 |
| fmo_target_fit | TYK2 | 16 | 0.300732 | 0.333333 | 1.45202 | 0.983302 |
| fmo_loo | CDK2 | 16 | 0.777686 | 0.633333 | 0.89339 | 0.614344 |
| fmo_loo | JNK1 | 21 | 0.338707 | 0.34845 | 0.968071 | 0.645219 |
| fmo_loo | P38 | 34 | 0.582179 | 0.40396 | 0.935037 | 0.635455 |
| fmo_loo | TYK2 | 16 | -0.273756 | -0.15 | 2.13752 | 1.41366 |

## Stratified bootstrap 95% intervals

| model | metric | low_95 | high_95 | repetitions |
| --- | --- | --- | --- | --- |
| boltz2_local | pearson_r | 0.512878 | 0.759762 | 10000 |
| boltz2_local | kendall_tau | 0.35772 | 0.589927 | 10000 |
| fmo_raw_tpie | pearson_r | 0.311129 | 0.630096 | 10000 |
| fmo_raw_tpie | kendall_tau | 0.168225 | 0.430556 | 10000 |
| fmo_target_fit | pearson_r | 0.46551 | 0.753041 | 10000 |
| fmo_target_fit | kendall_tau | 0.352792 | 0.587397 | 10000 |
| fmo_loo | pearson_r | 0.25246 | 0.539113 | 10000 |
| fmo_loo | kendall_tau | 0.187623 | 0.437859 | 10000 |

`fmo_target_fit` and `fmo_loo` are secondary SophosQM-style calibrated models and are not compared by delta to the paper's raw in-house FMO score.

## Interpretation

- `fmo_raw_tpie` is an unfitted enthalpic score and is never reported as a fitted affinity.
- `fmo_target_fit` follows the SophosQM per-target α·tPIE + β·clogP + γ regression, but is in-sample and optimistic.
- `fmo_loo` refits each target while withholding each ligand and is the main prospective self-test.
- Paper reference: Boltz-2 0.66/0.48; FMO 0.55/0.38 (Pearson/Kendall).
- Local Boltz-2 uses seed 0, shared target MSAs, 200 affinity sampling steps, and 5 affinity diffusion samples; exact protocol differences are recorded in the final audit.
