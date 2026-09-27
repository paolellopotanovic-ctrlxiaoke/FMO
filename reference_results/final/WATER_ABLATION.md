# Explicit-water ablation

All rows use the frozen 87-compound cohort and target-wise, compound-count-weighted Pearson/Kendall metrics.

| explicit_water_policy | model | pearson_r | kendall_tau | p38_pearson_r | p38_kendall_tau |
| --- | --- | --- | --- | --- | --- |
| retained_water | raw_tPIE | 0.479233 | 0.320055 | 0.568153 | 0.382416 |
| retained_water | tPIE+clogP in-sample | 0.608842 | 0.499633 | 0.654035 | 0.465003 |
| removed_from_fmo_pocket | raw_tPIE | 0.475757 | 0.315845 | 0.559260 | 0.371643 |
| removed_from_fmo_pocket | tPIE+clogP in-sample | 0.602041 | 0.488406 | 0.636632 | 0.436277 |
| removed_before_preparation | raw_tPIE | 0.413778 | 0.276553 | 0.400667 | 0.271102 |
| removed_before_preparation | tPIE+clogP in-sample | 0.557731 | 0.482793 | 0.523251 | 0.421914 |

## Paired P38 tPIE shift

| policy | mean delta vs retained water (kcal/mol) | median delta |
| --- | ---: | ---: |
| retained_water | 0.000000 | 0.000000 |
| removed_from_fmo_pocket | 20.034500 | 17.571000 |
| removed_before_preparation | 16.862647 | 17.076500 |

`removed_from_fmo_pocket` reuses the retained-water minimized complex and excludes waters only from the FMO pocket.
`removed_before_preparation` removes the 20 P38 production waters before PDB2PQR/PROPKA, then repeats OpenMM minimization and FMO.
Water PCM remains enabled in both ablations; this experiment isolates explicit crystallographic water, not implicit solvation.
