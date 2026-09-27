# SophosQM 校准消融

## 排名与误差

| model | pearson_r | kendall_tau | pairwise_mae_kcal_mol | mae_noncentered_kcal_mol | mae_centered_kcal_mol |
| --- | --- | --- | --- | --- | --- |
| raw_tpie | 0.479233 | 0.320055 | 10.2498 | 74.2409 | 7.311 |
| tPIE_only_in_sample | 0.479233 | 0.320055 | 1.02584 | 0.701761 | 0.701761 |
| clogP_only_in_sample | 0.481288 | 0.380918 | 1.03127 | 0.707227 | 0.707227 |
| tPIE_clogP_in_sample | 0.608842 | 0.499633 | 0.932494 | 0.630937 | 0.630937 |
| tPIE_only_LOO | 0.265808 | 0.176231 | 1.134 | 0.77536 | 0.775866 |
| clogP_only_LOO | 0.297842 | 0.253177 | 1.18087 | 0.8011 | 0.803654 |
| tPIE_clogP_LOO | 0.401952 | 0.330867 | 1.1565 | 0.771699 | 0.777047 |

## p38 官方系数迁移敏感性

| model | pearson_r |
| --- | --- |
| raw_tpie | 0.568153 |
| local_tPIE_clogP_in_sample | 0.654035 |
| local_tPIE_clogP_LOO | 0.582179 |
| official_p38_alpha_beta_gamma | 0.589925 |

## 解释

- SophosQM 原文是每个 series 用实验值拟合 `α/β/γ` 后，报告同一批化合物的 fitted prediction 与实验值相关性/误差；未声明 LOO、cross-validation 或外部测试集。
- `tPIE-only` 的正向仿射校准只改变绝对 kcal/mol 误差，不改变靶点内 Pearson/Kendall 排名。
- `tPIE+clogP` in-sample 是最接近 SophosQM 原文报告方式的口径，但提升不能当作外推能力；LOO 是本项目附加的严格检查，用于测试拟合系数是否稳定。
- 官方 SI 系数来自 SophosQM 六个原始 series；除 p38 目标名重合外，Boltz2 cohort 的配体、pose、几何和 QM 方法均不同，因此只做迁移敏感性。
- 本地 clogP 是 RDKit Crippen，原文是 VolSurf；本地 FMO 是 DFTB3，原文 SophosQM 主模型是 MP2/6-31G*。
- Boltz-2 论文的 in-house DFTB3 FMO baseline 未公开是否使用 raw tPIE 或 `tPIE+clogP` 校准；本表用于本地敏感性分析，不冒充其私有后处理。
