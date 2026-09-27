# 显式水消融实验

以下均为冻结 87 化合物队列；指标先逐靶点计算，再按化合物数加权。

| explicit_water_policy | model | pearson_r | kendall_tau | p38_pearson_r | p38_kendall_tau |
| --- | --- | --- | --- | --- | --- |
| retained_water | raw_tPIE | 0.479233 | 0.320055 | 0.568153 | 0.382416 |
| retained_water | tPIE+clogP in-sample | 0.608842 | 0.499633 | 0.654035 | 0.465003 |
| removed_from_fmo_pocket | raw_tPIE | 0.475757 | 0.315845 | 0.559260 | 0.371643 |
| removed_from_fmo_pocket | tPIE+clogP in-sample | 0.602041 | 0.488406 | 0.636632 | 0.436277 |
| removed_before_preparation | raw_tPIE | 0.413778 | 0.276553 | 0.400667 | 0.271102 |
| removed_before_preparation | tPIE+clogP in-sample | 0.557731 | 0.482793 | 0.523251 | 0.421914 |

## P38 成对 tPIE 位移

| 策略 | 相对保留水的 mean Δ(kcal/mol) | median Δ |
| --- | ---: | ---: |
| retained_water | 0.000000 | 0.000000 |
| removed_from_fmo_pocket | 20.034500 | 17.571000 |
| removed_before_preparation | 16.862647 | 17.076500 |

`removed_from_fmo_pocket` 复用保留水体系的 OpenMM minimized complex，只在构建 FMO pocket 时排除水；用于单独检验 explicit-water scorer 贡献。
`removed_before_preparation` 在 PDB2PQR/PROPKA 前移除 P38 的 20 个生产输入水，然后重做 OpenMM restrained minimization 和 FMO2-DFTB3/PIEDA；用于检验完整结构制备效应。
两个消融均保留 water PCM。实验只移除 explicit crystallographic water，不移除 implicit solvation。
结论：去水没有改善指标。scorer 级去水略降低 Pearson/Kendall；制备前去水明显降低 P38 和全队列 raw-tPIE 排序。生产协议继续保留 SOP-style crystal water。
