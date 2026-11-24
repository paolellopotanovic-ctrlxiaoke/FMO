# BACE：固定重原子、仅H的QM对照

完整36例，72项新QM独立验收。原受体、水、口袋、电荷、重原子坐标与数值阈值保持不变。

| 预先冻结读数 | 原始Pearson | Spearman | Kendall |
|---|---:|---:|---:|
| matched_tpie | 0.2319 | 0.1768 | 0.1103 |
| original_total_strain_score | 0.3847 | 0.3185 | 0.2205 |
| original_tpie_plus_hydrogen_component | 0.3373 | 0.2641 | 0.1757 |
| original_tpie_plus_residual_strain | 0.2842 | 0.1496 | 0.0863 |
| hydrogen_normalized_tpie | 0.3285 | 0.2579 | 0.1949 |
| hydrogen_normalized_tpie_plus_residual_strain | 0.3757 | 0.3478 | 0.2556 |

氢松弛与剩余局部应变分布（kcal/mol）：{'hydrogen_relaxation_energy': {'min': 8.476072298225915, 'median': 10.613836120181745, 'max': 13.509387913126501}, 'residual_local_strain': {'min': 8.10661417269339, 'median': 10.927097948518917, 'max': 16.606306555031697}}。

H_component + residual等于原总应变；不是新增两个独立的物理自由能项。
孤立配体的H优化不代表结合环境H优化；tPIE重算用于检验移植这些H的接触代价，不直接修改默认制备。
剩余项不等于纯重原子构象自由能，包含氢耦合、MM/QM模型差异和局部盆地依赖。负值也保留并解释，不截断或择优。
这些都是开发集诊断，组合系数固定1，不按实验挑几何，不把新分数称为纯tPIE；完整配对不确定性见JSON。
