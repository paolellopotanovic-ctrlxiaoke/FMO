# Galectin3：固定重原子、仅H的QM对照

完整26例，52项新QM独立验收。原受体、水、口袋、电荷、重原子坐标与数值阈值保持不变。

| 预先冻结读数 | 原始Pearson | Spearman | Kendall |
|---|---:|---:|---:|
| matched_tpie | 0.5106 | 0.4885 | 0.3375 |
| original_total_strain_score | 0.6387 | 0.6439 | 0.4675 |
| original_tpie_plus_hydrogen_component | 0.5215 | 0.5488 | 0.3994 |
| original_tpie_plus_residual_strain | 0.6722 | 0.6237 | 0.4551 |
| hydrogen_normalized_tpie | 0.3470 | 0.5614 | 0.4180 |
| hydrogen_normalized_tpie_plus_residual_strain | 0.5119 | 0.6466 | 0.4799 |

氢松弛与剩余局部应变分布（kcal/mol）：{'hydrogen_relaxation_energy': {'min': 4.017581089457782, 'median': 4.363052037838655, 'max': 10.54253441492759}, 'residual_local_strain': {'min': 10.881356272915106, 'median': 11.966248617491548, 'max': 22.492093790443825}}。

H_component + residual等于原总应变；不是新增两个独立的物理自由能项。
孤立配体的H优化不代表结合环境H优化；tPIE重算用于检验移植这些H的接触代价，不直接修改默认制备。
剩余项不等于纯重原子构象自由能，包含氢耦合、MM/QM模型差异和局部盆地依赖。负值也保留并解释，不截断或择优。
这些都是开发集诊断，组合系数固定1，不按实验挑几何，不把新分数称为纯tPIE；完整配对不确定性见JSON。
