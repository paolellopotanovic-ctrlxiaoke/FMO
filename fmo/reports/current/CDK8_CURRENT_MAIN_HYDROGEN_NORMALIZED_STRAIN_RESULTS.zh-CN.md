# CDK8：固定重原子、仅H的QM对照

完整10例，20项新QM独立验收。原受体、水、口袋、电荷、重原子坐标与数值阈值保持不变。

| 预先冻结读数 | 原始Pearson | Spearman | Kendall |
|---|---:|---:|---:|
| matched_tpie | 0.7072 | 0.6970 | 0.5111 |
| original_total_strain_score | 0.8001 | 0.7091 | 0.5111 |
| original_tpie_plus_hydrogen_component | 0.7174 | 0.6848 | 0.4667 |
| original_tpie_plus_residual_strain | 0.7862 | 0.7455 | 0.5556 |
| hydrogen_normalized_tpie | 0.6996 | 0.7091 | 0.5556 |
| hydrogen_normalized_tpie_plus_residual_strain | 0.7831 | 0.7455 | 0.5556 |

氢松弛与剩余局部应变分布（kcal/mol）：{'hydrogen_relaxation_energy': {'min': 0.4104846949482364, 'median': 0.7015257853022612, 'max': 1.1633875046853532}, 'residual_local_strain': {'min': 3.9630574187690613, 'median': 6.067943509333647, 'max': 7.446784195657223}}。

H_component + residual等于原总应变；不是新增两个独立的物理自由能项。
孤立配体的H优化不代表结合环境H优化；tPIE重算用于检验移植这些H的接触代价，不直接修改默认制备。
剩余项不等于纯重原子构象自由能，包含氢耦合、MM/QM模型差异和局部盆地依赖。负值也保留并解释，不截断或择优。
这些都是开发集诊断，组合系数固定1，不按实验挑几何，不把新分数称为纯tPIE；完整配对不确定性见JSON。
