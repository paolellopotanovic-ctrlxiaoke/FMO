# BACE：保留X-H方向的键长对照

完整36例独立核验；几何、输入/日志/二进制哈希、严格收敛、完整PIEDA和一致PCM三态均复核。

| 固定方法 | 原始Pearson | Spearman | Kendall |
|---|---:|---:|---:|
| matched_tpie | 0.2319 | 0.1768 | 0.1103 |
| original_total_strain_score | 0.3847 | 0.3185 | 0.2205 |
| radial_tpie | 0.2237 | 0.1704 | 0.1086 |
| radial_tpie_plus_residual_strain | 0.3935 | 0.3276 | 0.2300 |
| original_tpie_plus_radial_component | 0.2135 | 0.1717 | 0.1086 |
| original_tpie_plus_radial_residual | 0.3987 | 0.3290 | 0.2300 |

仅几何候选，不是QM优化驻点或完整自由能；不得按target选较好方法，不自动替换生产默认。
开发集配对bootstrap不能消除重复试验选择偏差，完整逐例证据与不确定性见JSON。
