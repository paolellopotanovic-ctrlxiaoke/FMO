# BACE：图触发共轭NH2修复

完整36例独立核验；几何、输入/日志/二进制哈希、严格收敛、完整PIEDA和一致PCM三态均复核。

| 固定方法 | 原始Pearson | Spearman | Kendall |
|---|---:|---:|---:|
| matched_tpie | 0.2319 | 0.1768 | 0.1103 |
| original_total_strain_score | 0.3847 | 0.3185 | 0.2205 |
| graph_repaired_tpie | 0.2867 | 0.2511 | 0.1693 |
| graph_repaired_tpie_plus_local_strain | 0.3500 | 0.2890 | 0.2045 |

仅几何候选，不是QM优化驻点或完整自由能；不得按target选较好方法，不自动替换生产默认。
开发集配对bootstrap不能消除重复试验选择偏差，完整逐例证据与不确定性见JSON。
