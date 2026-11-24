# CDK8：保留X-H方向的键长对照

完整10例独立核验；几何、输入/日志/二进制哈希、严格收敛、完整PIEDA和一致PCM三态均复核。

| 固定方法 | 原始Pearson | Spearman | Kendall |
|---|---:|---:|---:|
| matched_tpie | 0.7072 | 0.6970 | 0.5111 |
| original_total_strain_score | 0.8001 | 0.7091 | 0.5111 |
| radial_tpie | 0.7050 | 0.6970 | 0.5111 |
| radial_tpie_plus_residual_strain | 0.7975 | 0.7091 | 0.5111 |
| original_tpie_plus_radial_component | 0.7098 | 0.6970 | 0.5111 |
| original_tpie_plus_radial_residual | 0.7982 | 0.7091 | 0.5111 |

仅几何候选，不是QM优化驻点或完整自由能；不得按target选较好方法，不自动替换生产默认。
开发集配对bootstrap不能消除重复试验选择偏差，完整逐例证据与不确定性见JSON。
