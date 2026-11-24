# Galectin3：保留X-H方向的键长对照

完整26例独立核验；几何、输入/日志/二进制哈希、严格收敛、完整PIEDA和一致PCM三态均复核。

| 固定方法 | 原始Pearson | Spearman | Kendall |
|---|---:|---:|---:|
| matched_tpie | 0.5106 | 0.4885 | 0.3375 |
| original_total_strain_score | 0.6387 | 0.6439 | 0.4675 |
| radial_tpie | 0.5117 | 0.4909 | 0.3437 |
| radial_tpie_plus_residual_strain | 0.6368 | 0.6312 | 0.4551 |
| original_tpie_plus_radial_component | 0.5200 | 0.4926 | 0.3437 |
| original_tpie_plus_radial_residual | 0.6402 | 0.6456 | 0.4799 |

仅几何候选，不是QM优化驻点或完整自由能；不得按target选较好方法，不自动替换生产默认。
开发集配对bootstrap不能消除重复试验选择偏差，完整逐例证据与不确定性见JSON。
