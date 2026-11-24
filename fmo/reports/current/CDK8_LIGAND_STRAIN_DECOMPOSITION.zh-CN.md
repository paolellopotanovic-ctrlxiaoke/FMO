# CDK8：配体局部应变拆解

| 固定系数组合 | Pearson | Spearman | Kendall |
|---|---:|---:|---:|
| baseline_tpie | 0.2741 | 0.2364 | 0.1556 |
| matched_tpie | 0.2741 | 0.2364 | 0.1556 |
| baseline_tpie_plus_pcm_strain | 0.1836 | 0.1758 | 0.0667 |
| matched_tpie_plus_pcm_strain | 0.1836 | 0.1758 | 0.0667 |
| matched_tpie_plus_internal_response | 0.2229 | 0.1758 | 0.0667 |
| matched_tpie_plus_solvent_response | 0.2394 | 0.1152 | 0.0222 |

内部能响应含各自PCM极化，不等于真空形变能。溶剂响应为同一配体两种几何的PCM总项之差，不是整个结合过程的去溶剂化。
本分析在看到完整应变组合后设计，属于开发集诊断；不依据每个集合的指标挑选不同方法，不自动更改生产默认。
