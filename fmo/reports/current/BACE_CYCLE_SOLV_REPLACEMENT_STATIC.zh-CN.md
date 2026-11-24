# 匹配三态溶剂修正：替换配对SOLV而不是重复相加

完整队列，保留全部配体。以下是预声明单位系数诊断组合，不称为完整tPIE或结合自由能。

| 方法 | Pearson | Spearman | Kendall |
|---|---:|---:|---:|
| baseline_full_tpie | 0.2117 | 0.1682 | 0.1118 |
| electronic_pair_plus_cycle_electrostatic | 0.2659 | 0.1885 | 0.1374 |
| electronic_pair_plus_cycle_total_solvation | 0.1353 | 0.1270 | 0.1086 |

配对SOLV与整个结合过程的去溶剂化不是同一个量。该消融先从完整tPIE扣除配对SOLV，再加入A/B/AB一致PCM计算得到的溶剂增量一次。
这仍省略受体/配体内部重组、溶液构象人口与熵；诊断组合不是严格自由能循环，不因指标较高就自动纳入默认。
参数未拟合，但数据集已经反复用于开发；需跨集合回归与独立验证。完整bootstrap及逐例分解见JSON。
