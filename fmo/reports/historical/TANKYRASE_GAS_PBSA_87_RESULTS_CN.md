# 87 例全量实验：气相 FMO + PBSA 与完整诊断对照

日期：2026-09-28。**结论先行：原计划的“真正气相 FMO + 独立 PBSA”在现有 GAMESS FMO2-DFTB3 固定口袋协议下，只能得到 12/87 个合格分数；因此不存在该方法的合法 87-case 排序指标。本次给出的完整 87-case 指标仅属于另行标记的“PCM 极化量子电子能 + PBSA”诊断代理，不能冒充 Tankyrase 论文的气相 MFMO+PBSA。**

## 全量执行和原始日志

- 气相：`scripts/run_fmo_total_energy.py --no-pcm` 在冻结的四靶点 87 个配体上用同靶点共享 5 Å 口袋独立运行三态 A/AB/B，无 `$PCM`，完全同一 GAMESS FMO2-DFTB3/3OB-3-1 协议。记录 `runs/fmo87/independent_gas_pbsa/gas/`：CDK2 **11/16**，JNK1 **0/21**，P38 **1/34**，TYK2 **0/16**；两靶点的共享受体 A 发生 `FAILED_SCF_CONVERGENCE`，部分另两靶点 AB 的 dimer SCF 无法收敛。程序正常结束并不等于 SCF 成功，`scripts/score_gas_fmo_pbsa.py` 逐日志解析并拒绝任何失败能量。`runs/fmo87/independent_gas_pbsa/gas_pbsa_scores/summary.json` 记录 87 项身份、12 项合格、75 项失败原因。仅 CDK2 成功的 11 例相关性为 α=0 Pearson `−0.648`、Kendall `−0.491`；**这是选择性收敛子集，不是 CDK2 全 16 或总 87 指标**。
- PBSA：`scripts/run_fmo_pbsa.py --pbsa-only --all-cases` 对同一冻结 87 例独立运行 A/B/AB 三态，**87/87** 得到有限的极性和非极性差，`scripts/audit_fmo_pbsa.py` 逐日志、二进制哈希及 PQR 原子身份校验通过：`runs/fmo87/independent_gas_pbsa/pbsa/full/audit.json`。此步骤与旧 PCM 量子结果无依赖。
- 气相收敛修复尝试：旧 12 例 5 Å 试点 2/12；调大受体/复合物 `$FMOPRP MAXIT` 至 400 的单例未收敛；本次 JNK1 受体分别以 `MCONV(4)=514`（仅阻尼）、`576`（SOSCF）、`532`（位移+DIIS）三种统一设置试跑，均失败，原始日志在 `runs/fmo87/independent_gas_pbsa/convergence_probe/`。统一 3 Å 共享口袋在 12 例试点中 CDK2 0/3、JNK1 受体失败、P38 3/3、TYK2 受体失败，见 `gas_cutoff3_pilot.stdout`，无法产生完整队列；不能按各样本结果好坏挑选口袋并拼接成一种协议。

## 独立 87/87 诊断代理：不是气相法

在**已审计的 87/87 PCM FMO 固定结构**量子日志基础上，仅取受 PCM 极化的电子项，而**不是**其最终 PCM 溶剂自由能，再与**本次重新计算并审计**的 87/87 PBSA 三态分项相加，不依赖标签拟合权重。详见 `scripts/evaluate_fmo87_hybrid_diagnostic.py`、`runs/fmo87/independent_gas_pbsa/full87_diagnostic/summary.json` 和 `per_ligand.csv`。此诊断存在 PCM 极化密度与 PB 溶剂模型不一致、可能双计响应的物理问题，仅用于量化当前路径的结果；绝不叫做“原论文方法已做 87 例”。

| 87-case、各靶点内部算指标再按该靶点样本数加权 | Pearson | Kendall |
|---|---:|---:|
| PCM 极化电子能 + 本次独立 PBSA，α=0 | **0.00522** | **0.01736** |
| 上行 + RDKit 1 kcal/mol/rotor，α=1 | **−0.13068** | **−0.06500** |
| 原 tPIE，同队列参照，**不同物理量** | 0.47923 | 0.32006 |

| 靶点 | 样本 | PCM 极化+PBSA α=0 Pearson | Kendall |
|---|---:|---:|---:|
| CDK2 | 16 | −0.66549 | −0.48333 |
| JNK1 | 21 | 0.41179 | 0.21480 |
| P38 | 34 | −0.32319 | −0.12747 |
| TYK2 | 16 | 0.84020 | 0.56667 |

这些分靶点的大幅差异意味着单一 TYK2 的好指标不能替代全队列结果；与原旧代理的全队列结果吻合，并非新验证集。当前 87 例已经多次用于探索协议：即便某个参数改良同队列指标，也不能证明外部泛化。论文式**多层 HF/MP2 MFMO**、MOE/ONIOM 几何、RESP/AMBER PBSA 参数目前都未实现；不能用 DFTB3/PCM 模拟其等价输出。

## 重算与下一道门槛

```bash
PY=/path/to/fmo87/bin/python
$PY scripts/audit_fmo_pbsa.py --runs runs/fmo87/independent_gas_pbsa/pbsa/full
$PY scripts/score_gas_fmo_pbsa.py \
  --gas-root runs/fmo87/independent_gas_pbsa/gas \
  --pbsa-root runs/fmo87/independent_gas_pbsa/pbsa/full \
  --output runs/fmo87/independent_gas_pbsa/gas_pbsa_scores --strict
$PY scripts/evaluate_fmo87_hybrid_diagnostic.py \
  --pbsa runs/fmo87/independent_gas_pbsa/pbsa/full \
  --output runs/fmo87/independent_gas_pbsa/full87_diagnostic
```

第二条 `--strict` **理应退出非零**，并留下逐例审计结果；不能在它失败后声称计算流 87/87 完成。要真正满足原计划，必须建立让**全部 87 例在同一气相 Hamiltonian 和同一统一结构协议下收敛**的解决方案；可探查非 DFTB 的原生 GAMESS HF/MP2 MFMO、统一结构水/电荷与蛋白截断策略，但不能在看了实验标签后混用不同方法拼凑 87 例或把 PCM 能量改名气相能量。87/87 后再对日志逐项重算分数并做新的系列外验证。
