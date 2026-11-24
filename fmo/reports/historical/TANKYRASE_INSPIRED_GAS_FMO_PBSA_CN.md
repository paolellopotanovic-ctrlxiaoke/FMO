# 自主实现论文式气相 FMO + PB/SA 预测流程

日期：2026-09-28。**可以自己实现完整的三态能量计算流程；但它不是 Okimoto 2018 的逐步骤等价复现，也不是已验证的绝对亲和力预测器。**原论文的原始 23 例数值复核与未公开输入清单见 `TANKYRASE_ORIGINAL_REPRODUCTION_CN.md`。本页描述独立四靶点队列上的**新计算**，禁止借用 PCM 量子能量充作气相量子能量。

**后续实际全量运行结果与 87 例替代诊断：`TANKYRASE_GAS_PBSA_87_RESULTS_CN.md`；当前真正无 PCM 气相流程仅 12/87，不可按本页公式报告 87-case 指标。**

## 固定计算定义

每个靶点选择一个全队列共享的受体口袋 A；每个配体提取同一结合几何的 B、AB。现有 `scripts/run_fmo_total_energy.py --no-pcm` 分别对 A、B、AB 用完整 GAMESS FMO2-DFTB3/3OB-3-1（单片段 B 用 FMO1）计算**气相**总能量；`scripts/run_fmo_pbsa.py --pbsa-only` 独立对 PQR 三态用原生 AmberTools PBSA 求极性/非极性溶剂化。 `scripts/score_gas_fmo_pbsa.py` 只从成功的原始日志读取结果，计算

```text
E_gas_bind = [E_gas(AB) - E_gas(A) - E_gas(B)] × 627.5094740631
G_polar = G_PB(AB) - G_PB(A) - G_PB(B)
G_nonpolar = G_SA(AB) - G_SA(A) - G_SA(B)
score(alpha) = E_gas_bind + G_polar + G_nonpolar + alpha × N_rot
alpha ∈ {0, 1 kcal/mol/rotatable bond}; lower score = predicted stronger.
```

它不读取实验标签来确定系数；实验 ΔG 仅用于同靶点分组的回顾性 Pearson/Kendall 比较。没有三态成功日志的配体**保持失败**，不能使用 PCM 电子能量或仅成功子集给出“全队列指标”。PI 与论文原队列的相同统计实现见 `scripts/reproduce_tankyrase_okimoto2018.py`，但不同数据集的统计结果不可相互冒称等价。

## 操作与审计

先在目标环境中设定 `$PY`、`$GAMESS_DIR`、`$GAMESS_LIB`、`$PBSA`。气相试点命令示例：

```bash
$PY scripts/run_fmo_total_energy.py --no-pcm --only-target CDK2 --limit 3 \
  --gamess-dir "$GAMESS_DIR" --library-dir "$GAMESS_LIB" \
  --dftb-parameter-dir "$GAMESS_DIR/auxdata/DFTB/3OB-3-1" \
  --output runs/fmo87/independent_gas_pbsa/gas

$PY scripts/run_fmo_pbsa.py --pbsa-only --pbsa-executable "$PBSA" \
  --pilot-per-target 3 --output runs/fmo87/independent_gas_pbsa/pbsa

$PY scripts/score_gas_fmo_pbsa.py \
  --gas-root runs/fmo87/independent_gas_pbsa/gas \
  --pbsa-root runs/fmo87/independent_gas_pbsa/pbsa/pilot_3 \
  --output runs/fmo87/independent_gas_pbsa/scores --strict
```

PBSA 从 PQR 新计算三态，不再需要任何 PCM 分项文件。评分入口逐例校验 GAS 输入 SHA-256、无 `$PCM`、正常 SCF、完整 FMO 总能量及结果缓存；PBSA 输入/PQR SHA-256、原生日志分项及复合物 PQR 的 A+B 原子表。`--strict` 在任何所选气相或 PBSA 任务失败时退出非零，仍写 `summary.json` 提供失败理由。对已成功的不足 3 个配体的靶点不计算 Pearson/Kendall。需要新预测时不得使用实验标签选择 α 或修改协议；新系列先冻结制备和评分，再对照实验。

## 试点现况及非等价项

本地**之前已运行**的 12 例无 PCM 气相试点位于 `runs/fmo87/total_sensitivity/pilot_3/gas/`；本次对其和已完成的独立 PBSA 三态原始日志审计合流，输出 `runs/fmo87/tankyrase_gas_pbsa_pilot/summary.json`：**2/12** 样本气相计算三态成功，其余 **10/12** 因 SCF 不收敛无法评分；两例同属 CDK2，**不能报告排序相关性**。另外本次对新 `--pbsa-only` 入口运行独立的四靶点每靶点 1 例 PBSA 试点，**4/4 成功**，见 `runs/fmo87/independent_gas_pbsa/pbsa/pilot_1/summary.json`；该试点并没有补齐其余气相 FMO。绝不能把成功的 2 例指标外推。

与论文差别：本项目气相部分是 DFTB3 而非原文低层 HF/STO-3G、高层 MP2/6-31G(d) 的多层 MFMO；软件手册说明 FMO-DFTB 不支持多层 MFMO。当前蛋白用 PDB2PQR + 截断封端口袋，PBSA 的 PQR 使用未封端蛋白与 OpenFF NAGL 配体电荷/简化半径，而原文使用 AMBER99SB/GAFF RESP、Tan–Luo 半径、AMBER12 Molsurf；两套 A/AB 的原子集合**并非完全相同**。配体和蛋白皆为固定结合几何，不进行论文式 MOE 姿势搜索和 Gaussian09 ONIOM 优化；缺失游离态构象系综、严格熵和标准态热力学项。**即便三态都计算成功，这只是自建论文启发式 QM/半经验排序指标，不是准确实验 ΔG 的第一性原理自由能。**

要提高物理可比性，先解决气相 SCF；随后在相同原子/封端与电荷规则下重建 QM/PBSA 三态，逐个核对极性/非极性和姿势敏感性，再将 GAMESS 切换至非 DFTB 的多层 HF/MP2 并小分子逐例对照；新系列的外部验证须与已观察的 87 例严格分离。**不要因本次仅成功 2 例便声称流已完成生产验证。**
