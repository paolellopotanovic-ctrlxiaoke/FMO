# Tankyrase 2 MFMO + PB/SA 文献式方案：本地复算与预测检验

日期：2026-09-28。**严格区分原始论文的气相 MFMO+AMBER12 PB/SA 和本项目可计算的 PCM 极化 FMO+原生 AmberTools PBSA 代理。**两者不能混称。

## 研究问题及判定标准

- 原始论文 DOI `10.1021/acsomega.8b00175`：`ΔG_score = [E_gas(AB)−E_gas(A)−E_gas(B)] + [G_PBSA(AB)−G_PBSA(A)−G_PBSA(B)] + α N_rot`；`α=0/1 kcal/mol`。其论文中的 MFMO 多层能量、RESP/AMBER 力场、ONIOM 几何和针对柔性分子的构象筛选在本项目中**没有等价实现**。原文的相关性是单靶点内部同系列分析，不能直接与四靶点样本数加权的本项目结果比较。
- 上一轮真正**无 PCM**的 FMO2 气相试算仅 **2/12** 成功，其他失败的复合物或受体有真实 SCF 不收敛。本轮对 JNK1 受体和 CDK2/28 复合物分别把 GAMESS `$FMOPRP MAXIT` 从 200 延长至 **400**，两者均再次以 `FAILED_SCF_CONVERGENCE` 结束，日志在 `runs/fmo87/tankyrase_gas_recovery/`。本轮不会把读取 PCM 日志的 `Eunco+D` 冒充成功的气相计算。它是**受 PCM 极化的溶质电子能量**。因此本轮能检验的预定义假设是：**在既有 QM 固定结构能量基础上，独立原生 PB/SA 溶剂项能否改善同一 87 例的排序**；不能检验原论文原式是否适合新靶点。

真正气相 FMO 成功的两个重合例均在 CDK2：17 的 `gas ΔE=-15.0533`、`PB/SA ΔG_solv=+39.9426`、总 `+24.8893 kcal/mol`，实验 `-7.04`；29 对应 `-14.3800`、`+65.2446`、`+50.8646`，实验 `-9.88`。在此**一对**上原式方向错误，量级亦不可信；但**只有一对**、DFTB3 代替 MFMO、几何/半径不同，绝不可报告跨靶点原论文方法指标或由此否定它的原始系列结果。
- 指标只比较同成功样本、每个靶点单独计算 Pearson 与 Kendall 后按照样本量加权，**所有参数事先固定、不从 87 个亲和力标签回归**。本队列多次用于协议探索，任何提升也只是回顾性线索，必须新的系列预先锁定协议外推。

## 可运行实现和不可忽略的边界

进一步盘点发现 `fmo87` 环境中的原生 AmberTools 2022 `pbsa -pqr` 可运行，尽管高层 `MMPBSA.py` 因缺失 `MMPBSA_mods` 不可用。因此 `scripts/run_fmo_pbsa.py` **新建** `AB`、`A`、`B` 三种状态的 PQR 文件并**独立运行三次 PBSA**，不是简单从旧 PCM 日志拣选一列。单结构参数：`inp=1, ipb=1, istrng=0, epsin=1, epsout=80, space=0.5 Å, sprob=1.4 Å, cavity_surften=0.00542 kcal/mol/Å², cavity_offset=0.92 kcal/mol`；解析原生程序输出的 `EPB` 和 `ECAVITY`，其差值相加。旋转键可选 `α=1`，同时报告无该项的事先指定主比较。

准备：在同靶点全部配体的 5 Å 残基并集上选择 PDB2PQR 处理的蛋白、固定保留的结构水；逐原子验证 PDB/PQR 原子名和坐标一致，蛋白部分电荷/半径**来自 PDB2PQR**。配体从既有最小化复合物提取相同结合姿势，使用 OpenFF NAGL GNN AM1-BCC 赋电荷、原子元素和总电荷核验；配体半径预先设定为逐元素 Bondi 风格值。PBSA 受体与复合物共有完全相同的 PQR 原子、坐标、电荷和半径；配体游离态与复合物的配体部分亦完全相同，**没有采样游离构象或构象熵**。蛋白只取原 PDB2PQR 真实原子，**不带 GAMESS 为断端添加的 ACE/NME 量子封端**；两套 Hamiltonian 原子并非严格同构，应审视截断边界敏感性。带结构水的 P38 保留该口袋水；不计算水进出。

记录：每个 PBSA 作业的输入、输出、PBSA 二进制和 SHA256、PQR/输入哈希均保留在 `runs/fmo87/tankyrase_pbsa_proxy/`；`scripts/audit_fmo_pbsa.py` 从日志重算三项差，并验证 AB PQR 恰好为 A+B 的字节级相同原子表。每个配体两次新 PBSA，每个靶点共享受体 PBSA 结果。

**主分数明确叫“PCM 极化电子能量 + PBSA”，不是“气相 FMO + PBSA”。**其电子态仍受 PCM 极化；再加 PB 反应场可能存在不受控制的溶剂效应双计，绝不等价于论文气相 MFMO 加 PBSA。`PCM_full`、`PCM_polarized_solute`、`PCM_polarized_solute+PB_polar`、`PCM_polarized_solute+PB_polar+SA`、以及最后再加固定 rotatable-bond 惩罚均作为对照，测试这项代理是否给**本队列**更好的排序。

## 复现

```bash
PY=/path/to/fmo87/bin/python
PBSA=/path/to/fmo87/bin/pbsa
$PY scripts/run_fmo_pbsa.py --pbsa-executable "$PBSA" --pilot-per-target 3
$PY scripts/run_fmo_pbsa.py --pbsa-executable "$PBSA" --all-cases --workers 8
$PY scripts/audit_fmo_pbsa.py
$PY scripts/analyze_fmo_pbsa_pairing.py
```

### 12 例预试点

`runs/fmo87/tankyrase_pbsa_proxy/pilot_3/summary.json`：12/12 新 PBSA 三状态成功。目标样本数加权的 Pearson/Kendall：原始完整 PCM `0.110/0.333`、PCM 极化电子（非气相）`0.452/0.167`、加入 PB 极性项 `0.370/0.333`、再加 PBSA 非极性项 `0.381/0.333`、最后加 `1 kcal/mol/rotor` `0.122/0.167`。只有每靶点 3 例，数据**不足以判断新方法优劣**；也不能把 PCM 极化电子项在 12 例的高指标视作已经验证的气相方法。

### 87 例结果

`runs/fmo87/tankyrase_pbsa_proxy/full/summary.json` 与 `audit.json`：**87/87** 独立三状态 PBSA 全部完成，原始 PQR 和 PBSA 日志逐例重新解析，AB = A + B 原子表逐字节核查通过。所有系数均未使用这些标签拟合。

| 同队列目标内加权分数 | Pearson | Kendall | 说明 |
|---|---:|---:|---|
| 既有 tPIE（不是本文三项法） | 0.479 | 0.320 | 当前基线经验排序信号 |
| 原完整 PCM 结合能三项差 | -0.018 | 0.043 | 独立 PBSA 未加入 |
| PCM 极化电子差 | 0.227 | 0.209 | **不是**气相 MFMO；不包含最终 PCM 溶剂项 |
| PCM 极化电子差 + 原生 PB 极性三项差 | -0.009 | 0.003 | 单独加 PB 极性项后排序下降 |
| PCM 极化电子差 + 原生 PB/SA 三项差 | 0.005 | 0.017 | 主实验：基本无整体排序信号 |
| 上一行再加 `1 kcal/mol/rotor` | -0.131 | -0.065 | 固定旋转键惩罚更差 |

**分靶点的 PB/SA 混合主分数**：CDK2 `-0.665/-0.483`（16），JNK1 `0.412/0.215`（21），P38 `-0.323/-0.127`（34），TYK2 `0.840/0.567`（16）。靶点间差异很大，TYK2 的较好表现不应掩盖 CDK2/P38 的反序。12 例试点的初始正相关**没有在完整队列保持**，不能以试点结果代替全量。

`scripts/analyze_fmo_pbsa_pairing.py` 从同 87 例旧 tPIE 原始分数和已审计混合分数作**配对**、按靶点分层 3000 次配体 bootstrap（固定随机种子；`paired_comparison.json`）。混合分数减 tPIE 的相关性差：Pearson `0.005−0.479=−0.474`，95% bootstrap 区间 `[-0.699, -0.224]`；Kendall `0.017−0.320=−0.303`，区间 `[-0.504, -0.093]`。这支持在**该固定队列**上混合代理明显不如旧 tPIE，但 bootstrap 并不能纠正反复查看测试集的模型选择偏差，不能提供新靶点泛化保证。

## 能否断言论文方法更适合实验排序？

**不能。**已经算完的最接近、但不是等价实现的 PB/SA 混合代理在现有 87 例上并未改进排序，远低于旧 tPIE；因此不能支持“只要给原有三项法加 PB/SA 就会好”的假设。但**也不能据此否定原论文方案**：必须先解决 GAS FMO 对全部目标收敛、对接/局部 QM 几何、气相密度与溶剂的严格分离、配体及蛋白的统一力场电荷/半径、截断封端、自由态/结构水，以及外部系列验证。只在已看标签的 87 例上筛选 PCM/PBSA 组合不构成新方法的独立测试。
