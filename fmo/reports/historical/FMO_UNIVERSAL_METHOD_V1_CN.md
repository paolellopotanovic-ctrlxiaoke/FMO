# FMO Universal Method v1（当前证据冻结版）

日期：2026-10-04。

本文把当前已实现、已审计和仍未证明的内容分开。目标不是承诺一个对未知靶点必然高相关的单一能量公式，而是建立一套**标签无关、可复现、可迁移、失败时不伪造排序**的 FMO 计算方法。

## 1. 物理对象与分支关系

对同一个完整复合物 AB，先计算 FMO2/PIEDA 分解：

\[
tPIE = ES + EX + CT + DI + SOLV.
\]

当前四个候选分支是：

| 分支 | 含义 | 是否可与其他分支直接相加 |
|---|---|---|
| `tPIE` | 选定蛋白—配体 PIEDA pair 的直接相互作用总和 | 否；它已经包含下列子项 |
| `ES` | 静电 PIEDA 子项 | 否；是 `tPIE` 的子集 |
| `ES+SOLV` | 静电加溶剂响应的独立候选 | 否；不是新的独立自由能项 |
| `ES+CT+SOLV` | 静电、荷转移和溶剂响应候选 | 否；仍与 `tPIE` 嵌套 |

因此禁止使用 `tPIE + ES+SOLV`。如果要进行无量纲排序诊断，只能先在 target 内标准化：

\[
z_j = \frac{S_j - \operatorname{median}(S_j)}{1.4826\operatorname{MAD}(S_j)},
\qquad
z_{agg}=\frac1m\sum_j z_j.
\]

`z_agg` 是排序尺度，不是 kcal/mol，也不能被解释成热力学结合自由能。

## 2. 无标签决策树

### Gate A：身份和化学

冻结 cohort、ligand ID、canonical SMILES、形式电荷、SDF/PDB 原子顺序、蛋白链、二硫键、质子化假设和所有输入 hash。缺失、重复或原子身份不一致时直接 `FMO_FALLBACK`。

只允许自动补配体氢、修复明确的氢/水化学冲突；不自动移动严重非水重原子，不自动猜互变异构体或 microstate。系列间分子量、可旋转键和电荷差异只触发审计，不触发按系列后验归一化。

### Gate B：结构和水

先审计原始 pose，再审计实际 `minimized_complex.pdb`。记录最短重原子距离、vdW ratio、非水重原子 clash、水/氢 clash、配体 RMSD 和最大位移。

- 无严重冲突：运行统一 `k=10`。
- 有水冲突、修复或大位移：整个 target 注册 `k=1/10/100`，不按 ligand 选 restraint。
- 水删除只能 target-level 统一执行；保留 wet/dry 作为预注册诊断，不按亲和力切换。
- 独立 Boltz2 几何通过只能支持姿态方向，不能覆盖本地残余接触 warning。

### Gate C：QM 完整性

每个 target 必须有完整 case 覆盖、统一 receptor、统一 FMO fingerprint、SCF 成功、PIEDA pair 闭合和可解析的 selected pair。任何缺失或失败都不发布完整 ranking。

### Gate D：分支稳定性

所有注册 arm 都必须覆盖相同 ligand identity/order。计算 arm 间 pairwise Kendall-like 一致性和 median rank spread。稳定性只说明数值可重复，不说明与实验方向正确；Galectin3 是当前反例。

物理量聚合默认使用相同单位下的 median arm；robust-z mean 只作为无量纲诊断。禁止 per-ligand 取最低能、最好相关或最小误差 arm。

### Gate E：熵和构象采样

采样预算只由系统原子数、最大可旋转键数、最大配体重原子数和运行时校准决定。每个 state 使用固定 replica、block、seed、步数和原子身份。

只有同时满足以下条件才允许生成熵修正分支：

1. bound/free 覆盖完整且 protocol hash 一致；
2. 每个 replica 有效完成；
3. ESS、block variance 和重复 replica 通过；
4. split-half JS ≤ `0.15`；
5. 计算对象确实包含声明的 binding entropy，而不只是 torsion marginal proxy。

当前 E61 是 torsional bound/free proxy，不是完整 binding entropy。若 120 blocks 仍未通过，保留 uncertainty 和诊断，不把 `-TΔS` 加入任何生产分数。

## 3. target-specific 参数原则

参数差异必须由结构/数值观测触发，而不是由实验相关性触发：

| 结构状态 | 处理 |
|---|---|
| 原子身份、姿态和接触通过 | `k=10`，完整 FMO2/PCM，输出四分支 |
| 水/氢冲突 | target-level 统一 water policy，wet/dry 对照 |
| restraint 敏感或大位移 | 完整 `k=1/10/100`，median 物理聚合，失败则 review |
| 外部姿态方向未独立支持 | 保留所有候选分支，状态 `FMO_REVIEW` |
| SCF/PIEDA 不完整 | 不发布 ranking，状态 `FMO_FALLBACK` |
| 熵 JS/ESS 不通过 | 追加预注册 block；预算耗尽则 `ENTROPY_OMIT` |

当前生成的 target entropy plan 位于 `runs/fmo87/ablation/final_analysis/fmo_target_entropy_plans.json`；六个 target 均从同一 low-tier 起始计划开始，后续追加由收敛门决定。

## 4. 当前结果的诚实解释

- 原始 87：冻结 `ES+SOLV` 的开发集指标约 Pearson `0.634`、Kendall `0.508`；这是开发证据，不是外部普适性证明。
- Galectin3：当前 E68 26/26 QM 完成，tPIE 弱正向，PIEDA/溶剂分支反向；三种 restraint 均无严重重原子 clash，但仍有系统性 warning，因此 `FMO_REVIEW`。
- FXa：当前 E92 27/27 QM 完成，tPIE/ES/ES+SOLV 均正向；E92 是结构触发的当前 arm，E71 仅历史诊断。
- 当前三数据集的固定分支比较中，只有单独 tPIE 方向全部为正，但 Galectin3 Kendall 约 `0.174`，因此只能作为保守方向候选，不能宣称“优秀普适模型”。

完整证据分别见：

- `runs/fmo87/ablation/final_analysis/fmo_universal_decision_report.json`
- `runs/fmo87/ablation/final_analysis/fmo_cross_target_generalization_audit.json`
- `runs/fmo87/ablation/final_analysis/galectin_geometry_root_cause_audit.json`
- `runs/fmo87/ablation/final_analysis/fmo_goal_completion_audit.json`

## 5. AnewFEP 的迁移边界

AnewFEP 可迁移的是工程思想：把参数误差、结构响应、采样陷阱、数值稳定性和收敛性拆开，分别设置诊断和质量门；不能直接移植 REST2、λ-window 或 RBFE 公式到静态 FMO 端点能量。

因此本方法的“普适”定义是：

1. 输入标准化和 target-level 处理规则普适；
2. 失败模式可自动识别并 fail-closed；
3. 所有候选分支和消融可复现；
4. 任何新 target 都先经过同一结构/化学/QM/收敛门；
5. 最终准确性必须通过真正独立的数据集验证，不能由旧 87 case 的后验选择证明。

