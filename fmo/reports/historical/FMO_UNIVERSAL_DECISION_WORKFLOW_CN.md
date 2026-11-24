# FMO 通用决策流水线（当前冻结版）

更新时间：2026-10-04。

这份文档描述当前可以实际运行、可审计、标签无关的决策层。它不是声称 FMO 端点能量已经等价于完整自由能，而是规定：面对一个新 target，哪些结构处理可以自动做，哪些 QM 臂必须完整计算，哪些结果只能标记为 review。

## 1. 固定原则

1. **先结构，后能量。** 原始 PDB/SDF、氢、形式电荷、键级、蛋白受体、水和 ligand pose 先冻结并哈希；结构有严重非水重叠时 fail-closed，不静默移动重原子。
2. **target-level 统一。** 水、受体 pocket、几何 restraint 和 FMO 参数在一个 target 内统一；不逐配体根据实验结果挑 `k`、水版本或 score。
3. **所有候选分支并行保存。** 每例保留 `tPIE`、`ES`、`ES+SOLV`、`ES+CT+SOLV`。其中 `ES+SOLV` 是 PIEDA 分解中的子集，不能与 `tPIE` 相加。
4. **稳定不等于准确。** 多臂排序一致性只衡量数值/结构敏感性；Galectin3 证明稳定分支仍可能整体方向错误。
5. **熵不收敛就不入分数。** 当前 E61 的 96 blocks 未通过 ESS/分半 JS 门槛，故只能报告诊断，不加入 FMO ranking。

## 2. 自动阶段

### A. 输入和化学审计

- 校验 ligand 数量、唯一 ID、显式氢、元素、键级、形式电荷、立体化学和 SDF/PDB 哈希。
- 校验蛋白残基、链、二硫键、肽端封端、辅助因子和水的 identity。
- 把 pH/质子化不明确标为 `REVIEW`，不在单个配体上静默猜 microstate。

### B. 几何审计和修复

- 计算 ligand–receptor 最短接触和 vdW ratio；水冲突可以生成 target-level 去水对照。
- 只可自动重新生成 ligand 氢而保持重原子；严重非水接触不自动“修好”。
- 固定蛋白，仅 ligand relaxation，登记 `k=1/10/100`；所有臂都评分。
- 由源结构水几何门触发时，整个 target 使用同一个 union-water receptor；wet/dry 只能作为预注册对照。

### C. QM 完整性

- 每个 target 必须有完整 ligand 覆盖、统一 receptor hash、统一 FMO2/PCM/PIEDA deck 和成功的 SCF/PIEDA reconciliation。
- 有任何缺失或失败就不能生成完整 ranking，只输出 `FMO_FALLBACK` 和失败清单。
- `tPIE` 使用 selected protein–ligand PIEDA pair；不使用 full-log 中巨大 EX/CT 抵消项的宽格式误解析。

### D. 采样和熵

- 多构象必须固定生成数、seed、聚合规则，并报告每个构象；不得按实验标签 cherry-pick。
- 若做 bound/free entropy，要求同一原子身份/定义、完整轨迹、ESS、block variance、split-half JS 和重复 replica 同时通过。
- 未通过时延长预注册 block；达到预算仍不通过则 `OMIT`，不能把最后一段当成收敛自由能。
- `scripts/compile_fmo_target_entropy_plans.py` 根据每个 target 的蛋白原子数、最大配体可旋转键数和最大配体重原子数生成无标签采样预算，输出 `runs/fmo87/ablation/final_analysis/fmo_target_entropy_plans.json`。当前六个 target 均落在 low tier，说明这套预设复杂度门不需要人为制造不同参数；真正的自适应发生在 ESS/分半 JS 不通过时追加独立 blocks，而不是按实验结果调参。

## 3. 决策状态

| 状态 | 含义 | 可以输出 |
|---|---|---|
| `FMO_READY` | 结构、协议、QM 全覆盖，姿态方向有独立支持 | 冻结的 `ES+SOLV` 主 ranking；同时输出 tPIE/ES/CT/SOLV 诊断 |
| `FMO_REVIEW` | QM 完整但姿态方向或结构解释仍需独立验证 | 各候选分支和不确定性；不能宣称生产准确度 |
| `FMO_FALLBACK` | 标准受体臂缺失、QM 不完整或关键协议门失败 | Boltz2/其它独立模型作为单独结果；不能混加 FMO kcal/mol |
| `UNCERTAIN_ENTROPY` | 结构/QM 可用但熵未收敛 | 维持无熵 FMO 分数，熵只作诊断 |

## 4. 当前回放状态

统一机器报告：`runs/fmo87/ablation/final_analysis/fmo_universal_decision_report.json`。

- CDK2、JNK1、P38、TYK2：当前冻结原 87 协议下为 `FMO_READY`；主分数仍是预注册的 `ES+SOLV`，熵 `OMIT`。
- P38：源水几何门触发 target-level dry19 对照；dry19 `k=10` 使 tPIE 上升，但没有证明所有分支都会上升，因此 wet/dry 不按配体切换。
- Galectin3：修复后 26/26 QM，并通过独立 Boltz2 几何门；当前能量方向仍未建立，因此 `FMO_READY` 只表示结构/QM 可部署，不表示亲和力准确。
- FXa：E92 union-pocket 27/27 QM 完成，并通过独立几何门；8m 是 case-level 几何 outlier，GLN61 排除仍属于结构触发规则。旧 reference-pocket 结果继续只作诊断。

## 5. AnewFEP 的正确迁移

AnewFEP v2 的可迁移点是“自动化 + 误差源分解 + 收敛/稳定性质量报告”：局部力场参数、蛋白构象响应、rotamer/采样陷阱和数值稳定性分别做对照。不能照搬其 REST2、λ windows 或 RBFE 公式到静态 FMO 端点分数。对 FMO 的对应做法是：结构诊断触发统一对照臂；每个臂完整 QM；通过协议哈希、SCF、PIEDA 闭合、构象分布和收敛指标决定是否可信。

## 6. 不能做的事情

- 不将 `tPIE + ES+SOLV` 当成新物理量；这是重复计算。
- 不因 Galectin3 或 FXa 的历史相关系数而后验选择 score、water、pose 或 `k`。
- 当前统一执行协议由 `scripts/compile_fmo_universal_adaptive_workflow.py` 编译到 `runs/fmo87/ablation/final_analysis/fmo_universal_adaptive_workflow.json`；它只读取结构、化学、数值稳定性和收敛审计，不读取实验亲和力来选分支。
- 原 87-case 因素消融总表为 `runs/fmo87/ablation/final_analysis/fmo_factor_effect_summary.json`；它只用于解释 tPIE、PIEDA、PCM、PBSA、relax 和熵 proxy 的作用，禁止据此在同一数据集上后验换分数。
- 当前跨集审计为 `runs/fmo87/ablation/final_analysis/fmo_cross_target_generalization_audit.json`；它固定比较原始87、Galectin3 E68 和 FXa E92，并额外穷举四个原子分支的 15 个等权 robust-z 子集。当前只有单独 tPIE 在三组数据方向均为正，但 Galectin3 的 Kendall 仅 `0.174`，故这是保守方向候选而不是准确性证明。
- 不因某一子集达到 `0.55/0.38` 就声称新 target 普适；必须在冻结规则下完成外部验证。
- 不用未收敛熵、不完整 QM 或旧 receptor pocket 生成“完整 87/27 case 排名”。
