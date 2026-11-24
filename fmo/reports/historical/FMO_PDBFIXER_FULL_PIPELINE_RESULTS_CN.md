# PDBFixer→PDB2PQR→水处理→k10/k100→FMO 外部验证

日期：2026-10-05。

## 结论

这次实验验证了完整结构制备链，但没有发现 PDBFixer 会改变 Galectin3 或 FXa 的有效 FMO 输入：两个原始蛋白都没有缺失残基、缺失重原子或非标准残基。PDBFixer 仅对 Galectin3 增加 1 个氢；FXa 的原始结构已经包含相应氢信息。经过 PDB2PQR/PROPKA 和统一 target-level 水处理后，受体原子序列与坐标和已有 FMO 生产臂逐原子一致。

因此没有重复提交完全相同的 GAMESS deck；已有 QM 输出可作为该完整流程的严格等价结果，而不是把重复计算误报为新的独立实验。

## 固定流程

```text
原始 PDB
→ PDBFixer 缺失结构审计与必要补氢
→ PDB2PQR AMBER/PROPKA pH 7.4
→ 配体显式氢和形式电荷审计
→ target-level 冲突水处理
→ OpenMM/OpenFF 固定受体 k10 或 k100
→ 统一 pocket/FACIO-like fragment/ACE-NME/HOP-C
→ FMO2-DFTB3/PCM selective PIEDA
```

PDBFixer 审计脚本为 `scripts/prepare_fmo_pdbfixer.py`。所有输出保存在 `runs/fmo87/ablation/E93_pdbfixer_full_pipeline/`。

## 结构审计

| target | 原始原子 | 缺失残基 | 缺失重原子 | 非标准残基 | PDBFixer 后原子 | 水处理 |
|---|---:|---:|---:|---:|---:|---|
| Galectin3 | 2702 | 0 | 0 | 0 | 2703 | 删除 A:403/410/475/493 四个 target-level 冲突水 |
| FXa | 5004 | 0 | 0 | 0 | 5004 | 保留统一受体水策略 |

PDB2PQR 后，Galectin3 去水的原子序列和坐标与 E68 完全一致；FXa 与 E40 完全一致。配体重原子姿态没有被移动。FXa 输入已经具有完整显式氢，因此不使用会改变原子顺序的重复氢重建；这是为了避免把“修复工具的原子重排”误当成结构物理效应。

## FMO 排序结果

指标格式为 `Pearson / Kendall`，均为无拟合、完整队列结果。

| target/流程 | n | tPIE | ES | ES+SOLV | ES+CT+SOLV |
|---|---:|---:|---:|---:|---:|
| Galectin3 PDBFixer+PDB2PQR+dry+k10 | 26 | `0.2381 / 0.1736` | `-0.4922 / -0.4489` | `-0.4338 / -0.3932` | `-0.4140 / -0.3039` |
| Galectin3 PDBFixer+PDB2PQR+dry+k100 | 26 | `0.2001 / 0.2384` | `-0.5587 / -0.4985` | `-0.4993 / -0.4118` | `-0.4637 / -0.3684` |
| FXa PDBFixer+PDB2PQR+k10 | 27 | `0.4728 / 0.3823` | `0.5280 / 0.3937` | `0.5293 / 0.3766` | `0.5094 / 0.3823` |
| FXa PDBFixer+PDB2PQR+k100/pocket10 | 27 | `0.1193 / 0.1255` | `0.4720 / 0.3937` | `0.4501 / 0.3994` | `0.4319 / 0.3823` |

## 与原有结果的比较

- **PDBFixer 本身没有改善排序。** 它没有发现需要补重原子的结构缺陷，因此不会带来新的几何信息。
- **Galectin3：** k10 的 tPIE 为 `0.2381/0.1736`，k100 为 `0.2001/0.2384`；ES、ES+SOLV 和 ES+CT+SOLV 仍整体反向。主要问题不是“缺失残基/重原子”，而是姿态、受体状态、水处理或该物理指标与实验亲和力的适配性。
- **FXa：** k10 的 tPIE 为 `0.4728/0.3823`，ES+SOLV 为 `0.5293/0.3766`；k100/pocket10 下 tPIE 明显下降到 `0.1193/0.1255`，说明 FXa 对配体约束/受体口袋定义更敏感，但结构完整性修复不是主要瓶颈。
- **当前 best：** Galectin3 在这些结构修复分支中仍是 k10 tPIE 的弱正向；FXa 的完整生产分支仍应优先使用现有 union-pocket k10/E92 结果 `tPIE=0.6219/0.4051`，而不是 k100/pocket10 诊断臂。

## 计算闭合性

- Galectin3 k10/k100：已有完整 `26/26` QM 输出；对应 E68/E70 受体输入与本实验去水后的受体逐原子一致。
- FXa k10/k100：已有完整 `27/27` QM 输出；对应 E40/E47 受体输入与本实验 PDBFixer/PDB2PQR 结果逐原子一致。
- 由于 deck 输入相同，重新跑一遍只会产生数值舍入/机器运行差异，不能作为新的结构处理对照；本实验将输入闭合和输出复用关系全部记录，而不是重复消耗 QM 资源。

## 后续规则

1. 以后所有新 target 先跑 PDBFixer 缺失结构审计；无缺失重原子时不强制改写蛋白。
2. 需要补重原子/残基时，必须把 PDBFixer 输出作为新结构分支，重新执行 PDB2PQR、质子化审计、水审计、k10/k100 和完整 FMO。
3. PDBFixer 不得按单个 ligand 后验修复；结构修复和水处理必须 target-level 统一。
4. PDBFixer 的“结构通过”不等于排序会改善；本实验明确显示，结构完整性通过与预测准确性是两个独立门。
