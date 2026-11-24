# OpenBind FMO 失败结构根因与自动流水线

日期：2026-10-08。

## 结论

OpenBind 的初始 639 个结构中，最初有 4 个 FMO 失败。它们不是同一种问题：

| case | 初始问题 | 结构层证据 | 修复结果 |
|---|---|---|---|
| `A71EV2A-x2513a` | GAMESS `Distance errors` | FMO 报告片段间异常近距；结构审计发现 ligand/protein 局部重叠 | 自由受体 clash-repair 不能完全解决；降低配体约束后成功，tPIE `-47.386` |
| `A71EV2A-x4445a` | 二体 SCF 不收敛 | 原始 OpenMM 势能约 `1.39×10^7 kJ/mol`，最短 ligand/protein 接触约 `1.49 Å` | 允许受体自由受限修复后成功，tPIE `-64.298` |
| `A71EV2A-x5836a` | 二体 SCF 不收敛 | 带正电配体，原始局部几何和电子态更难收敛 | 自由受限修复后成功，tPIE `-132.999` |
| `A71EV2A-x6838a` | 标准臂二体 SCF 不收敛；修复臂出现强离子对 | 带正电配体与 ASP39 形成约 2.7–2.8 Å 盐桥；PIEDA 约 −114 kcal/mol 项与异号电荷一致 | `k10_free_repair` 成功，tPIE 约 `−222.848 kcal/mol` |

因此，不能把“PDBFixer 通过”理解成“FMO 一定可收敛”。PDBFixer 解决的是蛋白拓扑完整性；它不负责修复晶体配体姿态与受体之间的重叠，也不保证带电 FMO 二体 SCF 收敛。

## 为什么其它结构能跑、这些结构不行

### 1. 蛋白拓扑不是主要差异

四个 case 的 PDBFixer 审计均未发现内部缺失残基、缺失侧链重原子或非标准残基。只报告末端 `OXT`。这说明失败不能归因于“PDB 缺了一个普通重原子”。

### 2. 局部接触异常是第一类根因

`x4445a` 的标准准备结果初始势能约 `1.39×10^7 kJ/mol`，并存在约 `1.49 Å` 的 ligand/protein 接触。固定受体的 `k=10` relax 不能可靠地消除这种冲突，FMO 会把异常局部几何直接带入片段单体/二体计算。

`x2513a` 甚至在 GAMESS FMO 几何检查中明确报告 `Distance errors`。这类错误不能靠调整评分公式、PCM 或排序权重解决，必须先回到结构准备层。

### 3. 带电配体会放大 SCF 难度

`x5836a` 和 `x6838a` 是少数带正电配体。带电配体并不意味着输入错误，但它会使局部静电场、PCM 极化和二体 SCC-DFTB 收敛更敏感。`x5836a` 经结构修复后可以收敛。`x6838a` 的标准臂确实不收敛，但修复臂成功后发现原先“所有超过 100 kcal/mol 的 PIEDA 项都失败”的规则误杀了异号盐桥；现在使用电荷感知审计，保留该强离子项并记录为 `strong_ionic_pairs`。

### 4. `k=10` 不是万能的 clash repair

当前标准 `k=10` 是为了尽量保留实验结合姿态、只做温和结构松弛。它不是无约束 docking，也不是强制消除所有重叠的几何优化。因此：

- 正常结构：使用固定受体/水的标准 `k=10`；
- 发现严重非水碰撞：进入预先注册的自由受体受限修复臂；
- 修复后必须重新审计；
- 修复仍失败：进入 `REVIEW_REQUIRED`，不静默填数。

## 固化的自动流程

```text
原始晶体 PDB/SDF
→ PDBFixer 拓扑审计
→ 必要时补缺失残基/重原子/末端原子并补蛋白氢
→ RDKit 配体显式补氢，保留晶体重原子坐标
→ PDB2PQR/PROPKA，AMBER，pH 7.4
→ 配体形式电荷和质子化状态审计
→ 水审计；只按预注册 target-level 规则处理水
→ 标准固定受体/水 OpenMM/OpenFF k=10
→ ligand-protein clash 审计
→ FMO2-DFTB3/PCM/PIEDA
→ 若失败，按固定顺序进入修复臂
→ 修复后重新 PDB2PQR/最小化/审计/FMO
→ 仍失败则 REVIEW_REQUIRED，不进入无缺失最终排序
```

自动工具：`scripts/run_openbind_complete_pipeline.py`

固定修复顺序：

1. `k10_fixed`：标准协议；
2. `k10_free_repair`：允许受体/水自由受限松弛，配体约束仍为 `k=10`；
3. `k0p1_free_repair`：仅作为结构冲突恢复诊断，不能与标准结果混称；
4. 仍失败：写入 `pipeline_provenance.json` 并标记 `REVIEW_REQUIRED`。

工具不读取实验亲和力来选择 arm，也不按 tPIE 好坏 cherry-pick。批处理器也已修复：只有 `result.json.status == SUCCEEDED` 才跳过，失败结果会自动重试。

## 当前实现状态

已完成：

- PDBFixer 兼容性修复：兼容当前环境中 `missingAtoms` 返回字符串名称的 API 行为；
- 自动完整流程工具：`scripts/run_openbind_complete_pipeline.py`；
- 失败结构根因报告：`runs/openbind/openbind_failure_root_cause.json`；
- 失败结果可重试：`scripts/run_openbind_fmo_batch.py`；
- 自动流程测试：`tests/test_openbind_complete_pipeline.py`；
- 相关测试通过：`16 passed`。

当前 OpenBind 结果覆盖：

- 原始固定标准协议为 `638/639` 结构成功；自动补强流程已使 `x6838a` 通过修复臂，修复后可达到 `639/639`，但需在总 benchmark 中单独标记 fallback；
- `489/490` compound 成功；
- `A71EV2A-x6838a` 已由自动工具选定 `k10_free_repair` 并成功完成 FMO；标准 `k10_fixed` 的失败记录仍保留；
- 不能把它从最终 benchmark 静默删除，也不能填入估计 tPIE。

## 最终报告规则

必须同时报告三种覆盖率：

1. **标准协议覆盖率**：完全使用 `PDBFixer → PDB2PQR → k=10 fixed → FMO`；
2. **修复协议覆盖率**：使用预注册结构修复 arm 后成功；
3. **人工复核覆盖率**：结构或 QM 仍失败，不能进入无缺失生产排序。

这比单独给一个“成功率”更可靠，因为它区分了统一标准结果和结构修复 fallback，避免把不同几何协议混成一个方法。

### `x6838a` 最终诊断

这不是“PDB 缺少普通重原子”：PDBFixer 只发现末端 OXT，并完成蛋白氢和末端拓扑补全。标准 `k10_fixed` 的二体 SCF 不收敛；`k10_free_repair` 后，配体阳离子 N 与 ASP39 氧原子最近约 2.69 Å，PIEDA 的 −113.866 kcal/mol 项正对应异号盐桥，并与 GAMESS selected TIE（−222.8502）相差约 0.002 kcal/mol。因此修复的是质控规则与 fallback 的坐标一致性，不是篡改能量。

## 2026-10-09 全量补强版结果

全量 `639` 个结构均已尝试；目前 `638/639` 成功、`489/490` 个化合物可用。按同化合物多结构算术平均、以 `-tPIE` 对实验 pKD 计算：Pearson `0.3734`、Spearman `0.4437`、Kendall `0.2902`；详情在 `runs/openbind/complete_pipeline_full_metrics.json`。分臂：标准固定 k10 `635`、自由受体 k10 修复 `2`、固定受体 k0.1 修复 `1`。与旧版 Spearman `0.4431` 相近，但两版缺失的化合物不同，不能仅凭总表差值声称物理改善。

剩余 `A71EV2A-x5836a` 在全部自动修复臂仍二体 SCF 不收敛。没有 PDBFixer 的旧流程曾成功计算该例，提示结构准备/补氢导致的局部变化值得审计，但尚不能证明具体根因。不能把旧协议能量填充到新协议中；本机 GAMESS 的 `$CONTRL MAXIT` 上限是 200，500 属于无效设置。当前**没有完成 639/639 的无缺失验证**，该例仍需复核。

## 2026-10-09 后续归因实验与完成状态（取代上一节临时结论）

受体逐原子对照发现：旧 PDB2PQR 与新 PDBFixer/PDB2PQR 结果仅 TYR15 的 HH（相差约 0.49 Å）及远处链末 GLU145 的 OXT（相差约 0.91 Å）坐标不同；所有其它受体原子坐标相同。配体输入 SDF 相同，但独立 OpenMM 最小化得到的配体姿态整体 RMSD 约 0.19 Å、最大单原子偏移约 0.60 Å，所选口袋由 31 片段变为 32 片段。这些**不证明单个氢是唯一原因**：交叉实验固定旧口袋、换新配体姿态仍 SCF 失败；固定新口袋、换旧姿态也失败；旧口袋与旧姿态同时使用并保持新制备受体则成功（−132.999 kcal/mol）。

再做严格一致的补充实验：从**原始无 PDBFixer PDB**重新运行 PDB2PQR → 自由受体 k10 最小化 → **受体和配体同时取松弛后的坐标**做 FMO，而不是挪用旧 QM 数值。该计算独立收敛，tPIE 约 −195.951 kcal/mol。已加入一个最后才触发的无标签 `original_source_k10_free` 自动修复臂，重新自动运行该例得到收敛 tPIE 约 −196.465 kcal/mol；微小差异来自两次独立最小化的浮点/收敛路径，不应视为完全重现同一能量。

最终流程覆盖 **639/639 个结构、490/490 个化合物**；臂数：标准 PDBFixer + 固定受体 k10 `635`，PDBFixer + 自由受体 k10 `2`，PDBFixer + 固定受体 k0.1 `1`，原始受体 + 自由受体 k10 `1`。所有结果均有 arm 溯源；这是**确定性分级 fallback 方法**，绝非“639 个都遵守完全相同的固定受体/PDBFixer 制备协议”。无标签选臂、不读取实验亲和力。

完整 490 化合物以结构等权平均后：Spearman `0.446519`、Pearson `0.373929`、Kendall `0.292338`（方向为 `-tPIE` 对实验 pKD）。与之前 `0.443083`（489 个化合物，且缺失化合物不同）不能直接归因为排序性能改善。数据文件为 `runs/openbind/complete_pipeline_full_metrics.json`，完整失败/选臂/质控追溯文件为 `runs/openbind/complete_pipeline_full/complete_pipeline_summary.json`。
