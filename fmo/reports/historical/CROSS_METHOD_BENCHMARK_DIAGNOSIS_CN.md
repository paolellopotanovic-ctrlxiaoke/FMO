# FMO 与同 benchmark 物理基线对照

日期：2026-10-06。

## 核心结论

Galectin3 的失败不是这个数据集上所有物理方法的共同失败。官方 benchmark 中同一批 26 个 Galectin3 配体的 FEP+ 结果为 Pearson/Kendall `0.874/0.661`，而我们完整结构审计后的 FMO 最好仍是 tPIE `0.238/0.174`，ES、ES+SOLV 和 ES+CT+SOLV 均为负相关。因此 Galectin3 需要继续 debug FMO 的姿态/电荷/能量定义，不能用“物理方法都不适合这个集合”解释。

FXa Yoshikawa 27 例上，FMO 和 FEP+ 都是正向：FMO tPIE `0.622/0.405`，FEP+ `0.751/0.508`。在其它 FXa benchmark 子集上，两种方法的方向也大体一致：set3 两者都反向，set4 两者都正向，set6 两者都较强正向；set5 只有 2 例，不能作为有效验证。因此 FXa 的结果更符合“集合/结构系列本身难度不同”，而不是 FMO 单独崩溃。

## 同数据集指标

指标格式为 `Pearson / Kendall`；FMO 四列均为无拟合、同一套 FMO2-DFTB3/PCM selective-PIEDA 流程。FEP+ 为仓库 `21_4_results/ligand_predictions/` 中随 benchmark 发布的官方结果。

| 数据集 | n | FMO tPIE | FMO ES | FMO ES+SOLV | FMO ES+CT+SOLV | 官方 FEP+ |
|---|---:|---:|---:|---:|---:|---:|
| Galectin3 extra | 26 | `0.238/0.174` | `-0.492/-0.449` | `-0.434/-0.393` | `-0.414/-0.304` | `0.874/0.661` |
| FXa Yoshikawa | 27 | `0.622/0.405` | `0.652/0.399` | `0.643/0.371` | `0.621/0.371` | `0.751/0.508` |
| FXa set3 | 4 | `-0.879/-0.548` | `-0.935/-0.913` | `-0.948/-0.913` | `-0.912/-0.913` | `-0.428/-0.548` |
| FXa set4 | 11 | `0.724/0.455` | `0.754/0.527` | `0.710/0.491` | `0.707/0.491` | `0.646/0.455` |
| FXa set5 | 2 | `1.000/1.000` | `1.000/1.000` | `1.000/1.000` | `1.000/1.000` | `1.000/1.000` |
| FXa set6 | 6 | `0.957/0.867` | `0.862/0.733` | `0.892/0.867` | `0.890/0.733` | `0.903/0.333` |

set3、set5 的样本数太小，相关系数只作诊断，不能和 26/27 例结果等量齐观。

## 结构输入诊断

扩展 FXa 子集时，原始 benchmark PDB 本身暴露了明显的结构异质性：

- `fxa_set3`：PDBFixer 发现 THR244 缺失 6 个重原子和末端 OXT；PDB2PQR 原始流程无法修复，必须先 PDBFixer 补全。
- `fxa_set5`：PDBFixer 发现 GLN 的 OE1/NE2 缺失；同样需要先修复。
- `fxa_set4/6`：PDB2PQR 会输出异常的缺失行/水链拓扑，直接进入 OpenMM 时出现 HOH 模板错误；统一改用 PDBFixer 后去水的 label-free 分支才可计算。
- 四个子集最终均通过统一的 PDBFixer→去原始氢→PDB2PQR/PROPKA→无水 k10→FMO 流程，QM 成功率 `23/23`。
- 因此 set3 的反向结果不能单独解读成“FMO 物理失败”：其原始结构确实有缺失重原子；但 set4、set6 在相同修复流程下表现良好，说明该流程仍能在部分系列上工作。

## 对 Galectin3 的具体判断

Galectin3 与 FXa 的差异提供了最强的诊断对照：

1. Galectin3 已完成 PDBFixer/PDB2PQR/冲突水审计，26/26 QM 成功；不能再把主要问题归为缺失残基或 QM 不收敛。
2. Galectin3 k10/k100 的 FMO tPIE 都只是弱正向，而所有 ES 类三态/PIEDA 子集都反向；改变 k 不会恢复方向。
3. 同一 26 例官方 FEP+ 很强，说明实验标签、系列结构和基本构象信息并非完全不可预测。
4. 但 FEP+ 使用自己的力场、自由能路径、采样和质子化/构象处理，不能据此直接证明 FMO 的晶体姿态一定正确。它只能证明“存在一套物理计算流程可以在这组数据上工作”。

目前最合理的根因排序是：

| 候选根因 | 当前证据 | 判断 |
|---|---|---|
| 原始蛋白缺失重原子/残基 | Galectin3 已通过 PDBFixer 审计；FXa 主集也已修复闭合 | 不是 Galectin3 主因 |
| 单一 k10 参数 | Galectin3 k1/k10/k100 都没有恢复三态项方向 | 不是唯一主因 |
| 晶水冲突 | 去掉 4 个 target-level 冲突水后 tPIE 从更差变成弱正向 | 是污染源，但不是完整根因 |
| FMO 姿态/局部几何与 FEP+ 处理不同 | FEP+ 同集明显更好，FMO 多个能量项方向反向 | 高优先级待 debug |
| FMO 的 PIEDA 项与结合自由能不等价 | tPIE 弱正向但 ES/溶剂相关项反向 | 高优先级待 debug |
| 数据集物理方法普遍失效 | FXa set4/set6、Galectin3 FEP+ 均反例 | 已基本排除 |

## 下一步 debug 顺序

1. 对 Galectin3 逐配体对齐 FMO k10 pose、原始晶体 pose 和 FEP+ 输入 pose，比较重原子 RMSD、关键接触、氢键和芳香堆积；当前 k10 最小化的配体 RMSD 中位数约 `0.312 Å`，最大位移约 `1.37 Å`，所以不能仅凭“移动很小”排除局部接触错误。
2. 做 Galectin3 的 FMO-vs-FEP+ 同构象对照：固定同一套重原子姿态，只改变能量模型；这是区分“构象错误”和“FMO 能量定义错误”的关键实验。
3. 检查 Galectin3 每个配体的形式电荷、互变异构体、配体显式氢和 PDB2PQR 受体局部质子化，尤其针对 FMO 与 FEP+ 排名分歧最大的配体。
4. 在不读取实验标签的前提下，对关键接触做统一多构象采样；所有预注册构象都纳入，不能挑选最优构象。
5. 不要继续优先调 ES/SOLV 固定权重；Galectin3 已显示问题主要不是单一线性组合，而是姿态/微观态/模型边界的交互。

## 可复算产物

- 统一比较脚本：`scripts/analyze_cross_method_benchmark.py`
- 汇总 JSON：`runs/fmo87/ablation/E94_fxa_other_benchmark_sets/cross_method_summary.json`
- FXa 子集 FMO 输出：`runs/fmo87/ablation/E94_fxa_other_benchmark_sets/fixer_pipeline_runs/`
- 官方 FEP+ 输入结果：`vendor/public_binding_free_energy_benchmark/21_4_results/ligand_predictions/`
