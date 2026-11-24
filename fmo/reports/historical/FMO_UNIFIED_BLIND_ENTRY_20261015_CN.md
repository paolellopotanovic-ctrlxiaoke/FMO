# 统一盲测入口和故障边界（2026-10-15）

使用 **`scripts/fmo_workflow.py blind-run`** 作为新靶点一命令入口：内部依次调用 `intake` → `geometry` → `adaptive-run`。输入为无实验亲和力的配体表、受体 PDB、晶体坐标 SDF。`intake` 无标签选择相接触的蛋白链；`geometry` 对整靶点审计并统一处理严重冲突水。`adaptive-run` 对每个输入执行 PDBFixer、重置氢后 PDB2PQR/PROPKA、固定受体 k10、默认肽链编号修复、全配体 5 Å 共同口袋、FMO2-DFTB3/PCM 与 tPIE 质控。严格计算全部成功时直接输出严格结果，任何例数缺失均拒绝排序。

调用示例（路径需替换为实际输入和安装目录）：

```bash
python scripts/fmo_workflow.py blind-run --target TARGET --ligands cohort.csv \
  --protein receptor.pdb --poses crystal_poses.sdf --output runs/new_target \
  --gamess-dir software/gamess --library-dir /path/to/library \
  --dftb-parameter-dir software/gamess/auxdata/DFTB/3OB-3-1 \
  --workers 4 --ncores 4
```

默认无须事后手动重跑；`run` 仅是**严格臂/诊断入口**，不会启动回退，不应单独用作“新数据自动完成”的保证。

如果单晶体单配体 QM 失败，沿用已登记的逐例 SCF、重原子口袋和几何救援，并保留原计算日志。若同靶点**多配体**严格臂 QM 失败，先检查至少两个独立失败例：GAMESS 日志第一次失收敛时的 monomer 序号必须一致；在冻结碎片表中该序号映射到同一个非水残基；该残基仅因含氢距离而进入 5 Å 口袋，其最近重原子距离大于 5 Å。满足这些**无标签**条件时，冻结原来通过 PDBFixer→PDB2PQR、统一水审计与 k10 的整靶点受体和配体姿态，仅从**所有配体共享的 QM 口袋**剔除这个诊断片段，再对全队列重算。不会从受体 PDB 中删残基、不会把 FMO 变成不含氢计算。若无法定位相同单体，或最小修复没有让全部计算通过，才按已审计的重原子口袋差异尝试整批重原子口袋回退；否则拒绝排序。**优先级只看故障定位和结果覆盖，不看实验相关系数。**

## 完整运行验收

| 数据集 | 完整计算覆盖 | tPIE Pearson / Spearman / Kendall | 协议 |
|---|---:|---:|---|
| 原 87，四靶点内按样本量加权 | 87/87 | 0.6352 / 0.6272 / 0.4635 | 默认编号、严格臂全部通过 |
| OpenBind | 639/639 晶体、490/490 化合物 | 0.3791 / 0.4437 / 0.2905 | 默认编号；636 严格、1 重原子口袋救援、2 自由 k10 + DIIS 救援；每个化合物多个晶体取均值 |
| Galectin3 原始蛋白与晶体姿态，完整默认 `blind-run` | 26/26 | **0.56093 / 0.54567 / 0.38700** | 同一 `intake→geometry→PDBFixer→PDB2PQR→k10→FMO`；QM 全部严格通过；原始 SDF 实验属性已剥离，仅在算完后读取实验标签 |
| FXa 原始受体；旧的宽口袋回退 | 27/27 | 0.5450 / 0.5483 / 0.3595 | 默认编号；严格含氢口袋 27 例均 SCF 失败；重原子口袋选中 358 个原子 |
| FXa 原始受体；**现固化的故障片段优先臂** | 27/27 | **0.5833 / 0.5639 / 0.3880** | 失败日志共同定位 `monomer 4 = A:GLN:61`；仅从 450 原子 QM 口袋剔除该残基，保留其余口袋成员共 433 原子 |
| FXa 原始受体 A+B 二硫键自动联链；QM 选择性修复 | 27/27 | **0.58345 / 0.56694 / 0.39372** | 原始 A 链和 B 链由 A122–B132 二硫键连接；原始 188 个水均保留；GLN61 只从 433 原子 QM 口袋剔除 |

原 87 文件 `runs/fmo87/numbered_main_pipeline_20261015/complete_87_tpie.json`；OpenBind 文件 `runs/openbind/numbered_main_pipeline_20261015/complete_490_adaptive_tpie.json`；FXa 历史宽口袋回退文件 `runs/fmo87/numbered_main_pipeline_20261015/FXa/complete_27_adaptive_tpie.json`，新片段优先验证文件 `runs/fmo87/numbered_main_pipeline_20261015/FXa/complete_27_fragment_tpie.json`。旧的 `adaptive_status.json` 是运行时已冻结的宽口袋结果，**不能追溯改写**；新臂有独立 `fragment_repair_status.json` 和输出目录。现有数据复用了严格臂受体/k10 和此前已通过试点的 QM 缓存，尚未清空缓存从零执行一次新的默认 `blind-run`。原 87/OpenBind 没有该整靶点失败，计算分数不因添加回退规则改变。三套数值不可跨靶点直接合并为一个相关系数，原 87 为靶点加权、OpenBind 以 compound 为单位、FXa 为单靶点。参照标签仅供完整盲结果完成后计算统计量。

单晶体与整靶点回退的选择依据只包含结构、运行错误和冻结输入，不读取实验标签。此流程只能保证有足够可信输入和已知故障类型时的确定性处理；不能保证任意坏结构都能得到有效的 QM 排序。FXa 最新结构审计发现此前 A-only 流程**漏掉通过 A122–B132 二硫键连接的 B 链**：老 E92 受体有 4,972 个原子、188 个水，A+B 新受体有 4,989 个原子、188 个水。二者差异主要是老 E92 整体缺失连续肽链中的 A:GLN:61（17 个原子），但水氢、其他残基的氢/质子化也有变化；此前拿 A-only 的 145 个水归因于 PDBFixer 去水，是错误对照。旧 E92 Pearson 0.62190，最新 A+B 无标签最小 QM 口袋修复 Pearson 0.58345。**不能以实验指标更好为由直接删除连续肽链中的整个 GLN61**：TYR60–GLN61–ALA61A 的结构是连续的，旧 E92 缺 GLN61 时断开了主链。应先验证具有完整主链和一致质子化/水氢的可迁移受体处理规则，再按预冻结规则运行；历史受体只能作为诊断对照，不可伪称相同的完整结构制备。FXa/GLN61 的诊断规则是在已知该问题的历史案例上形成的，未来新靶点有效性仍需独立无标签外推检验。

### FXa 受体 × 姿态交叉消融（全部 27/27）

| 已冻结受体 | 历史 E92 姿态：Pearson / Spearman / Kendall | 新 A+B 流程姿态：Pearson / Spearman / Kendall |
|---|---:|---:|
| 历史 E92 受体（缺 GLN61） | 0.621897 / 0.582812 / 0.405136 | 0.621574 / 0.582812 / 0.405136 |
| 新 A+B 受体（保留 GLN61，仅 QM 口袋剔除） | 0.585941 / 0.566936 / 0.393724 | 0.583453 / 0.566936 / 0.393724 |

相同受体换配体 k10 姿态的 Pearson 差 <0.003；跨受体切换约 0.038。口袋附近的蛋白重原子与水氧坐标完全相同，但 5 Å 内有 10 个水氢和 TYR99 酚羟基氢移动超过 0.2 Å。**推断**：此差距以受体制备/质子化/水氢与 GLN61 的联动差异为主，不能将其确认为单一的 GLN61 效应；已经启动在新完整主链受体上分别仅变更水氢或 TYR99 氢的回顾性全 27 例实验，不能按已知标签挑 arm。完整输入/预测 SHA 与等长逐例核对见 `runs/fmo87/ablation/fxa_receptor_pose_factorial_20261015.json`，复算入口 `scripts/evaluate_fxa_receptor_pose_factorial.py`。旧回退 `adaptive_status.json` 保持不可改写；同一严格臂修复路径误判后另存的 `fragment_repair_status.json` 证明 27/27 选择性臂可自动完成两阶段，但**不是新独立冷启动默认入口验收**。默认 FXa 冷启动在 `runs/fmo87/unified_disulfide_fxa_default_retest_20261015`，尚在运行。Galectin3 已由原始受体、无实验属性的晶体坐标一次 `blind-run` 产生 26/26 严格 QM 成功；经 `scripts/evaluate_galectin3_unified_blind.py` 核对原始结构、全部 QM 输入/日志 SHA 和 k10 再解盲，正式指标在 `runs/fmo87/unified_galectin3_default_20261015/complete_26_default_tpie.json`；不能因 Gal 完成而声称四数据集的同版本冷启动都已完成。

### FXa 水氢差异的全量单因素定位与下一版默认入口

输入 FXa 晶体 PDB **已有全部 188 个水的 O/H1/H2**，历史 E92 的水氢坐标与原始 PDB 在逐原子精度完全相同；本轮 PDB2PQR 重置所有氢后改变了约 278 个水氢的方向，其中共用 QM 口袋内有 14 行水氢坐标不同。固定同一个 27 配体 k10 姿态、完整受体、QM 剔除同一失败 GLN61、相同非氢原子和 GAMESS 参数时：

| 更改项（27/27 同配对病例） | Pearson | Spearman | Kendall |
|---|---:|---:|---:|
| PDB2PQR 重定向水氢（此前默认） | 0.583453 | 0.566936 | 0.393724 |
| 仅恢复原始晶体 PDB 水氢方向 | **0.618882** | **0.582812** | **0.405136** |
| 仅替换 TYR99 酚氢（不改变水） | 0.586973 | 0.566936 | 0.393724 |
| 仅替换 TYR225 酚氢（不改变水） | 0.583486 | 0.566936 | 0.393724 |
| 源水氢 + 两个 TYR 酚氢全部恢复 | 0.621574 | 0.582812 | 0.405136 |
| 无标签 AMBER14/TIP3P 固定全部重原子、全受体氢最小化 | 0.620108 | **0.609678** | **0.427960** |
| 历史 E92 受体（主链缺失 GLN61；非可迁移结构规则） | 0.621897 | 0.582812 | 0.405136 |

证据文件：`runs/fmo87/ablation/fxa_receptor_hydrogen_diagnostic_20261015/water_hydrogens/blind_predictions.json` 与相应 `tyr99_hydrogen` 臂；全受体原始坐标交叉对照 `scripts/diagnose_fxa_hydrogen_geometry.py`。**结论**：口袋水氢取向解释了此前受体差异中大部分 *FXa 分数差距*，不是“删完整残基提高指标”；0.619 的条件是回顾实验，并非保证外推泛化。已实现的 `blind-run` 原水氢候选版本使用**不读取实验标签的规则**：原始水必须能通过完整 O/H1/H2、键长/夹角、同 ID 水氧坐标不移动的检查，才在 PDB2PQR 重建蛋白氢和 PROPKA 质子化后保留该水的输入取向；否则逐水记录回退到 PDB2PQR，若水氧位置漂移则 fail-closed。分别写入 `.water_hydrogens.json` 的源/产物 SHA。**重要反例：该候选协议 Gal 冷启动 26/26 完成，但 Pearson/Spearman/Kendall 为 0.37037/0.32877/0.24459，明显低于旧 PDB2PQR 水氢协议的 0.56093/0.54567/0.38700；不能发布其为统一最佳默认协议。**原87、OpenBind、Gal、FXa 必须按同版本完整验收；不能把旧指标移植给新版，也不能利用标签为 Gal/FXa 分别挑水氢处理。FXa/GLN61 救援依旧通过原子片段诊断自动执行。新冷启动目录 `runs/fmo87/unified_fxa_source_water_default_20261015` 与 `runs/fmo87/unified_galectin3_source_water_default_20261015`；计算未完前尚无最终四集统一指标。

完整逐例指纹、7 组全量配对和 27 个 deck 全文对照保存于 `runs/fmo87/ablation/fxa_hydrogen_full_ablation_20261015.json`；当全部恢复源水及 TYR99/TYR225 氢时，完整主链版本的 27 个 QM 输入与历史受体+新姿态版本**逐行完全一致**，说明历史受体的物理缺主链段并没有进入固定片段 QM deck；其余差距是制备氢坐标差异。无标签固定重原子氢最小化也有可比的 0.620/0.610/0.428，但它是**另一个候选协议**，不可因为 FXa 回顾分数好就事后按标签对不同数据集选择处理 arm；还需四集统一重跑及未知系列外部验证。已另以同一默认入口启动 `runs/fmo87/unified_source_water_87_20261015` 和 `runs/openbind/unified_source_water_20261015`；由于之前的数据用了不同水氢协议，历史 87/OpenBind 的分数不代替本轮结论。

入口：`scripts/fmo_workflow.py`；整批回退：`scripts/rescue_fmo_target_cohort.py`；单晶体回退：`scripts/rescue_fmo_new_target.py`；测试：`tests/test_rescue_fmo_target_cohort.py`。
