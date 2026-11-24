# Pearson 成因、全队列修复实测与边界

2026-10-13。复现工具：`scripts/audit_fmo_pearson_leverage.py`；逐靶点/逐例审计：`runs/fmo87/pearson_leverage_audit_20261013/evaluation.json`。分数统一为“越高预测越强”的 `−tPIE`，实验效力也统一为越大越强；原87各靶点先独立求相关再依例数加权。所有候选都在历史实验数据已知后提出，**是回顾性分析而不是新靶点盲测**。

## 排除单纯尺度问题

下面各非原始列只用当前 target 完整队列的预测分数计算中位数/IQR/秩，**不使用实验标签**；但不同批次增删配体会改变转换结果，不能与此前冻结的 OpenBind 参考尺度混为一谈。

| 测试集 | n | 原始 Pearson | 队列 asinh | 队列百分位 | 队列正态秩 |
|---|---:|---:|---:|---:|---:|
| 原87 wet k10* | 87 | 0.533 | 0.552 | 0.536 | 0.528 |
| 原87 P38 dry19 k10* | 87 | **0.601** | 0.606 | 0.584 | 0.580 |
| OpenBind | 490 | 0.374 | 0.430 | 0.469 | 0.431 |
| Galectin3 E68 | 26 | 0.238 | 0.262 | 0.267 | 0.239 |
| Galectin3 新结构/QM | 26 | **0.556** | 0.559 | 0.561 | 0.569 |
| FXa E92 | 27 | 0.622 | 0.639 | 0.648 | 0.619 |
| FXa 已审计编号新结构/QM | 27 | 0.586 | 0.606 | 0.629 | 0.584 |
| FXa 原始 PDB→重原子口袋/QM | 27 | 0.545 | 0.574 | 0.603 | 0.538 |
| Galectin3 重原子口袋/QM | 26 | 0.421 | 0.410 | 0.391 | 0.417 |
| FXa set4 | 11 | 0.724 | 0.724 | 0.656 | 0.640 |
| TNKS2（含重复微状态） | 33 | 0.370 | 0.348 | 0.333 | 0.389 |

*wet/dry19 为同一原87数据的两臂，并非独立验证。FXa set3 仅 4 例且 Pearson −0.879，不适合选择最优臂。单调尺度变换不能纠正排序颠倒，也未对所有数据集同时提高 Pearson。OpenBind 的 MW-only Pearson 为 0.542、原 FMO 为 0.374，但**线性控制分子量后** FMO 与实验效力的偏相关仅 0.082；等人数分子量四分箱内原始相关约 −0.096、0.036、0.369、−0.143。P38 dry19 的偏相关 −0.027，而 CDK2 为 0.844；以实验标签求偏相关只是诊断，不能冒充不拟合的生产校正。OpenBind 留一最高影响例移除后 Pearson 增幅约 0.013，故不是“一个异常值害了整体”，不能删数据换取好看相关。

## 全量结构处理得到真正排序改善

Galectin3 26/26 实际从相同晶体结构输入进入新的统一流程：根据**无标签**的严重配体–水冲突对整个 target 去 HOH A403/A475/A493，筛选受体链，PDBFixer，移除它已加的氢，PDB2PQR/PROPKA 重置氢/质子化，所有配体统一 k=10 受限最小化并确认严重接触 0，固定全队列共享口袋、FMO2-DFTB3/PCM 计算 tPIE。26/26 QM 成功；输出 `/tmp/fmo-goal-gal26-20261013/qm/full_target_fmo/blind_predictions.json` 和执行状态 `/tmp/fmo-goal-gal26-20261013/qm/run_status.json`。

| 相同 26 例、仅算完后才读标签 | E68 原结果 | 新结构+QM | 改善 |
|---|---:|---:|---:|
| Pearson | 0.238 | **0.556** | +0.318 |
| Spearman | 0.245 | **0.548** | +0.303 |
| Kendall | 0.174 | **0.381** | +0.207 |

配体配对 bootstrap 的 Pearson 改善区间 `[+0.093,+0.573]`（4000 次；**不包含**方法筛选、不同新靶点、实验噪声的不确定性）；两臂原始 tPIE 的 Pearson 为 0.826，故不仅是相同排序的非线性放缩。17 系列 18 例内部 Pearson 约 0.536→0.742；18 系列 8 例约 −0.017→−0.210，改善并非每个化学系列通用。水方案、氢、受体制备、k10 与口袋均变化，不能从本次整流程对比**单独归因**某个环节；后续应从相同结构起点做逐因素配对消融。Galectin3 新流程留一审计显示剔除 17m 使 Pearson *下降*约 0.113、剔除 18b 会提高约 0.092；不允许按标签删例。

在新修复的同一 26 例上，PIEDA 的 **ES-only Pearson −0.542，SOLV-only +0.464，ES+SOLV −0.228，ES+CT+SOLV −0.095**，全部不及整体 tPIE +0.556。SOLV 有正相关不意味着应该给 tPIE **重复添加** SOLV（tPIE 已包含该分项）；多分项的方向/抵消不是统一加权公式的证据。具体逐例及相关见输出 JSON 的 `component_metrics`。

## FXa 暴露的未修复结构门缺口

使用原始 FXa 受体的 27 例过了重原子结构门、做完 27/27 k10，但 5g 的 GAMESS 出现片段 `Distance errors`，涉及内部受体残基/水，不是简单的配体–蛋白碰撞，整体安全拒绝输出 27 例 Pearson（`/tmp/fmo-goal-fxa27-20261013/qm/run_status.json`）。历史 E92 曾根据 *QM 失败而非实验效力* 在全 FXa 排除边缘 GLN61，得到 27/27、Pearson 0.622；新工具尚未自动复现这条 target 级规则。

使用 E92 事先排除 GLN61 的受体重新进入新 PDBFixer→PDB2PQR 流程，27/27 完成 k10，但 5i 再次在 QM `Distance errors` 被拒绝（`/tmp/fmo-goal-fxa27-pretreated-20261013/qm/run_status.json`）。新旧 5i deck 的片段化**不同**：旧 E92 的 GLY218 为 7 原子，新流程的 GLY218 变为 15 原子、存在近 0.07 Å 的虚拟片段间原子接触。诊断后启用旧 E92 已审计的**全 target 统一** `--repair-peptide-numbering`，并加入 FMO deck 重原子/帽 `<0.5 Å` 强重叠的 QM 前安全门；FXa 全队列 27/27 GAMESS 成功、k10 后重原子冲突检查通过、FMO PIEDA 闭合，结果 `/tmp/fmo-goal-fxa27-numbered-20261013/qm/full_target_fmo/blind_predictions.json`。

同 27 例的 E92 → 新编号修复流程：Pearson **0.622 → 0.586**（差 −0.035，配对 case-bootstrap 条件区间约 `[−0.076,−0.0001]`），Spearman **0.583 → 0.563**、Kendall **0.405 → 0.388**。这不是在 FXa 上“提升了 Pearson”，也不能拿已知 GLN61 修复过的输入冒充原始 PDB 完全自动成功。新流程的 FXa ES-only Pearson +0.647、SOLV-only −0.155、ES+SOLV +0.617；与 Galectin3 的 ES-only −0.542、SOLV-only +0.464 **方向相反**。此处分析只对历史已标注数据进行回顾，不可按哪组标签相关高就挑目标特异 ES/SOLV 权重。

**新增同规则重原子口袋消融**：原始 FXa 受体（未排除 GLN61）用全 27 配体重原子到受体重原子的 5 Å 并集选**完整残基**，加统一 `--repair-peptide-numbering` 后，27/27 在 k10 和 GAMESS 成功；原先按含氢距离纳入 GLN61 的最近接触 3.91 Å，但只看重原子最近距离 5.19 Å，重原子规则使其自然出袋。相同规则从原始 Galectin3 PDB 起点完整跑了 26/26。`scripts/audit_fmo_pocket_selection.py` 输出完整逐残基距离/直接入袋或填序列缺口的原因、源文件哈希和所有 ligand ID。FXa 450→358 个口袋原子，差 6 个残基/水；Galectin3 442→335，差 9 个（其中 LEU147、THR175 是为桥接所选残基段而入袋，并非距离 <5 Å，不能误称接触）。因此本臂解决了 **SCF/覆盖率**，不是经实证的准确率提升：

| 固定方法、全部有效案例 | 相同结构条件旧臂 | 新重原子口袋臂 | 变化 |
|---|---:|---:|---:|
| FXa Pearson（已审计 GLN61 排除输入→原始 PDB 重原子并集） | 0.622 | **0.545** | −0.077，配对 bootstrap `[−0.141,−0.007]` |
| FXa Pearson（GLN61 预排除的含氢口袋→原始 PDB 的重原子口袋） | 0.586 | **0.545** | −0.042，配对 bootstrap `[−0.073,−0.001]`；**受体起点不同，不能归因于口袋** |
| Galectin3 Pearson（同一新制备和肽编号，含氢口袋→重原子口袋） | 0.559 | **0.421** | −0.138，配对 bootstrap `[−0.336,−0.006]` |

Galectin3 另外补跑了统一肽编号、含氢口袋的 26/26 单因素对照，Pearson 0.559（原含氢臂 0.556；差 0.003），所以与重原子口袋的 −0.138 差值不再受该编号变量混淆。对应 Galectin3 Spearman 0.550→0.373；FXa 不同起点的 Spearman 0.563→0.548。以上重采样只对历史相同配体进行条件不确定性估计，**不能把经过事后提案的 Galectin3 或 FXa 当作外推测试集**。重原子口袋作为需事先冻结的、无标签可执行的 **回退/诊断实验臂** 保留，不用 Pearson 选靶点独立开关。全量运行证据 `/tmp/fmo-goal-fxa27-raw-heavy-20261013/qm/run_status.json`、`/tmp/fmo-goal-gal26-heavy-20261013/qm/run_status.json`；逐例/分项/配对数据 `runs/fmo87/pearson_leverage_audit_20261013/evaluation.json`。

FXa 补齐了同已审计受体输入的完整 27/27 重原子口袋臂（Pearson 0.548；同输入含氢 0.586），但独立重跑 k10 的配体坐标不完全相同。**更严格的单因素验证**因此另外引用含氢口袋臂的 *同一个* `poses_root` 与同一个 `prepared_protein`（源 SHA 检查，绝不重新最小化），仅改变共享口袋选取方式，分别复算 Galectin3 26/26 和 FXa 27/27；唯一变化是重原子口袋规则，结果：

| 精确同姿态、同受体、同 QM 参数配对 | 原含氢口袋 Pearson | 重原子口袋 Pearson | 差值及配对 bootstrap 条件 95% 区间 |
|---|---:|---:|---:|
| Galectin3 26 | 0.559 | 0.420 | −0.139 `[-0.339,-0.007]` |
| FXa 27（同一 GLN61 已排除受体） | 0.586 | 0.548 | −0.039 `[-0.071,-0.001]` |

复算时明确使用 `--poses-root`、`--prepared-protein`、`--pocket-heavy-only`，评估工具 `load_frozen_poses_run` 检验新旧目标、SDF、受体 SHA、ligand ID、肽编号和完整成功案例一致；FXa 口袋更改不能再被误归因于原始 PDB 中 GLN61 的删除。结果 `/tmp/fmo-goal-gal26-frozenposes-heavy-20261014/blind_predictions.json` 和 `/tmp/fmo-goal-fxa27-frozenposes-heavy-20261014/blind_predictions.json`。以上阴性结果仍是研究中提出候选后的回顾配对，不等于证明所有将来靶点都不应尝试重原子口袋。

## 固化建议

优先修复结构的真实排序问题，而不是按实验指标拟合分数尺度。生产方法必须：QM 前扫描**受体片段内部**重原子和帽接触，冻结二硫键/片段化身份，遇到几何错误 fail-closed 并给出 target 统一修复证据；所有案例采用同一受体/水/口袋和参数，记录全部输入 SHA；再在原87、OpenBind、Galectin3、FXa、TNKS2 和**尚未看过标签**的新系列/靶点上同时评价 Pearson、Spearman、Kendall、MW-only 与失败率。现有 Galectin3 +0.318 与 FXa −0.035 的并存表明**统一流水线可以固化质量与失败门，但无法从回顾数据证明 Pearson 在每个 target 都更高**。目前不能以按标签选择最优水/能量子项的方式掩盖负例。

复现：`python scripts/audit_fmo_pearson_leverage.py --galectin-run /tmp/fmo-goal-gal26-20261013/qm/full_target_fmo/blind_predictions.json --galectin-numbered-run /tmp/fmo-goal-gal26-numbered-20261013/qm/full_target_fmo/blind_predictions.json --fxa-run /tmp/fmo-goal-fxa27-numbered-20261013/qm/full_target_fmo/blind_predictions.json --galectin-heavy-run /tmp/fmo-goal-gal26-heavy-20261013/qm/full_target_fmo/blind_predictions.json --fxa-heavy-run /tmp/fmo-goal-fxa27-raw-heavy-20261013/qm/full_target_fmo/blind_predictions.json --fxa-pretreated-heavy-run /tmp/fmo-goal-fxa27-pretreated-heavy-20261013/qm/full_target_fmo/blind_predictions.json --galectin-frozen-poses-heavy-run /tmp/fmo-goal-gal26-frozenposes-heavy-20261014/blind_predictions.json --fxa-frozen-poses-heavy-run /tmp/fmo-goal-fxa27-frozenposes-heavy-20261014/blind_predictions.json`。`/tmp` 是本机非持久路径，归档需要包括各臂 QM 原始输出和结构审计，而不只是相关系数。
