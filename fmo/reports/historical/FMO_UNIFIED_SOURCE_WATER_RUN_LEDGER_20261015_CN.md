# 原水氢候选臂：全量运行与验收账本

## 版本与计算边界

单入口 `scripts/fmo_workflow.py blind-run` 的**默认**受体水氢来自 PDB2PQR。仅显式指定 `--preserve-input-water-hydrogens` 的实验候选臂，才会对晶体 PDB 中**已有**的水 O/H1/H2 检查完整原子身份、水 O 不移动、O–H 和 H–H 几何；若都合格，就以原始水氢坐标替换 PDB2PQR 的重置水氢。无完整水氢的单个水保留 PDB2PQR 输出并记录；氧不匹配时 fail-closed。不能为不含氢的实验水假设方向，也不能事后按实验结果挑方向。每个实验臂输入的水氢报告与源 SHA 保存在 prepared PDB 同名 `.water_hydrogens.json`。

默认口袋仍为该 target 全部配体共用的 5 Å 并集；严格 SCF 故障后的整队列选择性失败片段排除逻辑不变。新方案不是简单恢复旧 E92 受体：原始完整主链和 A/B 二硫键必须保留，必要时只在 QM 片段口袋上做同队列诊断修复。方法效能来自回顾性开发，需要未来尚未见标签的独立外部靶点验证。

## 默认处理与可迁移改进的准入条件（2026-10-16）

**结构审计的“去水”不是“只删除冲突水氢”。** `decide_fmo_structure_treatment.py` 对全部配体逐一寻找严重水-配体原子冲突，取整靶点冲突水**残基**的并集，仅移除命中的完整 HOH/WAT/H2O 残基，给该靶点所有配体应用同一处理；其余水保留。另一独立步骤 `reset_prepared_hydrogens()` 在 PDBFixer 后移除**全部**现有氢（不移除重原子），由 PDB2PQR/PROPKA 统一重建，避免混用 PDBFixer 氢与 PDB2PQR 质子化或 HXT。可选原水氢臂在此后仅替换已通过几何和身份审计的水氢；**默认不启用**。

默认 `k=10 kJ mol⁻¹ Å⁻²` 是对配体各原子原位谐振约束的 OpenMM 最小化，受体和水均固定；仅登记的单例救援 `--repair-clashes` 才让受体/水一起活动，不能称默认均为自由全复合物最小化。OpenBind `A71EV2A-x5836a` 两次具有相同受体输入的自由 k10+DIIS 救援，最小化终态势能分别约 −17358.65 和 −17311.93 kJ/mol：前者 SCF 收敛，后者不收敛。固定受体 k10 两次几何差异很小；自由 k10 对局部最小值更敏感。只设置 OpenMM 单线程再次运行仍未使该实验臂收敛；**不能声称固定线程能保证同一局部极小值或 639/639**。原 OpenBind PDB2PQR 默认臂已重新独立验收其历史完整 490/490 化合物结果，Pearson/Spearman/Kendall 为 0.37911/0.44374/0.29048；原水氢臂 638/639 不发表全量排序。

可复用候选只能依据输入结构和无标签物理诊断决定：首先锁定受体链、质子化/水处理、最小化平台/线程/精度和随机种子，并记录输入和终态指纹；按预注册的失败类别触发多起始几何或 SCF 求解器救援，逐一验证严重碰撞、配体漂移、QM 收敛和全靶点处理一致性；对多个合格极小值用**事先定义的无亲和力物理准则**确定唯一选择或拒绝排序。不得把“哪个分数与该靶点实验相关更高”作为水氢或最小值选择依据。候选冻结后重新跑原87/Gal/FXa/OpenBind 四队列，并在未来未参与设计的靶点上检验；不能保证提高所有靶点。目前原水氢臂在 FXa 冷启动 27/27 Pearson 0.61976，但 Gal 26/26 Pearson 0.37037（默认 0.56093），所以不能直接提升为主流程。

## 原87 / Galectin3 / FXa / OpenBind 四队列

| 队列 | 本次统一版本计算入口 | 已知任务总量 | 复核脚本／验收输出 |
|---|---|---:|---|
| 原87 | `scripts/run_fmo87_unified_water_pipeline.py` 调用四次同一 `blind-run` | CDK2 16、JNK1 21、P38 34、TYK2 16 | `scripts/evaluate_fmo87_new_full_pipeline.py --runs runs/fmo87/unified_source_water_87_20261015` |
| Galectin3 | `scripts/fmo_workflow.py blind-run`，原始 26 个晶体坐标，脱实验属性 | 26 | `scripts/evaluate_galectin3_unified_blind.py --runs runs/fmo87/unified_galectin3_source_water_default_20261015` |
| FXa | 同一 `blind-run`，原始受体的共价 A+B 双链，脱实验属性 | 27 | `scripts/evaluate_fxa_numbered_main_pipeline.py --runs runs/fmo87/unified_fxa_source_water_default_20261015`；如尚未选中有审核日志的救援则停止 |
| OpenBind | `scripts/run_openbind_new_full_pipeline.py` 的每个晶体调用同一 `intake→geometry→adaptive-run --preserve-input-water-hydrogens` | 639 晶体，490 化合物 | `scripts/evaluate_openbind_new_full_pipeline.py --runs runs/openbind/unified_source_water_20261015 --allow-rescue` |

四个任务均已启动；所有输出在 `runs/`，受 `.gitignore` 屏蔽，**不会混入仓库源码**。任意“完成”须同时核对冻结源 CSV/SDF/PDB SHA、结构预审、每个 PDB2PQR 水氢 audit、协议 SHA、全队列 GAMESS 输入/日志、严谨的盲态覆盖、固定 `k=10`、救援优先级与正常退出码；只有最后才能解盲聚合 Pearson/Spearman/Kendall。原 87 的统计口径为四靶点内相关按样本数加权，OpenBind 为多个晶体的化合物级聚合；Gal 与 FXa 为单一靶点。

旧版本指标存于 `README.md` 的明确“prior PDB2PQR water-H protocol”表，不得代替上述新版本。FXa 27/27 已完成的**前置单因素对照**独立保存在 `runs/fmo87/ablation/fxa_hydrogen_full_ablation_20261015.json`；它不证明默认版 27/27 的冷启动或另外三个数据集的迁移表现。

## 冷启动运行记录与工程检查

- FXa：`runs/fmo87/unified_fxa_source_water_default_20261015`。188/188 个原始水氢通过完整性审计；严格臂 25 例 SCF 失收敛、2 例超时，依据独立失收敛日志一致定位 QM 口袋外围 `A:GLN:61`，整靶点 27/27 选择性修复冷启动成功。已验收 Pearson/Spearman/Kendall **0.61976/0.58495/0.40514**，见 `complete_27_adaptive_tpie.json`。既不删除真实受体主链，也不借助实验标签选择失败臂。
- Gal：`runs/fmo87/unified_galectin3_source_water_default_20261015`。159/159 个原水氢保留；26/26 严格 QM 完成，验收 Pearson/Spearman/Kendall = **0.37037/0.32877/0.24459**，见 `complete_26_default_tpie.json`。旧的 PDB2PQR 水氢版本是 0.56093/0.54567/0.38700：**此水氢候选规则跨靶点退化，不得以 FXa 回顾改善为由发布为通用最佳默认流程**。
- 原87：`runs/fmo87/unified_source_water_87_20261015`。四靶点由同一冻结 87 例无标签输入快照并行，**87/87 成功并通过输入/结构/水氢/k10/预测身份验收**；四靶点按样本量加权 Pearson/Spearman/Kendall = **0.63754/0.64629/0.48408**，见 `complete_87_tpie.json`。先前 PDB2PQR 水氢协议 0.6352/0.6272/0.4635；此候选改善原87，但不能抵销 Gal 退化。
- OpenBind：`runs/openbind/unified_source_water_20261015`。冻结 639 例且每例对应独立晶体受体，逐例生成无标签 CSV；**638/639 已完成，`A71EV2A-x5836a` 在原六个登记 SCF 救援臂仍失收敛；追加按口袋重原子诊断的 DIIS 救援亦未成功，完整统计不能发布**。旧 PDB2PQR 水氢流程中该例经自由 k10+DIIS 成功，但重跑时即使蛋白准备 PDB SHA 完全一致，重最小化的配体坐标不同，SCF 再现不稳定；不得拷贝旧值冒充新协议。
- 前期 Galectin3 26/26 的 `0.5609/0.5457/0.3870` 与前期原87及 OpenBind 的旧值为**不同水氢协议**，不可称为本次同协议已验收结果；FXa 仅源水氢的 0.61888 是固定姿态回顾实验，也不是此次冷启动已验收值。
- 初步结构检查：Gal 的 159 水与 FXa 的 188 水在两种受体处理下均保持氧原子，源水氢 O–H 几何均合法。Gal 26 个 tPIE 全部改变，改变量绝对值中位数约 10.48 kcal/mol。原 PDB 中带氢**不等于实验直接测定氢方向**；仅按原水氢存在、几何合格、粗碰撞或受体固定重原子的全局 AMBER 静态能量，不能保证同时选择 Gal 原 PDB2PQR 水氢与 FXa 源水氢。不允许依据目标实验指标事后给 Gal/FXa 指定不同默认氢规则；必须另提独立、物理预注册的规则，统一复跑所有队列与外部靶点。当前 `blind-run` 继续以 PDB2PQR 水氢为默认，已完成的来源水氢计算保留为 `--preserve-input-water-hydrogens` 显式实验臂。
- 直接对**两个已算臂的 tPIE 等权取均值**（不做拟合或选靶点）也不能消除反例：Gal 26/26 的 Pearson/Spearman/Kendall 为 0.46721/0.45330/0.33128，FXa 回顾固定姿态 27/27 为 0.60190/0.56694/0.39372；两者都比各自先前的单个较好臂低。不得将多构象平均称为已完成的跨靶点修复。
- 固定全受体重原子仅优化氢的物理候选也完成 Gal 26/26（`runs/fmo87/ablation/galectin3_hydrogen_full_ablation_20261015.json`）：Pearson/Spearman/Kendall = **0.53516/0.47520/0.34366**，低于旧 PDB2PQR 原臂的 0.56093/0.54567/0.38700。FXa 固定氢优化回顾实验虽约为 0.620/0.610/0.428，但不能保证其它靶点迁移，因此不作为默认。
- 软件已通过相关定点 pytest；构建轮子 `python -m pip wheel --no-deps --no-build-isolation -w /tmp/fmo87-package-audit .` 成功；总代码测试与四数据集真实跑完后方能更新为最终验收。

发布前仍须遵守 `SUBMISSION_CHECKLIST.md` 的 GAMESS/Boltz 授权和 `vendor/`/`software/`/`runs/` 排除规范；当前文件夹没有 `.git`，也没有代用户执行 git commit。
