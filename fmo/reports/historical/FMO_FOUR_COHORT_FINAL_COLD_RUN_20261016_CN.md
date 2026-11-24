# 四队列统一默认入口：2026-10-16 完整冷启动实测

原87、Gal、OpenBind 使用同一当前 `scripts/fmo_workflow.py` 无标签默认执行路径重新从原始输入运行；FXa 为已验收的 `blind-run` 冷启动。全部沿用 PDBFixer → PDB2PQR/PROPKA → 水冲突审计 → 固定受体 k10 → GAMESS FMO2-DFTB3/PCM tPIE；在确证失败时才走已经注册的无标签 QM 救援；未加入分子量、训练权重或逐例事后选分数。以下都是**完整样本**的结果（Pearson / Spearman / Kendall），不是部分运行的即时指标。

| 队列 | 旧默认完成值 | 本轮同版冷启动 | 完整覆盖 | 本轮验收文件 |
|---|---|---|---:|---|
| 原87：四靶点按样本量加权 | .635236 / .627219 / .463500 | **.636115 / .629477 / .464677** | 87/87 | `runs/fmo87/unified_auto_water_rescue_87_20261016/complete_87_tpie.json` |
| Galectin3 | .560934 / .545672 / .387004 | **.562159 / .542935 / .380812** | 26/26 | `runs/fmo87/unified_auto_water_rescue_galectin3_20261016/complete_26_default_tpie.json` |
| OpenBind：639 晶体，490 个化合物按晶体分数取均值 | .379106 / .443741 / .290483 | **.370362 / .443304 / .290311** | 639/639；490/490 | `runs/openbind/unified_auto_water_rescue_20261016/complete_490_adaptive_tpie.json` |
| FXa | .582673 / .566936 / .393724 | **.618578 / .582812 / .405136** | 27/27 | `runs/fmo87/unified_fxa_auto_water_rescue_final_20261016/complete_27_adaptive_tpie.json` |

**不能宣称所有指标都“零损失”**：FXa 和原87的三项均提升；Gal Pearson 略升，但 Spearman/Kendall 分别降低 .002737/.006192；OpenBind Pearson 降低 .008745，另两项降低 .000437/.000172。这里没有为了达标替换预测、删样本或按实验标签切换处理方法。运行入口：原87 `scripts/run_fmo87_unified_water_pipeline.py`（默认水氢、不传实验水氢开关）、Gal `scripts/fmo_workflow.py blind-run`、OpenBind `scripts/run_openbind_new_full_pipeline.py`（默认水氢），其验收入口分别为 `scripts/evaluate_fmo87_new_full_pipeline.py`、`scripts/evaluate_galectin3_unified_blind.py`、`scripts/evaluate_openbind_new_full_pipeline.py --allow-rescue`。FXa 独立验收见 `docs/FMO_FXA_REGISTERED_FAILURE_RESCUE_20261016_CN.md`。

## OpenBind 单例异常诊断（不得用于回溯删例）

两轮冻结的 639 个输入蛋白及配体姿态 SHA 一一一致。新运行严格完成 635、单例救援完成 4；旧运行严格完成 636、单例救援完成 3。`A71EV2A-x6838a` 的旧版固定受体 QM 收敛（tPIE −172.431），新版非常接近的固定 k10 姿态 QM 不收敛；固定姿态 DIIS、400 次迭代、charge-DIIS 和重原子口袋诊断均未修复。于是当前注册救援采用自由受体 k10 + DIIS 并成功，但其 tPIE 为 −223.225（变化 −50.794 kcal/mol）。相对新版严格姿态，救援**受体**全原子 RMSD 约 1.413 Å、794 个受体原子移位超过 1 Å，而配体 RMSD 约 0.337 Å。说明大的分数跳变与结构/SCF 救援切换相关，不能归因于新加的水氢规则：OpenBind 这些原始受体没有晶体水。

只用于定位问题、不作为正式指标：同时剔除该晶体所属第 297 号化合物后，剩余 489 个化合物旧/新 Pearson **.373911/.373751**，相差 −.000160；正式 490 例仍必须报告表中 .370362。这个单例说明“救援能成功”不等于“原物理姿态及同质分数可比”。下一轮需预注册一个**固定受体 SCF 稳定性/自由受体位移限幅**规则，在所有适用病例和未参与规则制定的数据集上重新完整验收；不能因为知道这例实验亲和力或为了恢复旧 .379 而挑选旧分数。
