# FXa 统一默认流水线：失败驱动水氢救援

## 不引入分子量、亲和力标签或靶点名称分支

统一入口仍是 `scripts/fmo_workflow.py blind-run`：原受体与配体审计 → PDBFixer → PDB2PQR/PROPKA → 全靶点冲突水审计 → 固定受体 k10 → FMO2-DFTB3/PCM tPIE。正常收敛时直接使用 PDB2PQR 氢，不运行水氢救援。只有整靶点严格 QM 失败，且至少两个独立配体的 SCF 日志定位**同一个**非水、重原子距口袋截断之外、仅由全原子接触选中的外围片段时，才考虑全靶点候选救援；不依据实验排序或分子量挑选。

首先核实原 PDB 中每个水都有完整、可解析的 O/H1/H2，且水身份与 PDB2PQR 产物一一对应；把原水的氢朝向转移到**同一 PDB2PQR 准备好的受体**，并对每个冻结 k10 配体姿态做严重氢碰撞检查。任一水缺失、身份不符、几何不安全、哈希不一致即拒绝此候选，继续原有外围片段或重原子口袋救援。救援后的相同受体/相同口袋修复应用于**整靶点全部配体**，必须全部 QM 收敛后才允许输出一个排序；旧严格臂与救援臂绝不逐配体混拼。水氢转移仅修改水氢，不改变受体重原子、配体姿态、能量模型或 tPIE 定义。

实现及审计：`scripts/rescue_fmo_target_cohort.py`、`scripts/evaluate_fxa_numbered_main_pipeline.py`、`tests/test_rescue_fmo_target_cohort.py`。`adaptive_status.json`、`source_water_prepared.water_hydrogens.json`、`source_water_post_k10_geometry.json` 与选中臂的 `frozen_protocol.json` 保存无标签触发原因和输入指纹。

## 已完成的比较与冷启动验收条件

**已经通过同一冻结姿态验收的整靶点重算**：在 `runs/fmo87/unified_disulfide_fxa_default_retest_20261015` 冻结的原严格失败输入上，按注册的 `validated_source_water_failed_fragment_cohort` 重新跑全 27 例；经哈希、水氢、姿态和 QM 27/27 门控后，`fragment_repair_status.json` 指向完整盲预测；运行 `scripts/evaluate_fxa_numbered_main_pipeline.py --runs runs/fmo87/unified_disulfide_fxa_default_retest_20261015 --adaptive-report fragment_repair_status.json` 得到 Pearson/Spearman/Kendall **0.61857/0.58281/0.40514**（结果：`complete_27_fragment_tpie.json`）。相对同冻结输入旧水氢版本 **0.58267/0.56694/0.39372**，三项均改善，Pearson 提高 **0.03590**，按配体配对 bootstrap（10,000 次，种子 20261016）的经验 95% 区间 **[+0.00337,+0.07704]**。统计不确定性估计基于开发集，不能当作跨靶点泛化保证。

独立完成的历史 FXa 候选试验（各 27/27）均用完整 A+B 受体：PDB2PQR 水氢为 Pearson/Spearman/Kendall **0.58267/0.56694/0.39372**；原结构水氢为 **0.61976/0.58495/0.40514**。**历史两次冷启动的 k10 冻结协议哈希不同**，所以这不是严格的单因素固定姿态比较；旧严格臂冻结姿态上的新救援重算另行验收。按 27 配体配对自助采样（10,000 次、种子 20261016），历史两臂 Pearson 增量点估计 **+0.03709**，经验 95% 区间 **[+0.00424,+0.07478]**。这些是同一开发集上的不确定性描述，**不是独立新靶点泛化保证，也不能将增益完全归因于单独一个氢处理步骤**。历史两个运行是不同原水处理臂；新的自动触发逻辑不能冒充它们已经验收。

第一轮完整冷启动 `runs/fmo87/unified_fxa_auto_water_rescue_20261016` 的严格失败及救援 QM 已完成 27/27，但后台进程在序列化修复前启动，生成 `adaptive_status.json` 时**缺失冻结接触审计文件**，因此尚不能使用严格验收器认可这个冷启动。修复序列化之后另起进程、重新跑完整默认入口的冷启动目录是 `runs/fmo87/unified_fxa_auto_water_rescue_final_20261016`：必须同时满足严格臂确诊、原水 188/188 验证、固定 27 姿态全部通过接触审计、救援 QM 27/27、同一计分程序逐例校验，最后由 `scripts/evaluate_fxa_numbered_main_pipeline.py --runs runs/fmo87/unified_fxa_auto_water_rescue_final_20261016` 计算指标。如果没有 27/27 或缺任何审计，不报告部分样本相关系数或宣称最终冷启动通过。

**最终冷启动已通过**：`runs/fmo87/unified_fxa_auto_water_rescue_final_20261016` 用同一 `blind-run` 默认入口从原 PDB/SDF 重新开始；严格 QM 在 5 例失收敛后终止、按日志一致确诊外围 `A:GLN:61`；188/188 原水完整、水氢对 27/27 固定 k10 姿态零严重碰撞；整个救援臂 GAMESS **27/27 成功**，冻结协议/输入/预测哈希及身份均通过评估器。`complete_27_adaptive_tpie.json`：Pearson/Spearman/Kendall **0.618578/0.582812/0.405136**。这个值和独立同姿态重算（0.618572/0.582812/0.405136）吻合。严格失败后的 5 个病例不能作为某种不完整的 PDB2PQR 新冷启动评分，比较基准是独立完成的上个 PDB2PQR 默认整队列结果（0.58267）。

在旧版本相同默认入口成功的原 87（87/87、加权 **0.6352/0.6272/0.4635**）与 Gal（26/26、**0.5609/0.5457/0.3870**）不触发整靶点失败救援，应该保持原有排序；OpenBind 既没有源晶体水，且按单例目标救援，不能进入这一整靶点水氢分支，旧默认完成 639/639 晶体、490/490 化合物，**0.3791/0.4437/0.2905**。后续应重算哈希与逐例分数核验不变，不能仅因为入口代码设计为条件触发就声称这些队列已经用新版本冷启动重跑。

**限制**：这个无标签触发条件是在查看 FXa 和 Gal 的开发集差异后确定的；不是盲选规则的跨新靶点验证，也无法保证任意新靶点不退化。真实未参与设计的靶点需要按相同冻结规则做完整覆盖和三指标检验。
