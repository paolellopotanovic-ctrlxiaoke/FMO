# FMO 工程化流程状态与证据边界

## 目的与结论

本文件区分三件事：已有代码、已经实际运行的实验、已经证明可跨数据集部署的规则。前两者不自动证明第三者。当前得到的是一套**可审计、无逐 case cherry-picking 的 FMO 分析/运行框架**，不是已证明对所有 target 都高准确的普适亲和力预测器。

主设计原则：先用标签无关的结构和计算诊断决定“输入能否可信、需要哪些预注册对照”；再完整运行同一 target 的全部注册 case；将 `tPIE`、`ES`、`ES+SOLV`、`ES+CT+SOLV`、三态差值和其它模型作为分开的候选；最后只在冻结规则后用实验标签回顾评估。排序不佳时进入 `FMO_REVIEW`，不在该 target 上临时换分数。

## 工程步骤与成熟度

| 阶段 | 当前可执行流程 | 状态与限制 |
|---|---|---|
| 数据身份 | 冻结 cohort、ligand ID、SMILES、电荷、SDF/PDB 原子对应关系、输入哈希 | 已运行；身份或原子顺序不一致应 fail-closed |
| 蛋白准备 | PDB 清理、PDB2PQR/PROPKA、target receptor、QM/PBSA 原子覆盖审计 | 已接入；protonation 仍是可审查假设，不是自动正确性的证明 |
| 配体化学 | 显式氢、电荷、原子顺序、3D 坐标和分子图审计 | 已接入；缺氢时只补氢；不自动改重原子或猜测互变异构体 |
| 原始几何 | 配体—蛋白、水、氢的接触/碰撞审计；最小化前后分别检查 | 已接入；非水重原子严重碰撞不自动“修到通过”，转 `REVIEW` |
| 晶体水 | 默认保留；按完整 target 汇总源结构严重配体—水碰撞；删除集合对该 target 所有配体一致 | 已接入；水是否保留是 target-level receptor policy，不按 ligand 或亲和力挑选 |
| 姿态最小化 | 固定蛋白，仅配体最小化；记录重原子位移、RMSD、力场和 restraint | 已运行固定蛋白 k100/k25/k10。它们是不同 pose 假设，不等于真实结合构象已确定 |
| 多起点/relax | 第二起点、multistart、游离配体和结合态 relaxation、strain 分析 | 已有试验；未形成可部署的自动构象选择规则。不得选实验排序最好的一条轨迹 |
| QM pocket/fragment | 全 cohort union pocket、统一 receptor/caps、FACIO-like fragments、ACE/NME、HOP、fragment 电荷/电子检查 | 已接入；生产解释受当前 fragmentation 近似和质子化状态限制 |
| QM/PIEDA | FMO2-DFTB3/PCM；原生/selected PIEDA；输出总 TIE 和 ES/EX/CT/DI/SOLV | 已运行；完整 PIEDA 组件只从通过验证的 selected-subsystem 表解析。full log 中错列/巨大抵消的 EX/CT 不作为分项 |
| 三态/溶剂 | A/B/AB 独立总能量、PCM、PBSA/capped-PBSA 和分解分析 | 已做多轮消融；公式、协议或构象不同的值不能混合相加。三态差值不等同 tPIE |
| 多构象 | 预注册构象集合、固定权重、所有构象同 fingerprint、统一平均/median/Boltzmann 诊断 | 聚合器已实现并在外部审计；仍是 `FMO_REVIEW`，尚未证明跨 target 增益。无 per-ligand 取最小值 |
| 熵 | 显式 torsional/conformer proxy 分块采样、误差/收敛门、失败计数 | E61 的 96 blocks / 576 结果完整、0 环境失败，但 proxy 未收敛且不是完整结合熵；未达排名资格。不能把非收敛值补进分数 |
| 结果完整性 | 逐 case 状态、input/executable/parameter hash、PIEDA 闭合、全部 cohort 覆盖、协议 fingerprint | 已接入；缺 case、输入不一致或结构门未过时不发布完整 ranking |
| 统计审计 | Pearson/Kendall、target 分层 bootstrap、固定对照、外部 target、Boltz2 单独对照 | 已有审计；小样本指标不等于普适泛化，不能拿同一验证集反复选参 |

## 已有 target 证据

| 数据集 | 当前结构/流程状态 | 已得到的结果 | 部署解释 |
|---|---|---|---|
| 原始87（CDK2/JNK1/P38/TYK2） | wet k100/k25/k10 全 cohort 已有；P38 另完成精确 dry19 × 三个 k 的 full-energy 与 selected-PIEDA | wet k10：tPIE `0.533/0.412`、ES+SOLV `0.634/0.508`、ES+CT+SOLV `0.644/0.504`；P38 dry19 k10 替换后 87-case：tPIE `0.601/0.439`、ES+SOLV `0.595/0.463`、ES+CT+SOLV `0.605/0.457` | wet k10 / ES+SOLV 是冻结的原87主候选；P38 dry19 只作结构触发的对照支路，不因 tPIE 升高就替代全部分数 |
| Galectin3（26） | 结构修复后几何和多臂计算完成；Boltz2 独立结构模型与 FMO 姿态通过预注册几何门 | repaired median-arm tPIE 约 `0.238/0.175`，ES+SOLV 约 `-0.442/-0.408`；几何门 26/26 通过，但当前 ES+SOLV 对该测试集方向仍反向 | 结构/QM 流程 `FMO_READY`，能量准确性 `NOT_ESTABLISHED`；不据标签在 tPIE/ES+CT+SOLV 间换分数 |
| FXa（27） | 修复结构、union-pocket 诊断及多臂分析完成；Boltz2 独立结构模型通过 target-level 几何门 | E92 27/27 QM 成功；几何门 26/27 通过，8m 为 case-level outlier；能量准确性仍未建立 | 结构/QM 流程 `FMO_READY`，保留 8m 几何复核和 GLN61 结构触发排除；不能把姿态支持等同于亲和力准确 |
| TNKS2 | 有 charge-aware、entity-level 的外部试验 | 部分 entity 级指标较 raw ligand-level 好，但构象/重复实体定义影响结果 | 仅作诊断证据；未进入已冻结的统一生产覆盖矩阵 |

所有相关原始指标、样本数和逐 target 指标以 `runs/fmo87/ablation/final_analysis/` 及对应 `E*` manifest 为准。上表四舍五入值只用于状态沟通。

## P38 dry19 QM 结项记录

- 触发依据是 label-free 的源结构水氧—配体重原子距离门：7 个 P38 source poses 有 <2.0 Å 水冲突；策略在 P38 target 层执行，不逐 ligand 换 receptor。
- 同一 target union pocket 中固定删除 19 个水残基、保留其余 pocket 水；使用湿态固定蛋白 pose 的 k100/k25/k10 分别计算。
- full FMO2 energy：`102/102` 成功；selected `MODMOL=1` PIEDA：`102/102` 成功。二者 pre-MODMOL AB deck hashes 对应，所有 tPIE 的差不超过 `0.02 kcal/mol`。
- 完整数据：`runs/fmo87/ablation/E83_P38_dry19_geometry/summary.json`（full-energy）；`runs/fmo87/ablation/E86_P38_dry19_selected_pieda/summary.json`（可靠 selected PIEDA）；`runs/fmo87/ablation/E87_P38_dry19_geometry_analysis/summary.json`（87-case 配对统计）；总协议门：`runs/fmo87/ablation/final_analysis/fmo_structure_geometry_protocol_audit.json`。
- k10 对照显示：结构触发去水提升 P38 的 tPIE，但 87-case ES+SOLV 和 ES+CT+SOLV 均降低；tPIE 的 95% 分层 bootstrap 区间跨 0。因此证据支持“源水 clash 是需要处理的输入风险”，不支持“去水总会改善每种分数”。

## 数学与选择规则

1. **跨分数不做物理量直接相加。** `tPIE = ES + EX + CT + DI + SOLV` 的 PIEDA 闭合是分项关系；`ES+SOLV` 本身已包含在 tPIE 分解内，`tPIE + ES+SOLV` 是重复计入，不是新自由能公式。
2. **同一 score 的正仿射变换不改变排序。** 若 `s' = a s + b, a>0`，Pearson、Spearman/Kendall 排序方向不因单位换算或中心化而改善。robust-z 只有在与其它不同信息的 score 组合时才可能改变 ranking；组合公式必须预注册并用外部 target 检验。
3. **多臂一致性只说明数值稳定，不说明方向正确。** Galectin3 是直接反例：部分结构/分数臂排名互相稳定，selected 方向仍与亲和力反向。
4. **不收敛的熵没有可解释的确定排序贡献。** 先改善有效样本量/采样相关性并达到预设收敛门，再讨论增量；若仍不收敛，输出 uncertainty/`OMIT`，不能把截断估计当作已收敛自由能。
5. **结构参数触发只能看标签无关观测量。** 可用 clash、原子/电荷完整性、pose displacement、局部几何、结构/序列覆盖、重复实体和 QM failure/SCF 指标决定 `REPAIR`、`REVIEW` 或注册哪些对照；不能用某分数的历史 target 相关系数来决定当前 ligand 的 restraint 或 score。
6. **鲁棒尺度只用于跨量纲诊断或预注册组合。** 对一列观测值先计算 `med(x)` 和 `MAD(x)=median(|x-med(x)|)`，再用
   `z_robust=(x-med(x))/(1.4826*MAD(x))`；若 `MAD=0` 才使用预先登记的标准差回退。它不能把 kcal/mol 变成新的物理能量，也不能单独改变排序；多臂物理能量优先在原单位下做 median/trimmed aggregation，robust-z 只作为平行诊断或明确冻结的无量纲组合。
7. **多构象/多臂聚合必须先过一致性门。** 对同一 ligand 的注册臂只能在原子身份、受体版本、FMO fingerprint 和收敛状态一致时聚合；默认报告 median，并同时报告 arm spread。若预注册的 pairwise Kendall 一致性低于 `0.80`、覆盖不完整或任一臂失败，则输出 `REVIEW`，不按单 ligand 选最佳臂。
8. **熵的入分数门是联合门。** 必须同时满足完整 bound/free 覆盖、有效样本量达标、重复或分半估计稳定、split-half Jensen–Shannon divergence `<=0.15`，否则熵只能写入 uncertainty/diagnostic，`ranking_eligible=false`。这一区分“算出了一个数”和“这个数足以用于排序”。

## 还不能声称已完成的工作

- k10 对原87有优势，但尚未证明对所有外部 target 普适；Galectin3/FXa 的姿态方向验证仍是关键缺口。
- 多构象聚合器解决了 cherry-picking 的软件风险，不等于采样物理正确、采样数够或能提高亲和力排序。
- E61 不是完整 binding entropy；按用户先前要求，96-block nonconverged audit 保留，120/128 扩展不作为当前结果。
- 当前 `FMO_READY` 覆盖已冻结的原87协议；Galectin3/FXa 已有完整外部 QM 臂，并通过独立 Boltz2 几何门。当前报告将二者标为结构/QM `FMO_READY`，但 `accuracy_status=NOT_ESTABLISHED`；Boltz2 仍独立报告，绝不与 FMO kcal/mol 混加。
- AnewFEP 已确认对应 ChemRxiv v2（DOI `10.26434/chemrxiv.15002526.v2`），并已按本机全文更新迁移审计；它是 REST2/显式溶剂 RBFE，不是 FMO 端点法。迁移边界和 HIF-2α、BACE1、FXa、JNK1 案例逻辑见 `docs/ANEWFE_VALIDATION_AUDIT_CN.md`。

## 下一阶段顺序

1. 用原始结构、最小化后结构和预测/晶体姿态做对齐与 orientation-independent pose 检查，优先解除 Galectin3/FXa 的标签无关姿态门。
2. 将 relaxation、多起点和 conformer ensemble 统一成预注册设计：固定生成数、重复 seed、固定权重、重复性/有效样本量、所有构象都评分。
3. 对 E61 式 torsional entropy 先诊断相关时间、block variance 和有效样本量，再设计成本受限的收敛扩展；先做 pilot，不对非收敛估计排序。
4. 在至少一个完全独立 target 批次锁定结构触发/score 公式之后再评估外部迁移；旧数据作为 development evidence，不能同时充当无偏测试。

## 2026-10-04 统一决策层更新

- 已核对 AnewFEP v2 正文，并将其局部参数、蛋白构象、rotamer/采样和数值稳定性案例映射到 FMO 的可验证对照臂；见 `docs/ANEWFE_VALIDATION_AUDIT_CN.md`。
- 已实现 `scripts/build_fmo_universal_decision_report.py`，其决策不读取实验亲和力标签，输出 `FMO_READY`、`FMO_REVIEW` 或 `FMO_FALLBACK` 及具体补救动作。当前机器报告为 `runs/fmo87/ablation/final_analysis/fmo_universal_decision_report.json`。
- 已为决策器增加测试并通过；该层明确将未验证姿态、未完成标准受体臂和未收敛熵 fail-closed。
- FXa 标准 `5 Å union pocket + k=10 + FMO2-DFTB3/PCM` 的 E88 臂已经证明在原始标准 pocket 上出现统一 SCF 失败；旧 E71 reference-pocket 结果继续只作诊断。

### FXa 标准 union-pocket QM 数值诊断（2026-10-04）

- E88 使用标准 `5 Å union pocket + k=10 + PCM + MAXIT=200`，前 8 个 QM 例中均未得到成功结果；失败类型为统一的 `FAILED_SCF_CONVERGENCE`，日志显示部分 FMO 单体/二体超过 200 次迭代。该臂已停止并保留在 `runs/fmo87/ablation/E88_fxa_standard_union_k10/`，不能用于排名。
- 试验性 E89 将 `MAXIT` 改为 400，但 GAMESS 原生输入限制 `MAXIT` 最大为 200，故 27/27 在输入检查阶段 `FAILED_GAMESS_RUNTIME`；E89 不是物理失败证据，不能把它当作有效 QM 臂。
- E90 固定同一结构、受体、PCM、`k=10` 和 `MAXIT=200`，只采用历史 FXa 数值诊断中已有的 target-level `DIIS + GAMESS SOSCF` 组合；该臂最终仍以同类 SCF 超迭代失败，输出目录为 `runs/fmo87/ablation/E90_fxa_standard_union_k10_diis_soscf/`，不能用于排名。
- E91 的无标签单例诊断显示，排除新增外围 `GLN61` 后同一 FXa ligand `5f` 从 SCF 不可收敛变为成功；E92 将该规则 target-level 固化，27/27 QM 全部成功。完整结果为 `runs/fmo87/ablation/E92_fxa_union_exclude_gln61_k10/`，结构规则 manifest 为该目录下的 `structure_treatment_manifest.json`。
- E92 四个分支指标为：`tPIE=0.622/0.405`、`ES=0.652/0.399`、`ES+SOLV=0.643/0.371`、`ES+CT+SOLV=0.621/0.371`（Pearson/Kendall）。与 E71 旧 reference-pocket 的配对比较只作 protocol 消融，不据标签切换生产分数。

### 独立 HF 三态 QM 队列结项（2026-10-04）

- 对原始 87 cases 运行了独立 gas-phase `FMO2-HF/STO-3G` 的 complex/ligand 队列；覆盖审计保存在 `runs/fmo87/ablation/final_analysis/hf_independent_coverage_audit.json`，不读取亲和力标签。
- 87 个 ligand 单体计算成功，但 87 个 complex 计算均在 `3600 s` 内未完成，因而完整三态差值覆盖为 `0/87`；这不是 SCF 收敛后数值异常，而是当前大体系原生 HF/FMO2 计算成本不可部署。
- 该队列已停止并保留输入、日志和超时 `result.json`，不得用失败臂生成排序。当前生产 QM 仍采用已闭合的 `FMO2-DFTB3/PCM` selected-PIEDA 主线；HF 队列仅作为精度/可行性验证证据。
- 这一结论不等于“HF 理论错误”：它只说明在当前 fragment、几何和硬件条件下，不能在可接受成本内完成 87 个复合物的独立三态验证。

### 熵当前状态（2026-10-04）

- 当前 E61 block audit 已有 `144/144` 结果、`0` 个环境失败，但 6 条 bound/free 轨迹中仍有至少一条 split-half JS 超过 `0.15`，故 torsional proxy 仍未收敛。
- 具体审计为 `runs/fmo87/ablation/E61_chunked_CDK2_17/block_entropy_audit.json`；bundle gate 的最新输出为 `runs/fmo87/ablation/final_analysis/adaptive_bundle_validation_latest.json`。
- 在满足收敛门以前不把 `-TΔS` 加入 tPIE、ES 或 `ES+SOLV`，也不生成熵修正排序；后续先做相关时间/分块方差诊断，再决定是否需要增加独立 blocks。
- 新增 `scripts/diagnose_fmo_entropy_convergence.py`，其当前诊断写入 `runs/fmo87/ablation/final_analysis/e61_entropy_convergence_diagnostic.json`：ESS 门全部通过，唯一失败是 bound 扭转分布的 split-half JS（最大 `0.15335`，阈值 `0.15`），因此结论是“采样分布仍有轻微漂移”，不是有效样本数不足。
- 24→48→72→96 blocks 的最大 JS 总体从 `0.3555` 降至 `0.1534`，说明增加独立 blocks 有效但尚未达到预注册门；后续若继续，只追加预注册独立 blocks 并重新审计，不能用当前 `1.3169 kcal/mol` 的 proxy 增量做排序。

### 外部靶点姿态证据审计（2026-10-04）

- 新增 `scripts/audit_fmo_external_pose_evidence.py`，报告为 `runs/fmo87/ablation/final_analysis/external_pose_evidence_audit.json`。
- Galectin3（E68，26/26）和 FXa（E92，27/27）均通过来源结构预检、化学一致性检查和固定 arm 数值稳定性门；这证明输入和数值流程可复现。
- 新增 `scripts/compare_fmo_boltz_pose_geometry.py` 和 `scripts/decide_fmo_external_pose_orientation.py`，用独立 Boltz2 结构只作几何参照，不使用亲和力标签选分数。门限为蛋白 RMSD ≤2.0 Å、配体质心偏移 ≤1.5 Å、接触残基 Jaccard ≥0.8、FMO 接触重叠 ≥0.8，且 target 通过比例 ≥95%。
- Galectin3 26/26 通过；FXa 26/27 通过，唯一失败为 8m。两者因此获得 `SUPPORTED_BY_INDEPENDENT_STRUCTURE_MODEL`，解除“姿态方向未验证”这一结构门，但这不等于能量排序准确。
- 决策层当前对两者输出 `FMO_READY` 仅表示结构处理、QM 覆盖和可复现部署门通过；profile 中明确 `accuracy_status=NOT_ESTABLISHED`。Galectin3 当前 ES+SOLV 仍为反向结果，不能通过切换 tPIE、ES 或组合项来掩盖该问题。

### 统一自适应协议编译（2026-10-04）

- 新增 `scripts/compile_fmo_universal_adaptive_workflow.py`，输出 `runs/fmo87/ablation/final_analysis/fmo_universal_adaptive_workflow.json`。该文件把每个 target 的结构处理、水策略、几何臂、QM 覆盖、姿态门、PIEDA/tPIE 分支、PCM/PBSA、熵和多臂聚合规则编译成单一可执行协议。
- 规则是标签无关且 fail-closed：P38 的水处理是 target-level 统一去水并保留 dry diagnostic；所有 target 默认冻结 `k=10 + ES+SOLV`，但 tPIE、ES、ES+CT+SOLV 作为独立审计分支，禁止相加；robust-z 只作 within-target rank diagnostic，不能当作 kcal/mol。
- 熵只有在 bound/free 全部收敛且 full binding entropy gate 通过时才可进入分数。当前 E61 的 96-block round 仍因 bound split-half JS `0.15335 > 0.15` 未通过，因此生产状态仍为 `OMIT`。
- `validate_fmo_adaptive_bundle.py` 已将该统一协议列为一致性门；当前 bundle validation 通过，且明确“通过不等于亲和力准确性证明”。

### 独立 HF 三态 QM（2026-10-04）

- 原 87 例的 HF/STO-3G 全量复合物任务在 `3600 s` 内全部超时，因此 HF 不作为生产主线。
- 已用同一完整协议启动 `CDK2/17` 长时限代表例：受体已成功收敛，复合物仍在运行，输出目录为 `runs/fmo87/independent_gas_pbsa/hf_sto3g_long_cdk2_17/`。该支线用于验证 native HF 三态与 DFTB3/PCM 的物理交叉检查，不阻塞已完成的生产 DFTB3/PCM 排名。
- `CDK2/17` 三态已经完整闭合：`E(A)=-14285.04505713 Eh`、`E(B)=-3572.436611976 Eh`、`E(AB)=-17857.498726245 Eh`，得到 `E(AB)-E(A)-E(B)=-10.703516 kcal/mol`。该数值只证明 HF 三态审计器和输入协议闭合，不足以代表 87-case 排序。
- 当前 HF 覆盖快照为 `runs/fmo87/ablation/final_analysis/hf_independent_current_coverage_snapshot.json`：CDK2/17 和 CDK2/29 两个代表性案例的 A/B/AB 三态均已完成；这仍是独立物理交叉检查，不是 87-case HF 全量排名，任何未完成三态都不会进入排名或协议选择。
- 新增 `--allow-unconverged-diagnostic` 到 `scripts/score_fmo_entropy_increment.py`：当熵覆盖完整但收敛门未通过时，可以导出 provisional 增量供物理敏感性分析，但强制写入 `ranking_eligible=false`。当前 CDK2/17 报告为 `+1.3169 kcal/mol`，输出 `runs/fmo87/ablation/final_analysis/e61_cdk2_17_provisional_entropy_increment.json`；它没有进入任何生产排序。

### 因素消融结论表（2026-10-04）

正式结果写入 `runs/fmo87/ablation/final_analysis/fmo_factor_effect_summary.json`。当前原始 87-case 的回顾性证据为：

| 因素 | Pearson | Kendall | 当前结论 |
|---|---:|---:|---|
| 直接 tPIE | 0.479 | 0.320 | 保留为直接接触诊断，不与其它能量相加 |
| 固定蛋白 PIEDA-ES | 0.563 | 0.449 | 当前回顾性最强分支，但不能因此后验替换冻结生产分数 |
| 熵 proxy | 0.359 | 0.271 | 不使用；且真实 bound/free 熵尚未收敛 |
| PCM 三态诊断 | 0.428 | 0.284 | 不显示稳定增益，保持独立分支 |
| PCM+PBSA | 0.177 | 0.117 | 当前明显不稳定，不做默认组合 |
| relax 差值 | -0.050 | 0.050 | 未显示可迁移增益，保留几何诊断 |

这些是同一 87-case 上的 retrospective audit，不是无拟合外推证明；因此统一流程仍冻结 `ES+SOLV`，并将各项作为平行审计分支。

### 当前三数据集跨集泛化审计（2026-10-04）

- 新增 `scripts/audit_fmo_cross_target_generalization.py`，输出 `runs/fmo87/ablation/final_analysis/fmo_cross_target_generalization_audit.json`。它固定比较原始87、Galectin3 当前 E68、FXa 当前 E92 三个生产臂；先冻结分支，再读取标签做审计，不允许按结果回改协议。
- 当前四个原子分支中，只有单独 `tPIE` 在三组数据上的 Kendall 方向都为正：原始87 `0.412`、Galectin3 `0.174`、FXa `0.405`；其最低值仍由 Galectin3 决定，因此只能称为“方向一致的保守候选”，不能称为已证明的优秀普适模型。
- `ES`、`ES+SOLV`、`ES+CT+SOLV` 在 Galectin3 上均为反向；对四个分支做的全部 15 个等权 robust-z 子集也没有修复该反向问题。输出 `scale_aggregation_current_E68_E92_all_subsets.json` 已明确记录这一结果。
- 因而当前最稳妥的工程决策不是把所有项相加，而是：`tPIE` 作为跨集方向审计/保守主候选，`ES+SOLV` 作为原始87的冻结开发协议，Galectin3 继续标记为 `accuracy_status=NOT_ESTABLISHED`，直到找到标签无关的结构或化学根因并重新验证。
- 熵预算已由 `scripts/compile_fmo_target_entropy_plans.py` 按 target 结构复杂度编译到 `fmo_target_entropy_plans.json`；当前六个 target 均为 low tier（3 replicas × 3 blocks × 10 ns/state 起始计划），后续只由收敛门决定是否追加 block，不由亲和力标签决定。

### 目标完成度审计（2026-10-04）

- `scripts/audit_fmo_goal_completion.py` 生成 `runs/fmo87/ablation/final_analysis/fmo_goal_completion_audit.json`，逐项检查统一协议、因素消融、bundle 一致性、外部 target 覆盖、多构象/多臂、熵和独立 HF 三态。
- 当前结果为 `goal_complete=false`，不是代码故障，而是证据诚实反映：E61 完整结合熵仍未收敛；Galectin3 结构/QM 可部署但 ES+SOLV 方向反向；HF 只有两个代表性 CDK2 案例完整，不能声称 87-case HF 验证。
- 因此当前交付物是“可运行、可审计、fail-closed 的自适应 FMO 工程流程”，而不是已经证明在所有数据集上都准确的最终物理模型。后续 QM 和熵任务完成后，必须重新运行该审计，不能仅凭报告文件存在就宣称目标完成。
