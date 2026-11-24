# 当前状态（2026-10-02）

## 自适应 FMO 流水线状态

- 原 87 例：固定、无拟合的 `ES+SOLV` 主分支在统一审计下为 Pearson
  `0.633634`、Kendall `0.507750`；四种 PIEDA 分支和结构/姿态消融均已保留。
- Galectin3：E68/E69/E70 在同一删水受体、补氢 SDF 和 fingerprint 下完成
  `k=10/1/100`，26/26×3 QM 成功且 PIEDA pair 闭合。独立 Boltz2 几何门 26/26
  通过；结构/QM 部署为 `FMO_READY`，但能量方向仍未建立，生产候选保持冻结的
  `ES+SOLV`，不能后验选择系列归一化。
- FXa：当前 E92 已按 target-level GLN61 结构触发规则完成 union-pocket 27/27 QM；
  tPIE/ES/ES+SOLV/ES+CT+SOLV 为 `0.622/0.405`、`0.652/0.399`、
  `0.643/0.371`、`0.621/0.371`。独立 Boltz2 几何门 26/27 通过，8m 为
  case-level outlier；结构/QM 部署为 `FMO_READY`，但能量准确性仍未建立。
  Boltz2 仍是独立结构/亲和力对照，不与 FMO kcal/mol 分数相加。
- 姿态选择器把跨 `k=1/10/100` 的秩稳定性作为数值复现门，并额外使用独立结构
  几何门；该门只能支持结构姿态，不等于亲和力准确性。
- E61 CDK2 ligand 17：24-block/state 正式第一轮已完成但未通过收敛门，bound/free
  min ESS=`1319.4/470.0`、max split-half JS=`0.355/0.263`，点估计
  `-TΔS=+1.4021 kcal/mol` 仅为诊断。按冻结的完整 protocol 已扩展至 48
  blocks/state；截至 2026-10-02 21:08 UTC，当前正式结果为 `171/288`，环境失败 `0`，
  bound replica 1 的 blocks 46–48、bound replica 2 的 blocks 25–27 已完成，blocks 28–30 正在滚动运行。当前协议的最终熵审计尚未生成，熵项继续保持
  `OMIT/PENDING`。
- E61 轨迹审计中，原始 `ligand_com_displacement_angstrom` 是未做周期最小镜像修正的
  绝对坐标差，不能解释为结合位姿漂移/解离；截至上述检查，已完成 bound 帧的
  最近蛋白接触距离中位数约 `2.74 Å`、95 分位数约 `2.87 Å`。后续不能用该原始
  COM 列作结构质量判据；当前 block 熵报告本身只用配体扭转角，且不代表完整结合熵。
- Galectin3/FXa 的外部标签来源与符号已新增自动审计：SDF `r_exp_dg` 与独立绝对
  ΔG CSV 的最大差分别为 `4.79e-6`/`4.79e-5 kcal/mol`；Galectin3 41 条、FXa 46 条
  相对边均满足 `ddG(1→2)=G2−G1`，最大残差分别为 `7.85e-10`/`0 kcal/mol`。
  因而当前 Galectin3 的反向 ES/溶剂分支不能归因于实验标签整体反号。审计产物：
  `galectin_label_provenance_audit.json`、`fxa_label_provenance_audit.json`；工具为
  `scripts/audit_fmo_external_label_provenance.py`。
- 部署构建器已加 fail-closed 覆盖门：必须逐例成功、四分数齐全、请求/完成数量吻合，
  且所有计算臂的 ligand 身份与顺序一致；空臂/部分失败时不输出聚合分数并降级
  `FMO_FALLBACK`。Galectin3 与 FXa 当前结构/QM 报告已通过独立几何门，但
  `accuracy_status=NOT_ESTABLISHED`；总 bundle validator 已优先检查 v2 部署报告。全套测试通过，最新 bundle `passed=true`。实现/产物：
  `scripts/build_fmo_adaptive_deployment.py`、`src/sophosqm_baseline/adaptive_protocol.py`、
  `galectin_adaptive_deployment_v2.json`、`fxa_adaptive_deployment_v2.json`。
- 熵审计统计已修正：E61 的独立 seeded block 不再首尾拼接后计算 lag-1 ESS；现按每个
  block 内估计 ESS、在 replica 内求和，并比较同一 replica 前/后半组 blocks 的扭转分布。
  这避免 restart 边界制造虚假自相关/收敛；最终 E61 审计将使用该修订统计，且仍只代表
  扭转熵 proxy，不代表完整结合熵。总审计也已识别并标记旧 12-block 指标为 prior
  checkpoint，不再误报为当前 24-block 质量门结果。全套测试通过。
- 最新 target 内 robust-z 尺度聚合审计结果为：原 87 的
  `tPIE+ES+SOLV` Pearson/Kendall=`0.476/0.391`，四分支聚合=`0.500/0.399`；
  FXa E71 为 `0.605/0.411`，但 Galectin3 E68 为 `-0.068/-0.065`（四分支
  `-0.286/-0.189`）。它只能处理分数尺度，不能修复结构/状态或方向问题，因此暂不作为
  已验证生产分数。

### 独立 TNKS2 盲态结构预检（2026-10-02 15:05 UTC）

开始审计公开 TNKS2 fullmap cohort，先剥离 ligand SDF 中所有实验属性，再检查配体几何、
蛋白接触和化学状态；未读取 `tnks2_fullmap_edges.csv`。确认共有 33 个结构（不是 53 个
配体；52 条边是配体间相对自由能比较），结构文件中的两条蛋白链 A/C 都与全体配体
处于接触范围，故不能仅凭链名任意排除其中一条。结构预检器已扩展为支持多链选择，且
新增单测。A+C 合并接触门后 27/33 通过，6/33 因配体正式电荷 +1 进入 review；不是
结构缺失或几何碰撞。结果为 `runs/fmo87/ablation/E72_tnks2_blind_preflight/`
下的 `structural_cohort.csv`、去属性 SDF 和 `label_free_preflight.json`。在为带电配体
冻结一致的质子化/电荷协议前，不对 27 个中性子集单独跑后宣称完整 target 排名，也不解盲。

### E61 当前运行更新（2026-10-02 15:05 UTC）

统一 OpenMM 8.2/CUDA 12.4 的 E61 扩展目录目前有 102/144 个正式 block：bound 66/72、
free 36/72；3 个 bound replica-3 worker 仍在运行，无环境失败 artifact。结果尚未达到
预期总数，不能重算或报告最终熵门结论；此前 12-block ESS/JS 仅属旧检查点。独立 block
的 ESS/分半统计修复已通过针对性测试。

统一方法矩阵：`runs/fmo87/ablation/final_analysis/fmo_method_matrix.json`；
按 target 的 readiness 矩阵：`runs/fmo87/ablation/final_analysis/fmo_target_readiness.json`；
物理因子决策表：`runs/fmo87/ablation/final_analysis/fmo_factor_decision_table.json`；
统一 target profile：`runs/fmo87/ablation/final_analysis/fmo_target_profiles.json`；
Galectin3 化学/QM 审计：`runs/fmo87/ablation/final_analysis/galectin_chemistry_audit.json`。

## 结论

Boltz-2 论文 Table 12 的四靶点 FMO cohort 已完整重建并通过严格身份审计：

| 靶点 | 化合物数 | v4 FMO | Boltz-2 |
|---|---:|---:|---:|
| CDK2 | 16 | 16/16 成功 | 16/16 成功 |
| JNK1 | 21 | 21/21 成功 | 21/21 成功 |
| P38 | 34 | 34/34 成功 | 34/34 成功 |
| TYK2 | 16 | 16/16 成功 | 16/16 成功 |
| 合计 | 87 | 87/87 成功 | 87/87 成功 |

最终产物位于 `runs/fmo87/final/`。

## 核心指标

| 模型 | 样本加权 Pearson | 样本加权 Kendall | 论文参考 |
|---|---:|---:|---:|
| local Boltz-2 | 0.651763 | 0.509589 | Boltz-2 0.66 / 0.48 |
| 自研 raw FMO2-DFTB3 tPIE | 0.479233 | 0.320055 | 内部 FMO 0.55 / 0.38 |
| SophosQM-style in-sample fit | 0.608842 | 0.499633 | 仅次要分析 |
| SophosQM-style LOO | 0.401952 | 0.330867 | 仅次要分析 |

Local Boltz-2 与论文发表值几乎一致，支持 87 化合物队列重建和指标实现正确。自研 FMO 点估计低于论文内部 FMO，但 bootstrap 95% CI 覆盖论文值：Pearson `[0.3111, 0.6301]`，Kendall `[0.1682, 0.4306]`。

## 已完成

- 精确冻结 87 化合物 cohort 及输入哈希。
- 四靶点 PDB2PQR/PROPKA pH 7.4 制备。
- 标准 restrained flexible-complex geometry。
- shifted Cα–C′ FACIO-like 分片、ACE/NME、`HOP_C` 和电子守恒审计。
- C 末端 OXT/NME 冲突修复。
- 87 个 FMO2-DFTB3/3OB-3-1/PIEDA 全部通过 SCF、PIEDA、TIE reconciliation 和输入哈希审计。
- 87 个 Boltz-2 输入哈希、canonical SMILES、MSA、checkpoint 哈希、affinity JSON/CSV 数值审计。
- 10000 次按靶点分层 bootstrap。
- Table 12 式样本加权指标和逐靶点差异分析。

## 2026-10-02 新增的可执行部署门

新增 `sophosqm_baseline.adaptive_protocol` 和
`scripts/build_fmo_adaptive_deployment.py`，将结构质量、QM 完整性、跨 restraint 臂
数值一致性、姿态方向验证和熵门合并为一个 label-free 报告。当前重生成结果：

- Galectin3：三臂 `ES+SOLV` 稳定性审计完成，独立几何门 26/26 通过；结构/QM 为
  `FMO_READY`，但当前 ES+SOLV 方向仍反向。
- FXa：E92 union-pocket 27/27 完成，独立几何门 26/27 通过；结构/QM 为
  `FMO_READY`，8m 为几何 outlier，不能据此宣称能量准确。
- 两者均不能仅凭跨臂稳定性发布为通用 FMO 生产规则。

产物：`runs/fmo87/ablation/final_analysis/galectin_adaptive_deployment_v1.json`、
`runs/fmo87/ablation/final_analysis/fxa_adaptive_deployment_v1.json`。

### 多臂聚合回顾性审计

新增 `scripts/audit_fmo_adaptive_ensemble.py`，仅用于在部署规则冻结后报告多臂聚合的
回顾性指标，不参与 score 选择。当前结果：

- Galectin3：三臂中位数 `ES+SOLV` 为 Pearson `-0.442`、Kendall `-0.408`；稳健尺度聚合为 `-0.500/-0.396`，没有修复反向排序。
- FXa：E71 单臂无法证明几何复现性；`ES+SOLV` 为 `0.606/0.371`，因此只保留为诊断结果。

产物：`runs/fmo87/ablation/final_analysis/galectin_adaptive_ensemble_audit.json`、
`runs/fmo87/ablation/final_analysis/fxa_adaptive_ensemble_audit.json`。

### E61 运行时状态（2026-10-02 03:00 UTC）

bound 状态的 9/9 个 block 已完成；free 状态启动第一个 block 时触发
`CUDA_ERROR_UNSUPPORTED_PTX_VERSION (222)`，因此剩余 9 个 block 尚未运行。CPU 探针
验证可运行，但 1000 步约 88 秒，不能作为当前 1,000,000 步 block 的实际回退方案。
调度器已支持 `FMO87_PLATFORM`、`FMO87_DEVICE`、`FMO87_CUDA_VISIBLE_DEVICES` 和
`FMO87_CUDA_COMPILER`，并会为运行时失败写入可恢复的环境失败记录；已完成 block 不会
重跑。该问题被单独记录为环境兼容性失败，不能解释为熵统计不收敛。

### CUDA runtime 已恢复（2026-10-02 03:43 UTC）

定位到原问题为 OpenMM 8.4 wheel 与当前 550.144.03 驱动/12.8 PTX 组合不兼容。
隔离安装 `OpenMM 8.2.0 + CUDA 12.4` 后，GPU1 最小 free-state probe 成功；E61 已用
该环境续跑，当前跳过 9 个 bound block，正在执行 free replica 1/block 1。该修复不改变
采样参数、随机种子、力场或收敛门。

### E61 provenance 统一重跑（2026-10-02 04:15 UTC）

完整 block 审计发现 bound 使用旧 `OpenMM 8.4/CUDA 12.9`、free 使用新
`OpenMM 8.2/CUDA 12.4`，因此严格 provenance 门拒绝把两态合并；这是正确拒绝，
不是统计失败。旧 bound 已备份到 `legacy_openmm84_bound/`，当前正在用统一
`OpenMM 8.2/CUDA 12.4` 重跑 9 个 bound block，已完成的 9 个 free block 保留并会在
重跑后重新审计。

### E61 熵扩展与环境修复（2026-10-02 07:30 UTC）

第一轮 6 blocks/state 完成后，ESS 和 split-half JS 仍未通过质量门，已冻结扩展理由并
追加 7–12 blocks/state。扩展首批后台任务曾因未传入 `openmmforcefields` 的 Python
路径而失败；该问题已定位为调度环境传递错误，不是 CUDA 或轨迹统计失败。调度脚本已
内置 `FMO87_ENV`/`PYTHONPATH` 组装，当前使用持久在线调度器重新运行，已有 36 个
完成 block 不会重跑。历史失败 artifacts 会在成功重试时清理，并在最终审计中单独计数。

### TNKS2 带电配体整组盲算启动（2026-10-02 15:28 UTC）

外层脚本原先无理由地拒绝所有非零总电荷；底层 FMO fragmentation 已支持 ligand fragment 正式电荷并向 `ICHARG` 输出，且电子数/自旋检查适用。现已移除该阻断，保留 SDF 的显式氢与电荷，冻结电荷直方图和“跨电荷原始分数仅诊断”限制，并增加 8 路独立 ligand 任务调度。对正电 `8a` 和中性 `8a neutral` 生成的 deck 审计确认 ligand fragment `ICHARG` 分别为 `1` 和 `0`。

TNKS2 E75 使用同一 PDB2PQR pH 7.4 受体和完整 33 例（27 个 0 电荷、6 个 +1），33/33 固定受体最小化成功；全原子 clash 审计无严重碰撞、无严重冲突水，门为 `PASS_AFTER_MINIMIZATION`。QM 于 15:26 UTC 启动；截至 15:28 有 17/33 成功、6 个运行中，其余排队，尚未解盲，不能报告性能。输入/协议/结构审计在 `runs/fmo87/ablation/E75_tnks2_charge_aware_parallel/`。

同一时刻 E61 扩展为 105/144 blocks（bound 69/72、free 36/72），无环境失败 artifact；仍未完成最终熵收敛审计。外层电荷态处理和多 worker runner 的完整 pytest 套件已通过。

### TNKS2 解盲与微观态归并结果（2026-10-02 16:18 UTC）

E75 的 33/33 GAMESS FMO2-DFTB3/PCM 结果已完成、QM coverage 为 `FMO_READY`，化学与结构审计通过；预测在查看标签前以 `blind_predictions_freeze.json` 冻结，包含协议/预测/每个 deck 与 log 哈希。之后才读取 SDF `r_exp_dg` 与 52 条边表；方向 `ddG(first→second)=G(second)-G(first)` 全部吻合，最大残差 `8.41e-5 kcal/mol`（门限 0.002）。

原始 33 记录中 tPIE Pearson/Kendall=`0.387/0.172`，ES=`-0.084/-0.113`，ES+SOLV=`0.214/0.053`；target 内 robust-z `tPIE+ES+SOLV`=`0.318/0.156`，四分支=`0.209/0.034`。因此在原始记录口径上，tPIE 最好，额外项与 robust-z 均未改善。

结构去氢/中和后发现 33 条记录实际含 27 个不同共价实体，其中 6 对中性/阳离子状态具有相同实验标签；52 条原始边含 6 条同实体态间边，归并后为 43 个实体对。给定 TNKS2 `macropka` 文件声明的 pH 7.7 和 6 个 pair 的 pKa（覆盖 12 个微观态），按主态做 label-free 路由后，27 实体的 tPIE 为 `0.558/0.245`，ES+SOLV 为 `0.429/0.154`，`tPIE+ES+SOLV` robust-z 为 `0.497/0.211`；43 个实体边上 tPIE 为 `0.455/0.200`，robust-z 组合为 `0.350/0.143`。这是有提示性的外部证据，但只有 6/27 实体有成对 pKa 状态、21 个仅一个输入态，且 pH 7.7 是否为真实 assay pH 尚待来源核实；不得将 pH 主态路由视为已验证生产方法。算术布居加权的 PIEDA 分数也仅是诊断值，不是严格热力学平均。

可复现产物：`runs/fmo87/ablation/E75_tnks2_charge_aware_parallel/unblinded_edge_analysis.json`；流程为 `scripts/freeze_fmo_external_predictions.py` + `scripts/analyze_fmo_external_edges.py`。冻结器现显式标记：混合 0/+1 电荷且未建模 assay pH 结合热力学时，即便 33/33 QM 覆盖完整也只能 `FMO_REVIEW`，不把 coverage `FMO_READY` 误当成科学部署就绪。

### TNKS2 几何 restraint × pocket 消融（2026-10-02 16:18 UTC）

结构质量门按预定规则触发整组 k=1/10/100。先跑自身姿态 union pocket 时，k=1 有 483 pocket atoms，k=10/100 是 464；因此另跑了以 k=10 结构冻结 pocket 的 k=1/100 控制臂。所有臂都是 33/33 GAMESS 成功、同一 receptor PDB SHA、无严重 clash/水冲突；k=1 固定口袋与自身口袋差异只有 19 pocket atoms。全套 k=1/100 都是在看过 k=10 实验标签之后启动，所以均严格记为 `posthoc_after_prior_unblinding`，不能作为独立验证或选择最终 restraint 参数。

| restraint/pocket | pocket atoms | 实体级 tPIE r/τ | 实体级 robust-z(tPIE+ES+SOLV) r/τ |
|---|---:|---:|---:|
| k10，自身 pocket（初始盲测） | 464 | 0.558/0.245 | 0.497/0.211 |
| k1，自身 pocket | 483 | 0.537/0.231 | 0.501/0.194 |
| k1，固定 k10 pocket | 464 | 0.540/0.240 | 0.495/0.211 |
| k100，自身 pocket | 464 | 0.594/0.314 | 0.523/0.234 |
| k100，固定 k10 pocket | 464 | 0.594/0.314 | 0.523/0.234 |

固定 pocket 后 k1 与 k10 的 tPIE 排名 Kendall=`0.973`，中位绝对能量差 `0.495 kcal/mol`；自身 pocket k1/k10 排名 τ=`0.969`，说明 19 个原子的 pocket 差异在本例不是主排序扰动源。k100 vs k10 固定 pocket 的 τ=`0.890`、tPIE 中位绝对差 `3.255 kcal/mol`；k100 的相关性更高是重要几何敏感性线索，但后验结果不能作为推荐 k=100 的证据。三个 restraint 下 robust-z 组合均未超过 tPIE。完整审计 `runs/fmo87/ablation/E80_tnks2_restraint_pocket_factorial_summary.json`；每臂的 blind score、审计与解盲数据保存在 E75–E79 目录。

截至 2026-10-02 16:18 UTC，E61 为 118/144 blocks（bound 72/72、free 46/72），无环境失败；熵审计仍未完成，熵仍不能纳入总分。TNKS2 电荷/实体分析、冻结/解盲工具及 restraint-pocket 全因子新增后，全套 pytest 已通过。

### 2026-10-02 16:38 UTC 进展更新：robust 聚合与 E61

重新生成 `adaptive_bundle_validation_latest.json` 并重算 E61 部分审计。robust-z 的当前
数据集对照已确认：原始 87 例 ES+SOLV 为 `r=0.634/τ=0.508`，robust-z(tPIE+ES+SOLV)
为 `0.476/0.391`、四分支为 `0.500/0.399`；Galectin3 E68 修复集分别为
`-0.434/-0.393`、`-0.068/-0.065`、`-0.286/-0.189`；FXa E71 修复集分别为
`0.606/0.371`、`0.605/0.411`、`0.605/0.394`。所以等权 robust 尺度组合没有跨数据集
一致增益，尤其不能修复 Galectin3 的方向问题；它仍是诊断比较，不是新默认分数。

E61 正式活动目录为 `125/144` block（bound `72/72`、free `53/72`）；另有 9 个旧
OpenMM 8.4 bound 结果位于 `legacy_openmm84_bound/`，不计入本轮。无当前环境失败
artifact。部分审计没有完整覆盖且未通过 split-half JS 收敛门，故熵项继续 `OMIT`；
后台调度器和 144-block 最终审计 watcher 仍在运行。此前的 finalizer 日志有旧轮次记录，
不得把日志中的旧 `completed=72` 当作当前计数；以活动目录中符合当前 protocol 的
`result.json` 和最终 protocol-hash 审计为准。

当前方法学边界：三个集合规模/来源不同，且 Galectin3/FXa 属外部诊断、TNKS2 含混合
电荷态。上述表现不能证明标签无关的普适精度；结构质控与身份覆盖是必要门，但并未解释
或纠正所有 target-specific 方向错误。现已定位 AnewFEP 预印本 DOI 和摘要，但 ChemRxiv
全文在当前环境不可访问；摘要级迁移及 FMO 不可照搬项已记录于
`docs/ANEWFE_VALIDATION_AUDIT_CN.md`，不推断文章未公开的具体实现。

### 2026-10-02 16:45 UTC 更新

E61 活动 block 计数推进至 `127/144`（free `55/72`，bound `72/72`）；3 个 free GPU
任务正在运行，环境失败数仍为 0，正式收敛审计等待全部 144 块结束。新增恢复计划回归修复：
`scripts/plan_fmo_entropy_recovery.py` 不再把恢复动作固定写成 `block_index=2`，而是
从当前协议 `blocks` 或审计 `protocol_blocks` 后追加下一块；测试覆盖从 24 追加到 25。
这避免熵未收敛时恢复计划重复指向已运行 block、无法真正增加独立样本的问题。
相关定向测试和全套 pytest 均通过。该修复只改恢复计划生成，不改变正在执行的 E61
protocol、seed、步数、环境或收敛阈值。

### 2026-10-02 16:58 UTC 熵恢复/审计稳健性修正

E61 当前为 `131/144`（bound `72/72`，free `59/72`），0 个环境失败 artifact。用显式
部分审计工具最新重算得到 bound min-ESS `1319.4`、max split-half JS `0.355`；free
min-ESS `309.7`、max split-half JS `0.798`。free 仍缺 13 个 block，且两态现有 JS
均未满足 `0.15` 门，故不能判断最终收敛，熵继续 `OMIT`。

发现并修复两项通用流程问题：一是熵采样预算曾把模拟 `ns` 直接当作 GPU 小时比例，现
要求提供同体系/设备实测 `--calibrated-hours-per-block` 才估算 GPU 小时；未知运行时或
预算不足均显式报告，不自动减少预声明 block。二是部分审计在缺整条 replica 时会把
后续 replica 错映射到前一个编号；现保留原 replica ID，并以缺失 replica 1 的测试验证
输出 `[0,2,2]`。相关测试通过，全量 pytest 再次通过。修改不影响 E61 当前冻结协议或
正在进行的 3 个 GPU 作业。

### 2026-10-02 17:04 UTC 采样更新

E61 活动 protocol 目录已推进至 `134/144`（bound `72/72`、free `62/72`），0 个环境失败；
free block 仍在 GPU 执行中。摘要级 AnewFEP 方法映射现已识别并加入
`docs/ANEWFE_VALIDATION_AUDIT_CN.md`：将力场/化学身份、蛋白构象控制和采样/数值稳定性
映射为 FMO 的独立审计层，同时明确 REST2、λ schedule 与静态 FMO 端点并不等价。

### 2026-10-02 17:15 UTC E61 自动续采已挂接

E61 活动目录为 `137/144`，7 个 block 尚未完成，无环境失败；基础 24-block 的最终
protocol-hash 审计完成后，持久执行会话 `77291` 中的控制器将按每轮 +24 block、保持原
seed/环境/步数/阈值续跑，直至 torsional proxy 通过或遇到完整性/运行时错误。每轮协议与
审计会归档到 `entropy_rounds/blocks_NNNN/`；即便 proxy 收敛也明确不具备完整结合熵，
不自动加入 affinity score。入口 `scripts/continue_fmo87_entropy_until_converged.py`；
控制器与续采计划测试及全套 pytest 已通过。该 E61 案例仅为 CDK2 ligand 17 熵采样诊断，
不能替代 87 配体队列各自的熵审计。

### 2026-10-02 17:32 UTC E61 首轮正式审计完成

E61 首轮活动 protocol 的 `24×3×2=144` blocks 已全部通过完整性和 provenance 审计，
protocol hash 与当前 24-block protocol 一致，0 个环境失败。正式结果：bound 每个
replica `24` blocks、min ESS `1319.4`、max split-half JS `0.355`；free 每个 replica
`24` blocks、min ESS `470.0`、max JS `0.263`。因此覆盖和 ESS 达标，但两态的 JS
均未达到 `0.15`，`converged_torsional_proxy=false`；临时 `-TΔS=+1.4021 kcal/mol`
仍不可用于排序，也不是完整结合熵。

控制器已归档首轮到 `runs/fmo87/ablation/E61_chunked_CDK2_17/entropy_rounds/blocks_0024/`，
并追加至 48 blocks（首轮结果不重跑，追加独立 block 25–48）。这正是“未收敛就增加
有效样本再审计”，不是调阈值或后验调权重；当前新增轮次正在运行。

### 2026-10-02 17:34 UTC 原 87 例 Boltz2 对照刷新

`runs/fmo87/final/boltz_production_audit.json` 已重新严格审计为 `87/87`，输入身份、共享
MSA、SMILES、checkpoint hash 和协议均通过。`evaluate_fmo87.py` 使用相同冻结 cohort
重新生成 `runs/fmo87/final/REPORT.md` 和 `metrics.json`：Boltz2 `r=0.6518/τ=0.5096`，
原始未拟合 FMO tPIE `0.4792/0.3201`；分靶点 Boltz2 为 CDK2 `0.760/0.600`、JNK1
`0.797/0.625`、P38 `0.438/0.361`、TYK2 `0.808/0.583`。这只是同 cohort 外部
对照，Boltz2 输出不是 FMO 能量，也不与 FMO 相加；FMO 的 zero-fit LOO 和分支审计仍
单独保留。

### 2026-10-02 17:40 UTC 熵增量不确定度审计增强

在不改变 E61 运行中 protocol 的前提下，熵审计新增独立 block bootstrap：每个状态按
block 重采样而不是按 frame 重采样，输出每个扭转边缘熵的 95% 区间，并将 bound/free
区间以保守边缘组合形式转换为 `-TΔS` 区间。该区间只表示有限扭转 proxy 的采样不确定度，
不表示完整结合熵，也不改变 `full_binding_entropy_available=false` 和“未过 JS 门不得
进排序”的规则。测试已覆盖输出结构，全套 pytest 通过；48-block 续采仍在运行。

### 2026-10-02 17:46 UTC 活动协议验证器修正

总体验证器此前把 48-block 活动协议简单报告为“prior sampling checkpoint”。现已增加
当前 partial audit 的 protocol hash 和 expected block count 判断：活动续采期间明确报告
`active block protocol has only a partial diagnostic audit`，只有真正 hash/覆盖完整的
final audit 才报告 `entropy_audit_current=true`。当前 48-block 协议期望 `288` 个结果，
已观察 `144` 个，partial audit 当前有效、最终审计未完成；全套相关测试通过。

### 2026-10-02 17:50 UTC 方法说明同步

更新 `docs/FMO_FINAL_ADAPTIVE_METHOD_CN.md`：修正过期 E61 数字为首轮正式 24-block
审计结果和当前 48-block 续采状态；修复步骤编号；补充“适配规则决定可执行性而非
无标签推导准确性”的数学边界，包括结构/QM gate、跨 arm 一致性、robust rank-scale
聚合和 prospective 方向验证。该表述明确了为什么不能保证同一 FMO 分支在每个 target
都同样好，也避免把稳定性误认为准确性。

### 2026-10-02 17:52 UTC 48-block 续采运行状态

E61 48-block 续采仍在严格执行完整 `1,000,000` production steps/block；当前已完成首轮
`144/288`，新增 block 正在写入轨迹，3 个 GPU worker 均正常运行。不能用未完成 block
的 frame 数提前生成最终收敛结论，也不缩短 block 或改变 replica/seed；这样会破坏与
首轮及 bound/free 状态的可比性。最终 runner 完成后会用当前 48-block protocol hash
执行严格 audit、bootstrap 区间和 bundle validator。

### 2026-10-02 17:55 UTC 已有熵/relax 消融的增量证据

在已完成的 87 例 E11 两个独立 conformer ensemble 诊断中，`T·S` proxy 的两次 seed
相关性为 Pearson `0.8014`、Spearman `0.7926`，但其直接排序相关仅约
`0.086/0.041` 与 `0.014/0.008`。把 proxy 加到 fixed tPIE 的回顾性增量分别为
`Δr=+0.0032/+0.0029`、`Δτ=+0.0163/+0.0163`；加入 fixed ES+SOLV 则为
`Δr=-0.0020/-0.0074`、`Δτ=+0.0090/-0.0004`。relax penalty 明显不利：加入
fixed tPIE 后 `Δr=-0.3149`、`Δτ=-0.1112`。这支持当前规则：即使有限熵 proxy
收敛，也必须单独做增量和 seed/区间稳定性审计，不能默认加到 ES+SOLV；relax
不应作为无标签生产修正项。

### 2026-10-02 17:56 UTC 48-block 续采开始产生新结果

E61 新轮次不再停留在启动状态：活动目录已从首轮 `144` 推进到 `147/288`，bound
replica 1 的 blocks 25–27 已完整结束，轨迹均达到 200 frames/result.json。调度器仍以
3 路并行继续执行其余新增 blocks；最终审计仍等待 `288/288`，暂不提前合并首轮和新轮次
的 entropy 统计。

### 2026-10-02 17:58 UTC 续采调度更新

E61 新增 block 调度已从 bound replica 1 的 25–27 进入 28–30；活动正式结果仍为
`147/288`，没有运行时失败。每个新 block 继续使用冻结的完整 production protocol，
不使用短 probe 结果替代正式样本。

### 2026-10-02 18:02 UTC 当前 bundle 校验

全局 `validate_fmo_adaptive_bundle.py` 已用 `runs/fmo87` 根目录重新执行并通过一致性门，
但结果明确包含未完成警告：E61 当前正式 block 结果为 `147/288`，`partial_entropy_audit_current=true`，
最终 protocol-hash audit 尚未生成，且 full binding entropy 仍不可用。这里的 `passed=true`
只表示文件、协议和部署状态没有出现新的结构性矛盾，不代表熵结果已收敛，也不代表可以发布
新的 entropy-adjusted ranking。控制器和 3 路 CUDA worker 继续按冻结的 48-block 协议运行。

### 2026-10-02 18:10 UTC 回归验证

在 E61 续采仍处于 `147/288`、最终熵审计尚未生成的状态下，使用项目环境重新运行完整
pytest：`100%` 通过（无失败）。这只验证流水线代码、协议审计和 fail-closed 规则没有
回归；不改变 E61 的 `OMIT/PENDING` 状态，也不把 partial audit 解释为收敛证据。

### 2026-10-02 18:21 UTC E61 续采推进与 FXa SCF 诊断

E61 的 blocks 28–30 已分别于 18:18–18:19 UTC 完成，活动正式结果由 `147/288`
推进到 `150/288`，调度器已自动开始 bound replica 1 的 blocks 31–33；当前仍无环境失败。
FXa 标准 pocket 的 `5f` 日志进一步显示 E62 在相邻 FMO iterations 间出现明显振荡
（例如 `maxD` 在约 `1.69` 与 `0.03` 间反复），E63 DIIS 和 E65 SOSCF 也没有形成
稳定下降轨迹，最终都在 layer-1 monomer SCF gate 失败。这支持“耦合片段/电荷响应
不稳定”而非单纯增加 `MAXIT` 的诊断；不把失败日志改写成可用能量。

### 2026-10-02 19:19 UTC E61 续采状态

E61 活动目录已推进至 `156/288` 个正式 block，环境失败 artifact 为 `0`；bound replica 1
的 blocks 37–39 正在运行。当前只存在首轮 `blocks_0024` 的旧审计文件，48-block
protocol-hash 对应的最终 audit 和 `entropy_until_converged_status.json` 尚未生成，故
本轮仍不能报告收敛熵或生成 entropy-adjusted ranking。剩余 `132` 个正式 block 由冻结的
3 路 CUDA worker 继续执行。

### 2026-10-02 20:16 UTC E61 block 续采

bound replica 1 的 blocks 43–44 已完成，活动目录推进至 `164/288`；环境失败 artifact 为 `0`，
调度器已启动 blocks 46–47。最终 48-block protocol-hash audit 尚未生成，因此熵项继续保持
`OMIT/PENDING`。

### 2026-10-02 20:24 UTC E61 计数同步

E61 活动目录实际正式结果推进至 `168/288`，环境失败 artifact 仍为 `0`；验证报告已重新生成并同步该计数。
bound replica 1 的 blocks 46–48 已完成，bound replica 2 的 blocks 25–27 正在 CUDA 运行；最终 48-block protocol-hash audit 尚未生成，熵项继续保持
`OMIT/PENDING`。
