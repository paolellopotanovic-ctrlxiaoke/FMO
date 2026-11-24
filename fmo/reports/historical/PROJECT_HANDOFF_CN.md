# 项目交接与下一阶段执行文档

## 1. 当前结论

截至 **2026-09-26**，本项目已经完成独立 GAMESS FMO87 baseline：冻结的 Ross/FEP+ 四靶点队列为 CDK2/JNK1/P38/TYK2 共 87 个中性化合物，FMO 与 Boltz-2 本地结果均可按严格身份审计闭环。

当前生产 baseline 不是 Boltz-2 论文作者的私有 in-house FMO 等价复现。baseline 的 target-average 指标为：

| 模型 | Pearson | Kendall | 论文参考 |
|---|---:|---:|---:|
| Local Boltz-2 | 0.6518 | 0.5096 | 0.66 / 0.48 |
| 本项目 raw FMO tPIE | 0.4792 | 0.3201 | 0.55 / 0.38 |

因此后续重点是协议敏感性和可证伪的全量验证，而不是反复调整 87 个样本去追逐单一目标数字。

## 2. 论文与实现差异

| 维度 | SophosQM 原论文 | Boltz-2 论文内部 FMO | 本项目公开 baseline |
|---|---|---|---|
| QM 主方法 | FMO2-MP2/6-31G* | 作者 in-house DFTB3 FMO，细节未公开 | GAMESS FMO2-DFTB3/3OB-3-1 |
| 轨道/校准 | tPIE + VolSurf clogP，按 series 拟合 | 论文表中只给 raw FMO 排名指标 | raw tPIE；校准只作次要消融 |
| pose | Glide SP，核心对齐，多 pose 选择 | Glide pose | 公开 Ross/FEP+ 结构输入；不能证明与作者 Glide pose 完全相同 |
| 几何 | OPLS3e/MacroModel 约束最小化 | 未公开 | OpenMM AMBER14/OpenFF/GNN AM1-BCC |
| 片段 | 残基/水/配体独立片段，FACIO 切口与 HOP | 未公开 | shifted Cα–C'、ACE/NME、HOP_C、配体/水独立片段 |
| 水与溶剂 | 保留 5 Å 内水并优化 H-bond 网络；具体 FMO 溶剂细节需看原始工作流 | 未公开 | 结构水 + 水 PCM；敏感性中测试 gas/PCM 与 UFF |
| clogP | VolSurf | 未公开 | RDKit Crippen，仅用于 SophosQM-style 消融 |
| 软件可复用性 | Schrödinger/FACIO/VolSurf 等组件不完全公开 | in-house FMO 不公开 | GAMESS、队列、审计和编排代码可复用 |

关键判断：原论文的 SophosQM 回归思想不能直接替代 Boltz 表 12 的 raw FMO；公开材料也不足以恢复 Boltz 私有 scorer 的所有输入与数值选项。

## 3. 已实现内容

- `src/sophosqm_baseline/fmo_input.py`：`FMODeckOptions` 显式化 PCM、DFTB dispersion、MODESD、MODGAM、NACUT 等协议参数；GAMESS `$DFTB` 输入超过 80 列时提前拒绝，避免生成会触发 scanner/段错误的 deck。
- `src/sophosqm_baseline/jobs.py`：水策略和 deck options 进入 job/structure/fragment manifest，便于审计。
- `scripts/run_fmo87.py`：支持 `--no-pcm`、`--dftb-dispersion`、`--dftb-modesd`、`--dftb-modgam`、`--drop-waters`。
- `scripts/run_fmo_sensitivity.py`：支持固定 12-case 或全 87-case 协议敏感性实验，保留失败原因、GAMESS 结果和 target-weighted 指标。
- `scripts/evaluate_fmo_protocol.py`：把敏感性结果物化为标准 `fmo_features.csv`，并以相同的 target-average 口径检查完整性、成功率和排序指标。
- `scripts/analyze_fmo_components.py`：从真实 GAMESS 结果汇总 ES/EX/CT/DI/solvation，并比较可解释分量组合。

## 4. 已完成敏感性实验

### 4.1 12-case结果

`runs/fmo87/sensitivity_12/summary.json` 和 `runs/fmo87/sensitivity_12_default_modes/summary.json` 是机器可读结果。

| 协议 | 成功数 | 观察 |
|---|---:|---|
| current PCM | 10/12 | 现有生产协议可运行 |
| PCM + UFF，默认 GAMESS 分支 | 10/12 | 稳定性接近 baseline；12-case 上未显示明确排名提升 |
| gas + UFF，默认分支 | 3/12 | 多数 SCF 不收敛 |
| gas，无 dispersion，默认分支 | 3/12 | 多数 SCF 不收敛 |
| gas/PCM + `MODESD=2, MODGAM=1` | 0/12 | `lumpesd_` 段错误或 SCF 失败，不可作为科学比较结果 |

这里的失败不是“指标差”，而是 GAMESS 数值/实现稳定性失败。特别是当前 GAMESS build 在 FMO DFTB 的 `MODESD=2` 路径中出现 `lumpesd_` segmentation fault；该候选不能继续用于 87-case 追分。

### 4.2 已完成的全量实验

`runs/fmo87/pcm_uff_default_full/` 已完成 87-case 的 PCM + UFF + 默认 MODESD/MODGAM 协议。结果为 79/87 成功，失败集中在 CDK2 3 例和 P38 5 例；成功子集 target-average Pearson/Kendall 为 `0.2068/0.3324`。CDK2/P38 排序明显恶化，虽然 TYK2 提升到 `0.8892/0.7833`，但不能抵消其他靶点损失。因此 PCM+UFF 不升级为生产协议，baseline 仍保持冻结。

## 5. 指标解释

- raw tPIE 是配体-环境 pair interaction 的焓性 proxy，不是绝对结合自由能。
- 正向共同 affine 变换不会改变靶点内排序；SophosQM 的 `α·tPIE + β·clogP + γ` 需要按 series 拟合，in-sample 相关性不能当作外推能力。
- 本地 baseline 的 `tPIE+clogP` in-sample 为 0.6088/0.4996，但 LOO 为 0.4020/0.3309；这支持“校准可改善尺度、不能证明复现私有 FMO 排序”的判断。
- target-average 指标必须先按四个 target 分别计算，再按化合物数加权；不能用所有 87 个点直接全局相关替代论文口径。
- Boltz-2 的标准 YAML 输入不包含复杂 FMO 制备，但 exact reproduction 仍依赖 sequence、共享 MSA、canonical SMILES、Boltz 版本/权重、seed、recycling、affinity sampling、diffusion samples 和 affinity-head ensemble；因此本地 `0.6518/0.5096` 与论文 `0.66/0.48` 的小差异不能归因于 FMO 式结构制备。

## 6. 下一步执行顺序

1. 不把 `pcm_uff_default_full` 纳入生产 baseline；它是已完成的负向协议对照。
2. 停止继续跑 gas/`MODESD=2` 变体：前者主要 SCF 不收敛，后者在 `lumpesd_` 路径崩溃。
3. 若继续缩小差距，优先做小规模、可解释的 pose/几何/fragment-water 敏感性，而不是在固定 87 个样本上调后处理系数。
4. 任何新候选必须同时记录：deck hash、GAMESS executable hash、参数目录、结构/片段 manifest、成功率、按靶点 Pearson/Kendall、失败栈和运行时间。

## 7. 重要入口

- 论文/方法差异：`docs/METHOD_ANALYSIS_CN.md`
- 最终 baseline 报告：`docs/FINAL_REPORT_CN.md`
- 生产复现：`docs/REPRODUCTION_CN.md`
- 协议实验：`scripts/run_fmo_sensitivity.py`
- 协议评估：`scripts/evaluate_fmo_protocol.py`
- 完整分析计划：`docs/ANALYSIS_PLAN_CN.md`
- 冻结队列：`data/fmo87_cohort.csv`

## 4.3 水、cutoff 与几何对照

12-case 结构协议和同成功样本 paired 评估如下：

| 协议 | 成功数 | 候选 Pearson/Kendall | 同成功 baseline | 结论 |
|---|---:|---:|---:|---|
| 去水 | 11/12 | -0.0260 / 0.0909 | 0.0921 / 0.2727 | 恶化 |
| cutoff 4 Å | 10/12 | -0.1112 / 0.0000 | 0.3525 / 0.4000 | 恶化 |
| cutoff 6 Å | 10/12 | -0.0599 / 0.0000 | 0.3525 / 0.4000 | 恶化 |
| fixed receptor | 12/12 | -0.0084 / 0.0000 | 0.4349 / 0.5000 | 恶化 |
| weak restraint 10 | 9/12 | 0.2001 / 0.2500 | 0.1906 / 0.2500 | 无实质改善 |
| strong restraint 1000 | 9/12 | 0.0502 / 0.2500 | 0.1906 / 0.2500 | 恶化 |

几何结果位于 `runs/fmo87/sensitivity_geometry_12/`，去水/cutoff 结果位于 `runs/fmo87/sensitivity_structure_12/`。几何变体使用原生产 OpenFE 环境；base Python 缺少 `openmmforcefields`，不能把该环境错误误报为 QM 失败。

## 4.4 PIEDA 分量诊断

`runs/fmo87/final/component_analysis/component_metrics.csv` 显示 raw tPIE 为 `0.4792/0.3201`；ES 为 `0.4398/0.3069`，DI 为 `0.4051/0.3294`，ES+CT+DI 为 `0.4668/0.3213`，tPIE−EX 为 `0.4932/0.3236`。简单替换 PIEDA 分量不能达到论文 FMO `0.55/0.38`，所以当前主要差异更可能来自 pose、片段边界和私有 scorer 定义。

## 4.7 Pose 实验边界

公开 JNK1 SDF 含有 alternate `flip` pose。对 18636-1、18629-1、18660-1 的初步运行表明，直接换 pose 会同时改变自动 pocket 与片段数，不能解释为纯 pose 效应；固定 baseline pocket 的尝试又被当前 fragment-bond 一致性检查拒绝。因此未将该结果纳入指标，也未扩展 21 个 JNK1 样本。若继续研究，必须先实现固定 fragment membership、仅替换 ligand coordinates 的底层接口。

对 Boltz-2 的实际输入/输出进一步核对：`runs/fmo87/boltz_inputs/JNK1_18636-1.yaml` 只提供 protein sequence 和 ligand canonical SMILES，没有 non-flip 或 flip 坐标，因此 Boltz-2 不是在两份公开 pose 中二选一，而是自行生成 ligand pose。将 Boltz 输出 CIF 与公开 protein 坐标对齐后，Boltz 预测相对公开 `non-flip` / `flip` 的 ligand heavy-atom RMSD 为：18636-1 `0.714/2.032 Å`、18629-1 `0.785/1.767 Å`、18660-1 `1.248/2.673 Å`；三例均更接近 non-flip。该结果支持“Boltz 输出更接近 non-flip”，但不等价于证明 Boltz 内部使用了 Ross/FEP+ 的 non-flip 坐标。

## 4.5 下一阶段优先级

1. 片段边界：shifted Cα–C′、HOP_C、ACE/NME、sequence-gap bridging 的单变量对照。
2. pose provenance：检查 source/reference/fmp 坐标、重复 pose 和 ligand heavy-atom RMSD。
3. 局部 pose 扰动：只对预先声明的少量高残差样本测试，不使用实验标签选择 pose。
4. 若片段与 pose 仍不能改善，则将差异正式归因于公开信息无法辨识的 Boltz 私有 FMO scorer 协议。

## 4.6 FACIO 边界规则对照与质量门修复

- `no_hop` 12-case：0/12 在 deck 前检查失败。关闭 HOP 后，`A:GLU:8` 片段出现奇电子（57 e⁻）却仍要求 singlet；说明当前 shifted Cα–C′ 切法不能删掉 HOP 而维持同一电子结构。
- `no_gap_bridge` 初次对照：仅 1/12 通过旧质量门；P38 2ee 出现 ligand-environment PIEDA pair `126.777 kcal/mol`，由 `EX=125.808 kcal/mol` 主导，且总 tPIE 达 `42.339 kcal/mol`。
- 将 runner 的单 pair 绝对能量阈值从 1000 收紧至 100 kcal/mol 后，以新输出目录重跑，no-gap-bridge 变体为 0/12；该 P38 2ee 结果被正确分类为 `FAILED_UNPHYSICAL_PIEDA`。
- 冻结 87-case baseline 的最大 pair 绝对值为 28.019 kcal/mol，因此该修复不影响现有 baseline 数值，只能加强后续异常碎片边界拦截。

代码入口：`src/sophosqm_baseline/fragmentation.py` 的 `use_hop`、`src/sophosqm_baseline/jobs.py` 的 `bridge_gaps`，以及 `src/sophosqm_baseline/runner.py` 的 100 kcal/mol pair 质量门。
