# FMO 自适应方法最终说明

日期：2026-10-01。本文是当前已完成证据的工程交付说明；它区分“已验证的
生产规则”“探索性候选”和“仍未完成的自由能项”，不把端点 FMO 分数称为
完整的结合自由能。

## 1. 固定的处理顺序

1. 冻结 target、配体身份、PDB/SDF 哈希、质子化、净电荷、晶水策略和口袋定义。
2. 审计所有 receptor 原子和配体显式氢，删除全队列共同冲突水，重新生成配体氢。
3. 对实际最小化后的 complex 再做硬碰撞和位移审计；若修复或系列尺度异质性
   明显，整组运行 `k=1/10/100`，不能逐 ligand 后验改参数。
4. 在每个 arm 中保存四个独立候选：
   `tPIE`、`ES`、`ES+SOLV`、`ES+CT+SOLV`。
5. 做 PIEDA pair、SCF、输入哈希和能量闭合检查；失败 case 不得静默丢弃。
6. 若需要组合多个候选，先在 target 内用预先冻结的 robust-z 或 rank 聚合处理
   尺度，再做等权聚合；这只是排序数学，不是把不同物理模型的 kcal/mol 直接相加。
7. 只有完整 bound/free 采样通过身份、ESS、split-half JS 和重复 replica 门时，
   才允许产生熵分支；否则追加预声明 sampling block，预算耗尽则 `OMIT`。
8. 输出候选分数、跨 arm 稳定性和不确定状态。稳定性不能替代准确性验证。

## 2. 分数的物理关系

设蛋白–配体 PIEDA pair 的五项为 `ES, EX, CT, DI, SOLV`：

\[
tPIE=ES+EX+CT+DI+SOLV,\qquad ES+SOLV\subset tPIE.
\]

因此 `ES+SOLV` 不能再和 `tPIE` 相加；若将两个候选混合，必须预先写成
独立的数学模型并验证，而不是声称能量项相加。PCM、PBSA、三态差值和熵也
只有在参考态、溶剂模型和重复计算边界一致时才可构造新的自由能代理；当前
它们保留为独立诊断分支。

## 3. 无标签自适应规则

- 结构硬碰撞或最小化后异常：`REPAIR_REQUIRED`，修复后整组复算。
- 修复、重原子 review、位移异常或系列平均重原子数差异至少 15%：运行三种
  restraint，并进入 `FMO_REVIEW` 候选流程。
- 跨 restraint 的秩稳定性只用于判断数值复现；若稳定性通过，仍标记
  `UNVERIFIED_LABEL_FREE`，除非已有独立新靶点证据证明方向规则。
- 若四个零拟合候选都缺乏可靠方向证据：输出 `FMO_UNCERTAIN`，允许独立
  Boltz2 fallback，但绝不把 Boltz2 分数与 FMO kcal/mol 相加。
- 任何系列特化归一化、后验权重或按实验标签选分支，都必须进入明确的训练/验证
  协议，不能作为零拟合生产方法。

## 4. 当前证据矩阵

| 数据集 | 结构/QM | 当前最可信状态 | 结果 |
|---|---|---|---|
| 原 87 例 | 87/87 成功，PIEDA 审计通过 | `FMO_READY` | 固定 `ES+SOLV`：Pearson `0.633634`，Kendall `0.507750` |
| Galectin3 | E68/E69/E70 统一删水、补氢、最小化和 QM 输入；26/26×3 QM/pair 闭合 | `FMO_REVIEW` | 三臂稳定性仅作为数值审计；`ES+CT+SOLV` 保留为 review candidate，未验证方向前生产候选仍为 `ES+SOLV`，不能部署系列特化归一化 |
| FXa | E71 在删冲突水、补氢、统一最小化后，使用显式记录的旧 E40 参考口袋完成 27/27 QM；当前受体结构门为 `PASS_AFTER_MINIMIZATION` | `FMO_REVIEW` | tPIE/ES/ES+SOLV/ES+CT+SOLV 的 Pearson/Kendall 分别为 `0.570/0.386`、`0.607/0.342`、`0.606/0.371`、`0.589/0.365`；E71 是参考口袋诊断臂，标准新口袋仍需独立验证，不能直接冻结为生产协议 |

Galectin3 的 18 系列不是输入损坏：它与 17 系列均为中性且 QM 成功，但平均
重原子数约为 `35.25` 对 `29.61`，所以被标记为系列异质性 review。三臂稳定性
选择器仍可能选到真实排序反向的分支，证明稳定性不是准确性。

### 4.1 当前因素消融的直接结论

以下数值来自同一 `fmo_method_matrix.json`，只用于解释已经完成的解盲审计，不能
反过来作为无标签选分支的依据：

| 因素/分支 | 原 87 例 | Galectin3 修复三臂 | FXa E71 参考口袋 | 当前解释 |
|---|---:|---:|---:|---|
| `tPIE` | `r=0.533`, `τ=0.412` | `r≈0.20–0.24`, `τ≈0.17–0.24` | `r=0.570`, `τ=0.386` | 直接蛋白–配体 PIEDA 相互作用；Galectin3 中是唯一保持正方向的 FMO 分支，但强度弱 |
| `ES` | `r=0.566`, `τ=0.477` | `r≈-0.49–-0.69`, `τ≈-0.45–-0.50` | `r=0.607`, `τ=0.342` | 静电项在 Galectin3 中方向失配，不能单独作为通用主分数 |
| `ES+SOLV` | `r=0.634`, `τ=0.508` | `r≈-0.43–-0.54`, `τ≈-0.39–-0.39` | `r=0.606`, `τ=0.371` | 原 87 例当前最佳冻结分支；Galectin3 中溶剂响应没有修复方向，不能普适外推 |
| `ES+CT+SOLV` | `r=0.644`, `τ=0.504` | `r≈-0.41–-0.49`, `τ≈-0.30–-0.37` | `r=0.589`, `τ=0.365` | 解盲后看似略优于原 87 的 `ES+SOLV`，但不能因此改写零拟合生产规则 |
| 多 restraint (`k=1/10/100`) | — | 数值重复性通过 | 当前 E71 不是同一标准 pocket | 稳定性只能证明复现性，不能证明实验方向 |
| robust-z 聚合 | `r=0.476/τ=0.391`（三分支） | `r=-0.068/τ=-0.065`（tPIE+ES+SOLV） | 约 `r=0.605/τ=0.411` | 只消除尺度差异；不能修复方向、结构或参考态问题 |
| 熵 proxy | 未进入生产 | 未完成 | 未完成 | 只有 bound/free 完整收敛后才可做增量消融 |

这张表给出的不是“Galectin3 应该切换到 tPIE”的部署规则，因为那会使用实验
标签进行后验选择。它支持的工程结论更窄但更可靠：`tPIE`、`ES` 和溶剂相关分支
确实量度不同物理内容；加入溶剂项在原 87 例有明显增益，但在 Galectin3 不能
消除方向失配。因此当前通用流水线应输出全部分支和 review 状态，而不是假设一个
分支在所有靶点都同样有效。

无标签多臂审计现在还会记录跨 arm 的 pairwise Kendall-like 一致性和 robust-z
偏差。Galectin3 的 `k=1/10/100` 结果为：`tPIE` 的中位一致性 `0.691`、偏差
`0.276`；`ES` 为 `0.588/0.301`；`ES+SOLV` 为 `0.815/0.165`；
`ES+CT+SOLV` 为 `0.870/0.181`。但同一队列中后两者对实验排序仍为反向，说明
“多臂稳定”只能作为数值复现门，不能作为准确性门。FXa 目前只有 E71 参考口袋
一个 arm，pairwise 一致性按规则报告为不可估计，而不是伪造稳定性证据。

## 5. 熵状态

E55 的显式两态采样未通过收敛门，当前熵不进入生产分数。E61 的 CDK2 ligand 17
在统一 `OpenMM 8.2.0 + CUDA 12.4` runtime 下已完成第一轮 `24×3×2=144`
blocks；正式审计的 bound/free min ESS=`1319.4/470.0`，max split-half
JS=`0.355/0.263`，仍未通过 `JS≤0.15`。点估计 `-TΔS=+1.4021 kcal/mol` 不能进入
排名。控制器已追加独立 blocks 25–48，当前 48-block 期望总数为 288，尚未完成，
故熵保持 `ENTROPY_OMIT/PENDING`。最终审计还会报告按 block bootstrap 的保守区间。
即使扭转 proxy 通过数值门，它也不等于完整结合熵，不能单独宣称为 `-TΔS_bind`。

## 6. 适配判定的数学边界

对一个新 target，结构/QM/姿态门只决定“是否可信地执行和报告”，不能从无标签输入
推导实验亲和力方向。具体地：

1. 若结构硬碰撞、化学身份或 QM coverage 失败，输出 `STRUCTURE_REVIEW`/`FMO_FALLBACK`，
   不输出完整排名。
2. 若多个 restraint/pocket arms 的 rank concordance 低于预设门，保留多臂不确定度，
   不按单个配体选 arm；这测的是数值敏感性，不是准确性。
3. 若 arms 稳定但没有独立方向验证，仍只能输出冻结的候选分支和 `FMO_REVIEW`；
   稳定性不能证明 Galectin3 式反向排序是正确的。
4. 若预先声明的候选需要组合，只能在 target 内做
   \(z_j=(S_j-\operatorname{median}S_j)/(1.4826\operatorname{MAD}S_j)\)，再做
   等权 rank aggregation；这不是 kcal/mol 加法，也不能用实验标签选权重。
5. 只有独立 prospective/外部方向验证才能把 review candidate 升级为 production
   branch；否则使用独立 fallback，并明确 FMO 与 fallback 不相加。

## 7. 可复现入口

- 主流程：`docs/FMO_ADAPTIVE_PIPELINE_CN.md`
- 判定树：`docs/FMO_ADAPTIVE_PRODUCTION_DECISION_CN.md`
- AnewFEP 边界：`docs/ANEWFE_VALIDATION_AUDIT_CN.md`
- 完整方法矩阵：`runs/fmo87/ablation/final_analysis/fmo_method_matrix.json`
- 按 target 的 readiness 矩阵：`runs/fmo87/ablation/final_analysis/fmo_target_readiness.json`
- 按 target 的物理因子决策表：`runs/fmo87/ablation/final_analysis/fmo_factor_decision_table.json`
- 跨数据集因素/评分总表：`runs/fmo87/ablation/final_analysis/universal_factor_summary.json`
- 统一 target profile：`runs/fmo87/ablation/final_analysis/fmo_target_profiles.json`
- 生产 Runbook：`docs/FMO_ADAPTIVE_RUNBOOK_CN.md`
- Galectin3 化学/QM 审计：`runs/fmo87/ablation/final_analysis/galectin_chemistry_audit.json`
- 自适应选择压力测试：`runs/fmo87/ablation/final_analysis/galectin_adaptive_selection_audit.json`
- 熵采样：`runs/fmo87/ablation/E61_explicit_two_state_20ns_CDK2_etkdg_free/`
- 熵增量门控：`scripts/score_fmo_entropy_increment.py`（仅在完整覆盖且通过收敛门后运行）
- 尺度聚合消融：`runs/fmo87/ablation/final_analysis/scale_aggregation_audit.json`

E68 与旧 E57 的配对重算显示，四个分支的最大单例变化仅约 `0.125 kcal/mol`，
整体 tPIE 变化约 `-0.012 kcal/mol`；因此把 E57 的受体/SDF provenance 修正为
一致输入并没有消除 Galectin3 的反向 ES/溶剂分支。当前最合理的解释已从“单纯
最小化受体错配”收窄为：Galectin3 的 pose/系列化学状态与这些端点相互作用分支
之间仍存在未解决的方向失配，不能靠继续加权掩盖。

当前最重要的结论不是“一个分数适合所有靶点”，而是已经固化了可迁移的输入
处理、质量门、结构修复、分支审计、熵收敛门和失败策略；对于纯 FMO 在某个
外部靶点没有方向证据的情况，流程会明确降级，而不是用后验调分掩盖失败。

### 多构象 ensemble 的固定规则

多构象采样允许作为降低单一静态 pose 敏感性的显式消融，但不能从多个结果中挑选
一个最有利构象。构象生成器、数量、随机种子、受体/配体 fingerprint、权重和聚合
公式必须在解盲前冻结；每个 ligand 的全部预注册构象都必须通过 QM 质量门并纳入。
`scripts/aggregate_fmo_conformer_ensemble.py` 默认使用固定 weighted mean，并同时
报告 median、Boltzmann 诊断和逐例 min/max 区间。缺失构象、部分 QM 失败、混合
Hamiltonian 或结构 fingerprint 时直接 `FMO_REVIEW`，不能用幸存构象补齐。该规则
已经实现但尚未在各数据集完成运行，因此不能提前宣称会提高跨集合指标。
FXa 的复核暴露了一个重要流水线问题：旧结构门只检查了部分受体信息，新全原子
门发现 HOH306/309/372/379 与配体存在硬冲突。E62 已删除这些整水残基、重新
生成 27 个配体氢并完成 27/27 最小化；修复后的标准新口袋在原 solver、
DIIS/charge-DIIS、GAMESS-SOSCF、PCM、无 PCM 和 `k=100` 试验中曾未形成可用 QM
结果。E64 的 `MAXIT=500` 被 GAMESS 合法性门拒绝，并未进入科学计算。随后 E71
使用 `frozen_protocol.json` 中明确记录的旧 E40 参考口袋完成 27/27 QM，证明修复
后的结构和 QM 本身可以闭合，同时把问题收窄为“新姿态/新 union-pocket 拓扑与旧
口袋不一致”的结构性诊断问题。E71 因此只能标记 `FMO_REVIEW`：不能冒充标准
生产协议，也不能继续把 FXa 概括成“完全无法计算”。

进一步对 E62/E63/E65 的第一个失败样本 `5f` 做了日志级诊断：E62 实际已经生成
27 个最小化 complex（`post_minimized_audit.json` 为 `27/27`），所以失败不能归因于
“忘记最小化”。失败发生在标准 5 Å union-pocket 的 FMO iteration 201，日志明确报告
`layer 1` 的 monomer SCF 未收敛。改变 solver 后仍未形成可用结果：E62 SOSCF 的
`maxD=0.007528`，E63 DIIS 的 `maxD=1.9350988`，E65 GAMESS SOSCF 试验的
`maxD=8.5409213`；gas/k100 pilot 也在同一类 monomer-SCF gate 失败。由此当前
可排除“仅由受体最小化遗漏造成”的解释，问题更具体地收窄为标准 pocket 的局部
片段环境/电荷响应与当前 DFTB3/PCM SCF 的不稳定组合。下一步应定位失败 fragment
及其局部接触，再做结构级或片段边界级修复；在此之前不能把 E71 旧 pocket 的
成功结果升级为标准 FXa 生产证据。

## 8. 独立 HF 三态复核（2026-10-05）

为检验气相 `FMO2-HF/STO-3G` 三态差值，原 87-case 队列以同一 target 共用一个
receptor 态，并分别计算 ligand 与 complex：

\[
\Delta E_{3state}=E_{AB}^{FMO2-HF/STO-3G}-E_A^{FMO2-HF/STO-3G}
 -E_B^{FMO2-HF/STO-3G}.
\]

仅当 A、B、AB 三个态均成功且输入 manifest 协议一致时，case 才算完整。排序指标
仅对某 target **全部 cohort case** 都完整时计算；部分完成的 target 只报告覆盖率，
不按幸存样本计算相关性。这个规则避免 QM 超时把易算样本选择性留下而造成偏差。

截至 2026-10-05 02:39 UTC，HF 三态覆盖为 `24/87`：CDK2 `16/16`，JNK1 `8/21`，
P38 `0/34`，TYK2 `0/16`。CDK2 是当前唯一具备完整靶点覆盖的 HF 子集，其
未拟合回顾性指标为：

| 分数 | Pearson r | Spearman ρ | Kendall τ |
|---|---:|---:|---:|
| HF 三态 ΔE | 0.398 | 0.418 | 0.333 |
| tPIE | 0.797 | 0.753 | 0.567 |
| PIEDA ES | 0.852 | 0.871 | 0.700 |
| PIEDA ES+SOLV | 0.875 | 0.868 | 0.700 |
| PIEDA ES+CT+SOLV | 0.879 | 0.874 | 0.700 |

这里的 `ES+SOLV` 与 `ES+CT+SOLV` 是同一复合物姿势中 PIEDA 分量的直接和，
不是三态 HF 分解项，也没有与 tPIE 再相加。固定种子、5000 次配对 bootstrap
得到 HF 三态 Kendall τ 的 95% 区间约 `[-0.063, 0.704]`；它相对 tPIE 的配对
Δτ 为 `-0.233`，95% 区间约 `[-0.613, 0.144]`。区间较宽且跨 0，因此当前证据
只支持“CDK2 这 16 例上，三态 HF 点估计没有超过 PIEDA ES 分支”，不能判定总体
优劣，也不能外推到其它靶点。所有标签只用于事后评估，没有参与结构、参数、分支
或权重选择。

计算与审计入口：`scripts/run_hf_fmo87.py`（可断点续跑）、
`scripts/audit_fmo_hf_coverage.py`（包含共享 receptor 状态）、
`scripts/audit_fmo_independent_three_state.py`、
`scripts/analyze_fmo_hf_rankings.py` 和 `scripts/monitor_fmo_hf_queue.py`。
实时证据位于 `runs/fmo87/ablation/final_analysis/hf_sto3g_live_coverage_audit.json`、
`hf_sto3g_live_three_state_audit.json`、
`hf_sto3g_complete_target_ranking_audit.json`；87-case QM 仍在运行，不能将本节
结果表述为完整 HF 排序或普适精度结论。
