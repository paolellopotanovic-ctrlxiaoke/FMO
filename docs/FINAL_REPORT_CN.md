# FMO87 / Boltz-2 完整差异分析

日期：2026-09-26  
队列：Boltz-2 论文 Table 12 四靶点 87 个中性化合物  
生产目录：`runs/fmo87/final`

## 摘要

本项目的自研 official-like GAMESS FMO2-DFTB3/PIEDA 工具链已在精确 87 化合物队列上完整运行并通过严格审计。所有 87 个 FMO 和 87 个 Boltz-2 结果均成功，没有跳过样本，也没有用 MP2 替代 DFTB3。

Local Boltz-2 得到样本加权 Pearson/Kendall `0.6518/0.5096`，与论文 `0.66/0.48` 几乎一致。这有力支持队列重建、标签映射、单位转换和指标实现是正确的。

自研 raw FMO tPIE 得到 `0.4792/0.3201`，低于论文内部 FMO `0.55/0.38`，差值为 `-0.0708/-0.0599`。但分层 bootstrap 95% CI 覆盖论文值，因此不能把这解释为统计上显著的劣化。

## Table 12 对照

| 模型 | Local Pearson | Paper Pearson | Δ Pearson | Local Kendall | Paper Kendall | Δ Kendall |
|---|---:|---:|---:|---:|---:|---:|
| Boltz-2 | 0.651763 | 0.66 | -0.008237 | 0.509589 | 0.48 | +0.029589 |
| 自研 raw FMO tPIE | 0.479233 | 0.55 | -0.070767 | 0.320055 | 0.38 | -0.059945 |

Boltz-2 论文只说明 FMO 行来自作者 in-house DFTB3-based FMO code，pose 为 Glide docked poses；它没有公开该分数是否为 raw tPIE、是否加 solvent/water/geometry 约定，或是否使用 SophosQM 式 `α·tPIE+β·clogP+γ` 校准。因此本表的 FMO 对照是公开可审计 raw-tPIE reconstruction 与论文私有 scorer 的方法级对照，不是同一后处理实现的 bit-level 复现。论文没有给内部 FMO 的 PMAE/MAE/命中率，因此不对这些项做虚假直接比较。raw tPIE 的绝对 kcal/mol 误差只作为诊断；它是未校准 enthalpic score，不是绝对亲和能。

为避免把不同口径混在一起，主对照按下面的参考系解释：

| 对照 | 含义 | Pearson / Kendall |
|---|---|---:|
| Boltz-2 paper FMO | 私有 in-house DFTB3-FMO 参考值；水处理与后处理未知 | 0.55 / 0.38 |
| 本项目 raw tPIE，保留水 | 公开、无实验标签泄漏的直接 FMO 分数 | 0.479233 / 0.320055 |
| 本项目 raw tPIE，制备前去水 | 完整无 explicit-water 分支；重新制备、重新 minimization、重跑 FMO | 0.413778 / 0.276553 |
| 本项目 fitted tPIE+clogP，保留水 | SophosQM-style 逐靶点 in-sample 拟合，乐观 | 0.608842 / 0.499633 |
| 本项目 fitted tPIE+clogP，制备前去水 | 同一 in-sample 拟合用于完整无 explicit-water 分支，乐观 | 0.557731 / 0.482793 |

`raw tPIE` 不应被称为论文私有 FMO 的“正统”形式；它只是本项目可审计、无后处理泄漏的公开 scorer。`fitted tPIE+clogP` 也不能被解释为 Boltz-2 论文 FMO 的已知实现，因为论文未披露其校准方式。

## Local Boltz-2 完整指标

| Metric | Local | Paper |
|---|---:|---:|
| Pearson | 0.651763 | 0.66 |
| Kendall | 0.509589 | 0.48 |
| PMAE | 0.884080 | 0.85 |
| MAE non-centered | 0.900475 | 0.75 |
| MAE centered | 0.625489 | 0.59 |
| Within 1 non-centered | 67.82% | 69% |
| Within 1 centered | 82.76% | 83% |
| Within 2 non-centered | 93.10% | 97% |
| Within 2 centered | 97.70% | 98% |

除 P38 和 TYK2 non-centered 命中率外，多数指标与论文高度一致。TYK2 存在近似常数偏移，centering 后 1 kcal/mol 命中率从 31.25% 升至 81.25%，说明主要是绝对尺度偏差而非排序错误。

## 逐靶点 raw FMO

| 靶点 | n | Pearson | Kendall | Centered MAE |
|---|---:|---:|---:|---:|
| CDK2 | 16 | 0.796505 | 0.566667 | 9.118148 |
| JNK1 | 21 | 0.347779 | 0.186158 | 5.935048 |
| P38 | 34 | 0.568153 | 0.382416 | 9.468846 |
| TYK2 | 16 | 0.145538 | 0.116667 | 2.724375 |

总差异主要由 JNK1 和 TYK2 驱动。CDK2 的 raw tPIE 排序强于论文内部 FMO 平均参考；P38 中等；TYK2 最弱。

## Bootstrap

10000 次按靶点分层重采样：

| 模型 | Pearson 95% CI | Kendall 95% CI |
|---|---:|---:|
| Boltz-2 local | [0.512878, 0.759762] | [0.357720, 0.589927] |
| raw FMO tPIE | [0.311129, 0.630096] | [0.168225, 0.430556] |
| FMO target fit | [0.465510, 0.753041] | [0.352792, 0.587397] |
| FMO LOO | [0.252460, 0.539113] | [0.187623, 0.437859] |

论文 Boltz-2 和论文 FMO 的参考值都落在对应 local bootstrap CI 内。因此点估计差异不能被夸大为统计显著性差异。

## 计算时间

| 阶段 | 结果 | 口径 |
|---|---:|---|
| FMO scoring | mean 39.62 s/case | runner wall，不含制备和几何 |
| FMO scoring | median 36.59 s/case | runner wall |
| FMO scoring | max 61.46 s/case | runner wall |
| FMO 87 cases | sum 57.45 min | 串行 runner wall |
| FMO 87 cases | about 30 min 35 s | 四靶点并行实际 wall |
| Local Boltz-2 | 3 h 37 m 56 s | 四 target 共享 GPU 的原始推理 wall |
| Boltz-2 paper | 20 GPU sec | Table 12 方法级单预测参考 |

FMO 的 87-case 总 CPU 估计为 50.75 CPU min。Boltz-2 的本地 3 h 37 m 56 s 不能直接换算为论文 `20 GPU sec`：本地四个 target runner 同时共享 GPU 1、逐 ligand 重启进程、加载 checkpoint，并受 GPU contention 影响。日志中可见的 Boltz prediction progress 均值约 65.28 s/case，median 62 s，range 18–160 s。

## SophosQM-style 次要模型

| 模型 | Pearson | Kendall | MAE | 解释 |
|---|---:|---:|---:|---|
| raw tPIE | 0.479233 | 0.320055 | 74.240920 | 未校准 enthalpic ranking |
| in-sample α·tPIE+β·clogP+γ | 0.608842 | 0.499633 | 0.630937 | 有拟合自由度，乐观 |
| leave-one-out | 0.401952 | 0.330867 | 0.771699 | 更接近 prospective 表现 |

LOO 在 TYK2 上 Pearson 为 `-0.273756`，说明该校准形式对该系列不稳定，不应把 in-sample 结果解释为泛化能力。

### 原文与本项目口径

SophosQM 原文明确是逐 series 用实验值拟合 `α/β/γ` 后报告 fitted prediction 与实验值的相关性/误差；SI 给出各 series 拟合系数，但没有声明 LOO、cross-validation 或外部测试集。因此 `tPIE_clogP_in_sample` 是最接近 SophosQM 原文报告方式的本地结果，而不是泄漏可控的泛化估计。`tPIE_clogP_LOO` 是本项目额外加入的严格消融，用于检验同一校准形式对未见配体是否稳定，不能称为原文协议。

Boltz-2 Table 12 的 FMO 行没有公开校准细节。不能断言它是 raw tPIE，也不能断言它完整采用了 SophosQM 的 `tPIE+clogP` 回归；只能说它是 in-house DFTB3 FMO scoring baseline。因此本项目把 raw tPIE 作为主 FMO 对照，把 fitted/LOO 版本作为敏感性分析。

## 校准误差消融

官方 SophosQM SI 给出六个原始 series 的 `α/β/γ`：DNA ligase `0.074/-0.456/1.856`、MUP-I `0.076/-2.079/-3.437`、HSP90 `0.092/0.045/0.837`、p38 `0.030/-1.107/0.223`、JAK2 `0.0840/-0.643/0.773`、MCL-1 `0.0054/-0.9418/-1.31116`。这些系数不能直接移植到 Boltz2 四靶点 cohort：原始方法还使用不同配体系列、Glide/MacroModel 几何、VolSurf clogP 和 FMO2-MP2/6-31G*。

可复现消融见 `runs/fmo87/final/CALIBRATION_ABLATION_CN.md`：

| 消融 | Pearson | Kendall | 结论 |
|---|---:|---:|---|
| raw tPIE | 0.479233 | 0.320055 | 未校准 FMO 排序 |
| tPIE-only in-sample | 0.479233 | 0.320055 | 正向共同仿射校准不改变排序 |
| clogP-only in-sample | 0.481288 | 0.380918 | clogP 含靶点内信息 |
| tPIE+clogP in-sample | 0.608842 | 0.499633 | 拟合自由度带来乐观提升 |
| tPIE+clogP LOO | 0.401952 | 0.330867 | 系数外推不稳定 |
| 官方 p38 系数迁移到本地 p38 | 0.589925 | — | 只比 raw p38 0.568153 小幅提高 |

因此，与 Boltz2 论文 raw FMO 的差距主要不在“缺少 SophosQM 拟合参数”；更可能来自 pose、几何、碎片/溶剂协议和私有 DFTB3 FMO scorer。SophosQM 校准能改善绝对 kcal/mol 误差，但 in-sample 排名提升不能外推，LOO 反而低于 raw tPIE 的 Pearson。

## 显式水消融

为检验 explicit water 是否解释与论文内部 FMO 的差距，项目完成两个互补分支；两者均为 87/87 strict complete，并保持 water PCM 开启。

| 策略 | raw tPIE Pearson | raw tPIE Kendall | P38 Pearson | P38 Kendall |
|---|---:|---:|---:|---:|
| 保留 SOP-style crystal water | 0.479233 | 0.320055 | 0.568153 | 0.382416 |
| 只从 FMO pocket 去水 | 0.475757 | 0.315845 | 0.559260 | 0.371643 |
| 制备前去水并重新优化 | 0.413778 | 0.276553 | 0.400667 | 0.271102 |

SophosQM-style `tPIE+clogP` in-sample 结果同样下降：保留水为 0.608842/0.499633，scorer 级去水为 0.602041/0.488406，制备级去水为 0.557731/0.482793。

生产输入中只有 P38 含 20 个 selected crystallographic waters；其他三个靶点没有水。P38 成对 tPIE 相对保留水的变化为：

| 策略 | mean Δ(kcal/mol) | median Δ(kcal/mol) |
|---|---:|---:|
| 保留水 | 0.000000 | 0.000000 |
| 只从 FMO pocket 去水 | +20.034500 | +17.571000 |
| 制备前去水并重新优化 | +16.862647 | +17.076500 |

解释边界必须区分两个分支。`removed_from_fmo_pocket` 复用保留水体系的 minimized complex，只在 5 Å pocket/fragment 构建时移除水，因此检验 scorer 中 explicit-water 的直接贡献；`removed_before_preparation` 在 PDB2PQR/PROPKA 前移除水，随后重做 OpenMM restrained minimization，因此同时检验质子化网络、氢键网络和几何变化。二者都没有关闭 water PCM，所以不是 implicit-solvation 消融。

结论：explicit water 不是当前与 Boltz-2 论文 FMO 指标差距的主要原因。去掉水没有改善任何主指标；制备前去水还明显破坏 P38 raw-tPIE 排序。生产协议应继续保留已有 crystal water 且不新增水。本地 Boltz-2 输入没有 explicit water，但论文未披露其私有 FMO baseline 是否保留水，因此不能把双方水协议假设为一致。

## 与论文内部 FMO 差异的合理来源

1. **Pose 来源不同**：论文内部 scorer 使用 Glide docked poses；本项目使用公开 Ross/Chen reference poses。端点相互作用对 pose 极其敏感。
2. **几何工具不同**：论文描述的私有管线不可得；本项目用 OpenMM AMBER14/OpenFF restrained minimization 替代 MacroModel/OPLS3e。
3. **Scorer 非同一软件**：本项目是独立 GAMESS FMO2-DFTB3/PIEDA 实现，不使用 Boltz 作者或 Schrödinger 私有 FMO scorer。
4. **水的处理不同**：自研 FMO 协议保留已有晶体水、不新增水，并对 DFTB3 dimer 使用 water PCM 以稳定极性片段；生产审计显示 34/87 个 case 有 explicit water（全部在 P38，共 403 个，每个 P38 case 10–16 个），CDK2/JNK1/TYK2 为 0 个。显式水消融显示去水不改善指标，因此它不是主要差距来源。本地 Boltz-2 输入只含 protein sequence、ligand SMILES 和 MSA，不含 explicit water；Boltz-2 论文没有公开其内部 FMO baseline 是否保留 crystal water。
5. **绝对尺度未校准**：raw tPIE 与实验 ΔG 之间有几十 kcal/mol 偏移；论文 Table 12 的 FMO 比较只报告排序相关。

## 结论

这套自研工具链已经完整 work：结构制备、标准几何、化学分片、DFTB3 FMO2/PIEDA、质量门、Boltz2 复现实验、精确队列匹配、Table 12 式指标和显式水消融全部闭环。相对于论文内部 FMO，raw tPIE 点指标偏低，主要暴露在 JNK1/TYK2，但 bootstrap 不支持显著劣化结论。Local Boltz2 的强复现说明差异更可能来自 FMO pose/scorer 协议，而不是队列、评估器或 explicit-water 处理错误。
