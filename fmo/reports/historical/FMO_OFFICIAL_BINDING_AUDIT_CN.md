# GAMESS 官方 FMO 结合能算法对照与 87 例根因审计

日期：2026-09-27。本报告检验实现是否符合 GAMESS 的结合能定义，并分析评分协议的局限；以下各项**没有根据实验标签拟合参数**。

## 官方定义：能量、几何、子体系

- 仓库内 GAMESS 官方手册 `software/gamess/docs-input.txt:15730`：`MODMOL=8` 将片段分成子体系，额外加 `16`（即 `MODMOL=24`）可计算 AB 结合能；`MOLFRG=(1,…,1,2)` 指定 A 和 B。`$ISAXYZ/$ISBXYZ` 可指定游离 A/B 几何；**不指定时直接沿用复合物中的几何**。因此“固定结合构象 B”和“优化后的游离 B”是两个不同协议。
- `software/gamess/docs-references.txt:10265`：选择性 FMO 跳过部分二聚体，不提供完整体系总性质。不能拿旧 `MODMOL=1` tPIE 日志的总能量做 AB−A−B；本次 AB 使用完整 FMO2，B 使用独立 FMO1。
- `software/gamess/source/fmoio.src:7810`：官方原生子体系输出分别给出单体状态变化、去溶剂化、相互作用以及总 *binding energy* 和独立的 `Check binding energy`。`software/gamess/source/fmoio.src:4162`：PCM 完整溶剂项被加至能量，故本协议取 `EunD+so` 而非只含受溶剂极化溶质能量的 `Eunco+D`。`software/gamess/docs-input.txt:13176` 定义 `IFMO=-1` 的 PCM<1> 模式；`software/gamess/docs-input.txt:15812` 建议 FMO-DFTB 分析梯度用 `MODGRD=32`。

## 官方 GAMESS 原生数值交叉验证

`scripts/validate_fmo_official_subsystem.py` 根据原 AB deck **独立调用** GAMESS 的 `MODMOL=24`，直接读官方 `Check binding energy`；以下是和现有三项独立运算 `E(AB)−E(A)−E(B)` 的对比：

| 靶点/配体 | 独立差值 kcal/mol | 官方子体系 kcal/mol | 绝对差 |
|---|---:|---:|---:|
| CDK2 / 20（原结合构象 B） | −34.214492 | −34.214495 | 0.0000034 |
| P38 / 2aa（原结合构象 B） | −41.287344 | −41.287348 | 0.0000045 |
| TYK2 / jmc_30（原结合构象 B） | −33.678184 | −33.678188 | 0.0000043 |
| CDK2 / 20（`$ISBXYZ` 优化游离 B） | −6.073357 | −6.073363 | 0.0000057 |
| P38 / 2p（`$ISBXYZ` 优化游离 B） | +137.467734 | +137.467740 | 0.0000063 |
| CDK2 / 32（`$ISBXYZ` 优化游离 B） | +77.219431 | +77.219430 | 0.0000015 |

原生分项例：CDK2/20 的单体项 −9.705、去溶剂化项 +58.181、相互作用项 −82.691，合为 −34.214 kcal/mol；P38/2aa 对应 −5.213、+52.396、−88.470，合为 −41.287。上表支持的结论是：**这些极端值不是减号、片段选择或能量解析脚本写错。**原始日志与对照 JSON 位于 `runs/fmo87/relative_total_energy_official_check/` 和 `runs/fmo87/relative_total_energy_official_relaxed_check/`。

另外检查准备结构：`runs/fmo87/final_fmo_v3_dftb3_pcm/minimized/CDK2_20/minimization.json` 明确记载固定蛋白/水、约束配体；比对全部 87 个样本，蛋白重原子与参考准备蛋白的坐标一致。CDK2/20、P38/2aa 的最近蛋白—配体重原子距离约 2.57、2.81 Å；没有证据表明本批数据是柔性蛋白姿势误拼接造成的严重碰撞。

## PCM 分项对排序的影响

`scripts/analyze_fmo_total_energy_components.py` 从原始 178 份结果日志中，分别对官方打印的三种量做 **全 87 例、各靶点内** 无拟合排序：

| 同一批 FMO/PCM 日志中的量 | 加权 Pearson | 加权 Kendall |
|---|---:|---:|
| `Eunco+D`：受 PCM 极化的溶质，但尚未加 Gsol | 0.227 | 0.209 |
| `EunD_es`：加入 PCM 静电溶剂项 | 0.047 | 0.102 |
| `EunD+so`：加入完整 PCM 溶剂项（原方法） | −0.018 | 0.043 |

**第一行不是单独的气相计算**；它的分子电子状态是在 PCM 存在时计算的。详细逐例数据在 `runs/fmo87/relative_total_energy/component_diagnostics.csv`，汇总在同目录 `component_diagnostics.json`。这些数据证明 PCM 溶剂项足以改变排序，**不能**将其中指标较高者回选为事先定义的独立测试模型。

## 官方允许的自由 B 构象：全 87 例检验

之前 B 直接从 AB 提取，仍处于结合构象。`scripts/run_fmo_ligand_relaxation.py` 从这个姿势出发对 87 个游离 B 进行原生 GAMESS FMO1-DFTB3/PCM `RUNTYP=OPTIMIZE`、`MODGRD=32`、`OPTTOL=0.001`、最多 160 步。仅当日志**明确显示 EQUILIBRIUM GEOMETRY LOCATED**才接受；全部 87/87 收敛。`scripts/audit_fmo_ligand_relaxation.py` 独立核查 87 份日志、元素/原子顺序、最终优化能量、输入/可执行文件/DFTB 参数哈希，以及整个队列身份：`runs/fmo87/relative_total_energy_ligand_relaxed/audit.json`。注意：这仅是从原姿势出发的**单一局部极小值**，不是游离配体热力学系综。

| 靶点 | 原结合构象 B Pearson / Kendall | 优化 B Pearson / Kendall | 中位能量放松 kcal/mol |
|---|---:|---:|---:|
| CDK2 | −0.111 / −0.100 | −0.611 / −0.333 | 51.84 |
| JNK1 | 0.096 / 0.196 | −0.463 / −0.310 | 13.97 |
| P38 | −0.366 / −0.206 | −0.382 / −0.246 | 18.84 |
| TYK2 | 0.668 / 0.517 | 0.666 / 0.383 | 16.76 |
| **87 例加权** | **−0.018 / 0.043** | **−0.251 / −0.162** | — |

CDK2/20 放松后的官方 `$ISBXYZ` 从原 −34.214495 变为 −6.073363 kcal/mol（比实验 −8.72 更接近），但是**87 例整体更差**。P38/2p 的放松代价约 140.89 kcal/mol，导致官方原生结合能 +137.47 kcal/mol；CDK2/32 的代价约 97.86，导致 +77.22。对 87 个优化结构按 SDF 原键检查，未发现重原子共价键长度超过 2.2 Å；P38/2p 的重原子（未对齐）RMSD 约 0.42 Å。这并不排除不正确的自由配体质子微态、构象极小值或某些几何与 DFTB3/PCM 能量模型不一致。

## 经过排查后的边界

1. **排除：**当前代码把旧 selective FMO 能量当总能量、取错 `EunD+so`、AB/A/B 减法符号错误、原始蛋白姿势互不一致。真实异常样本的官方原生子体系结果也支持同一个极端能量值。
2. **确认的协议缺口：**原方法省略了 B 从结合构象到游离极小值的变形代价；一旦以该 QM/PCM 等级计入，这一项往往非常大且有化学系列依赖，但简单补全反而降低本队列排序能力。由此**不能**宣称已通过游离 B 优化修复亲和力预测。
3. **尚未排除：**AB 与 A 也未在一致的 QM/PCM 哈密顿量下重整，受体截断/cap/PCM 空穴敏感，缺少水置换、多构象、自由配体微态和结合/未结合熵。当前不足以判断哪个单独问题是根因，更不能把静态 `ΔE` 解释成实验 `ΔΔG`。如要改进，应先锁定小规模外部验证集，在相同 QM 等级审计 AB、A、B 的几何与溶剂项，再评估成本和统计稳定性；不应在 87 个测试标签上追选某个已打印的能量分量。
