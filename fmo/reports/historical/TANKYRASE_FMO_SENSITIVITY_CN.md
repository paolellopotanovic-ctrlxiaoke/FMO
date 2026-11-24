# Tankyrase 2 MFMO 论文对照与固定几何 FMO 三项差敏感性

日期：2026-09-27。本文区分**原文确实做过的计算**、**本地完成的计算**和**尚未完成的计算**；严禁把已经反复查看的 87 例当成独立外部测试集。

## 1. tPIE 是否更有亲和力的物理意义？

不是。tPIE 是选定结合姿势上的配体–口袋片段相互作用之和；实验亲和力是溶液中结合/未结合系综的自由能差。当前 tPIE 的 87 例目标内样本数加权 Pearson/Kendall 为 `0.479/0.320`，这是**此数据集上的经验排序关联**，不表示 tPIE 已是可靠的跨靶点或全新化学系列自由能预测器，更不意味着 tPIE 比完整 `AB−A−B` 在热力学上更接近自由能。某个物理项加进静态分数后能否提升经验相关性，取决于该项是否可靠、参考态是否匹配以及误差抵消。

## 2. 原始论文到底算了什么

来源：*Use of the Multilayer Fragment Molecular Orbital Method to Predict the Rank Order of Protein–Ligand Binding Affinities: A Case Study Using Tankyrase 2 Inhibitors*，ACS Omega 2018，DOI `10.1021/acsomega.8b00175`，正文 `Materials and Methods / Scoring Function Used for the Prediction of Protein–Ligand Binding Affinity`。原文单个优化结构上的 MM-PB/SA 风格评分为：

`ΔG_pred ≈ ΔE_bind,gas(MFMO) + ΔG_solv,PB/SA − TΔS_bind,gas`

其中 `ΔE_bind,gas = E_gas(AB)−E_gas(A)−E_gas(B)`；`ΔG_solv,PB/SA = G_solv(AB)−G_solv(A)−G_solv(B)`，极性项用 AMBER12 MM-PB（内/外介电常数 1/80，网格 0.5 Å，指定 AMBER 原子半径/电荷），非极性项用 SASA（`SURFTEN=0.00542 kcal/mol/Å²`、`SURFOFF=0.92 kcal/mol`、探针 1.4 Å）；旋转键熵代价 `+α×Nrot kcal/mol`，作者对照 `α=0` 和 `α=1`，后者值来自既有文献，不是精确统计力学熵。**论文不是直接把一次 PCM 计算的 `EunD+so` 叫 PB/SA，也没有同时使用 PB/SA 和 PCM 溶剂项。**

前处理同样重要：AMBER/TIP3P 显式水中先做 MM 最小化，移除水后对配体和部分位点做 ONIOM 局部 QM 几何优化；对柔性 class B 还做多初始构象尝试。打分阶段每个候选使用单一优化结构，**并没有完成自由态/结合态热力学系综采样**；原文也承认其绝对亲和力明显高估。原文对 class A 报告优化前后单独 MFMO 气相结合能 Pearson `0.557→0.713`，结合 PB/SA 的普通 FMO 在其靶点上 `R=0.900, PI=0.887`；普通 FMO+PCM 排除非极性项 `R=0.870, PI=0.870`，完整 PCM `R=0.507, PI=0.251`，其作者认为该体系 PCM 色散项影响很大。**PI 不是 Kendall，也不能与我们的四靶点指标直接对号入座**。这些是同体系方法选择的关联性结果，不是冻结参数后的外部新靶点验证。

## 3. 我们能做到什么、不能冒充什么

我们当前是 GAMESS FMO2-DFTB3/3OB-3-1、固定蛋白/配体姿势、保留指定口袋水、GAMESS 水 PCM `EunD+so`。原文使用 MFMO、气相 QM 结合能、独立 MM-PB/SA 溶剂、ONIOM 几何，因此当前实现既**不是原论文复现**，也不是比它更完整的自由能模型。截至本报告第一次实验，没有验证论文式 AMBER12 `MMPBSA.py`/`pbsa`/ONIOM 工作流；**2026-09-28 后续盘点发现环境中存在可运行的原生 AmberTools `pbsa -pqr`**（但 `MMPBSA.py` 缺失 Python 模块），进一步开展独立 PQR-PBSA 单结构试验，详见 `TANKYRASE_PBSA_FOLLOWUP_CN.md`。不能拿现有 PCM 分项**更名**成 MM-PB/SA。

我们可以严格地做三类研究：①已审计 87/87 PCM 能量分项和固定 `+1 kcal/mol/严格 RDKit 可旋转键` 的熵代理；②独立重跑原子/参数哈希隔离的 PCM 仅静电、真正无 PCM 气相、无显式水、6 Å 口袋、UFF 色散等预先指定变体；③记录每一变体的失败率，只在同成功样本上与原 PCM 分数做配对比较，不准把 3/12 的成功指标伪装成 87/87。固定 `+1` 仍然是一个引入先验的经验模型；RDKit 的严格可旋转键定义不等同论文的原子类型规则。

## 4. 冻结的本地实验矩阵与解释边界

原有的 tPIE 协议敏感性试验在 `docs/FMO_SENSITIVITY_CN.md`；**不得代替新三项能量差试验**。本轮不按标签挑化合物：试验集为 `data/fmo87_cohort.csv` 中每个靶点前 3 个配体（共 12 例）；原则上只有数值稳定且不混淆模型的方案才值得在另一队列进行全量验证。五个独立变体只改一个旋钮：

| 变体 | 相对原始三项 PCM 能量差的唯一改动 | 意义 |
|---|---|---|
| `pcm_polar_only` | PCM `ICAV=0 IDISP=0` | 真正重跑关闭空腔/色散项的独立 deck，不是读同一日志的一列 |
| `gas` | `use_pcm=False` | 真正重跑无 PCM；收敛失败必须保留 |
| `no_waters` | 不保留蛋白侧结构水 | 显式水处理，不是切换介电常数 |
| `cutoff_6` | 全队列共有口袋半径 5→6 Å | 截断体系改变，蛋白 A 仍同靶点共享 |
| `dftb_uff` | DFTB `DISP=UFF` | QM 色散更改，不等于 PCM `IDISP` |

`scripts/run_fmo_total_sensitivity.py` 根据固定队列在新目录分别调用三项全 FMO2 runner；独立参数和选中配体 manifest 不容混用；`scripts/evaluate_fmo_total_sensitivity.py` 逐成功计算校验 deck、GAMESS 二进制、全部 DFTB 参数哈希和原始日志能量，只允许同靶点至少 3 例计算相关性；所有负向/未收敛/不完整结果保留。

执行方式（`$PY` 为已有 RDKit/NumPy/SciPy 环境，`$LIB` 为 GAMESS 编译依赖库目录）：

```bash
$PY scripts/run_fmo_total_sensitivity.py \
  --gamess-dir software/gamess --library-dir "$LIB" \
  --dftb-parameter-dir software/gamess/auxdata/DFTB/3OB-3-1 --workers 8
$PY scripts/evaluate_fmo_total_sensitivity.py
```

本轮全量通过 `scripts/run_fmo_total_energy.py --workers 12` 在**独立的** `runs/fmo87/total_sensitivity/full/<protocol>` 路径运行；之后用 `scripts/run_fmo_total_sensitivity.py --all-cases --summarize-only --protocols pcm_polar_only no_waters cutoff_6 dftb_uff` 登记和核查覆盖度，最后 `scripts/evaluate_fmo_total_sensitivity.py --runs runs/fmo87/total_sensitivity/full` 重新审计所有成功日志。气相未达到试点准入门槛，不能从已经失败的气相子集选用结果。

### 87 例已有、完整且不拟合的复算（2026-09-27）

数据：`runs/fmo87/total_sensitivity/pilot_3/evaluation.json` → `baseline_87`。下表为**目标内** Pearson/Kendall 按目标样本数加权；已用标签观察了结果，不能以其中最好一项作为新的独立测试模型。

| 评分（已有 PCM 日志，无新 GAMESS） | Pearson | Kendall | 注意 |
|---|---:|---:|---|
| 完整 PCM 三项差 | -0.018 | 0.043 | 原预先指定主协议 |
| 仅 PCM 静电项的中间输出 | 0.047 | 0.102 | **不是**独立关闭非极性的新计算 |
| PCM 极化溶质，未加 `Gsol` | 0.227 | 0.209 | **不是**真正气相 |
| 完整 PCM + 固定 RDKit 旋转键惩罚 | -0.152 | -0.075 | 固定 `α=1` 反而变差 |
| PCM 静电项 + 固定 RDKit 旋转键惩罚 | -0.082 | -0.012 | 未拟合，但不是 PB/SA |
| 完整 PCM，游离配体单次 QM 优化 | -0.251 | -0.162 | 不是游离态构象系综 |

### 独立 12 例运行结果

运行记录在 `runs/fmo87/total_sensitivity/pilot_3/protocols.json`、`run_summary.json`、`evaluation.json` 及逐个 GAMESS 日志。若 12 例有失败，只计算各靶点**同成功子集**与原协议的配对比较；`3/3` 各靶点才允许四靶点都参与排序指标。每个样本的三项能量差保持 `AB−A−B`，不读取实验标签用于模拟或系数调整。

| 12 例变体 | 成功 | 变体 Pearson/Kendall | **相同成功样本**原 PCM Pearson/Kendall | 判读 |
|---|---:|---:|---:|---|
| PCM 仅静电独立重跑 | 12/12 | 0.336 / 0.333 | 0.110 / 0.333 | Pearson 上升；每靶点仅 3 例，全量见下 |
| 真正气相 | 2/12 | 不可评价 | 不可评价 | JNK1/TYK2 受体及 P38 三个复合物 SCF 未收敛 |
| 去掉显式水 | 12/12 | -0.233 / 0.167 | 0.110 / 0.333 | 双指标下降 |
| 共有口袋 6 Å | 12/12 | 0.092 / 0.333 | 0.110 / 0.333 | 无明确改善 |
| DFTB UFF 色散 | 12/12 | 0.389 / 0.167 | 0.110 / 0.333 | Pearson 上升、Kendall 下降；不能称改进 |

12 例的目标内 Pearson 每组仅 `n=3`，不稳定；所有候选均为已经看过标签的探索性结果。上述气相受体未收敛、其他未成功作业不参与相关性分母；不准把气相 `2/12` 误报为验证结果。

### 独立全量实算的四个协议对照

`runs/fmo87/total_sensitivity/full/pcm_polar_only/audit.json` 审计了 **87/87** 个与原协议只差 `ICAV=0 IDISP=0` 的独立 GAMESS FMO2 作业：四靶点样本数加权 Pearson/Kendall `0.0473/0.1017`；各靶点 CDK2 `-0.358/-0.217`、JNK1 `0.540/0.425`、P38 `-0.369/-0.167`、TYK2 `0.691/0.567`。逐样本与完整 PCM 原日志直接提取的 `EunD_es` 差不超过 `3.6×10⁻¹¹ kcal/mol`：这里非极性 PCM 项在该 GAMESS 协议中作为后处理加项，不改变本批电子分项。**它提升全队列平均相关性很有限，CDK2/P38 仍反序；并非 Tankyrase 论文中 PB/SA 的复现。**

`runs/fmo87/total_sensitivity/full/no_waters/audit.json` 独立审计 **87/87**：去水 Pearson/Kendall `0.0226/0.0628`，原始为 `-0.0177/0.0432`；主要只改变 P38（`-0.366/-0.206 → -0.263/-0.156`）。共享口袋的 atom count：CDK2 `535→535`、JNK1 `433→433`、P38 `635→578`、TYK2 `421→421`，其他三个靶点本来就没有被选中的结构水。与 12 例的负向变化不同，说明不能用小试点外推全量的符号；该全量结果改善也仍不足以达到可靠排序。

`runs/fmo87/total_sensitivity/full/dftb_uff/audit.json` 独立审计 **87/87**：GAMESS DFTB `DISP=UFF` 的 Pearson/Kendall 为 `0.0281/0.0998`，各靶点 CDK2 `-0.528/-0.350`、JNK1 `0.425/0.368`、P38 `-0.370/-0.160`、TYK2 `0.909/0.750`。TYK2 上很强不代表整体改善，CDK2/P38 依然反序；DFTB UFF 的修正与 PCM 色散 `IDISP` **不是同一个项**。

`runs/fmo87/total_sensitivity/full/cutoff_6/audit.json` 独立审计 **87/87**：共有口袋从 5 Å 扩到 6 Å 后 Pearson/Kendall `-0.0470/0.0347`，各靶点 CDK2 `-0.126/-0.100`、JNK1 `-0.026/0.100`、P38 `-0.361/-0.185`、TYK2 `0.671/0.550`，没有提升整体排序。

### 全队列归纳

| 固定几何三项差协议 | 成功/队列 | 目标样本数加权 Pearson | Kendall | 与原 PCM 相比 |
|---|---:|---:|---:|---|
| 原始 PCM＋结构水、5 Å、DFTB3 | 87/87 | -0.0177 | 0.0432 | 原始预先冻结 |
| PCM 仅静电独立 deck | 87/87 | 0.0473 | 0.1017 | 两项略升，但仍很弱 |
| PCM、去结构水 | 87/87 | 0.0226 | 0.0628 | 两项略升，但仍很弱 |
| PCM、DFTB UFF 色散 | 87/87 | 0.0281 | 0.0998 | 两项略升，但仍很弱 |
| PCM、共有口袋 6 Å | 87/87 | -0.0470 | 0.0347 | 两项下降 |
| 真正气相（试点） | 2/12 | 不可评价 | 不可评价 | 计算不稳定，非完整对照 |
| 原始 PCM＋单次放松游离配体 | 87/87 | -0.2508 | -0.1618 | 游离态处理后反而恶化 |

注意：所有指标都是**同一批 87 个已观察实验标签**上的探索性观察；不能对 4 个新变体再选一次最好的分数并宣称获得经验证的零拟合新方法。所有自带 kcal/mol 单位的三项差在这里仅当排序代理，不是实验相对结合自由能的无偏数值估计。

## 5. 未完成与升级门槛

本轮不声称跑过未落盘的 GAMESS；未完成的 PB/SA、ONIOM、蛋白/复合物 QM 几何松弛、微态/多构象、水重排、熵自由能需要独立工作流。只在本 87 队列里继续调整评分权重会造成数据泄漏。任何完整 87 候选均需相同成功样本、输入身份、成功率、对照指标与**新靶点/新化学系列外部验证**；计算是否可行与预测是否泛化是两个问题。
