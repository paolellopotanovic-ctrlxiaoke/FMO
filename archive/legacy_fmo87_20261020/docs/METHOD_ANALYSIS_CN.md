# SophosQM 方法复现分析

## 论文中实际做的事情

来源：Guareschi et al., *SophosQM: Accurate Binding Affinity Prediction in Compound Optimization*, ACS Omega 2023, 8(17), 15083–15098, DOI `10.1021/acsomega.2c08132`, PMID `37151542`。

论文把问题拆成两个近似项，而不是跑显式分子动力学或严格计算结合自由能：

1. 以 FMO2 下配体与口袋各片段的相互作用和 `tPIE` 近似焓项；
2. 以配体 `clogP` 作为非焓项/熵效应的宏观代理；
3. 每个相似结合模式的 series 拟合 `ΔG_sim = α·tPIE + β·clogP + γ`；
4. 论文对每个化学 series 用实验值拟合 `α/β/γ`，然后报告同一批化合物的 fitted prediction 与实验值相关性/误差，并强调固定几何、同系列/相似结合模式的适用边界。

FMO 的计算是片段嵌入场下的量子化学：先求片段单体，再求片段二聚体并以 pair correction 构造二体展开；配体的 `tPIE` 是 ligand–environment pair interactions 之和，不是复合物总能量，也不是把所有 FMO dimer energies 全部求和。GAMESS 负责底层量子化学和 FMO/PIEDA 求解；我们的代码应实现的是可审计的计算协议、输入/输出编排与 benchmark，而不是重写积分/SCF/MP2 数值引擎。

论文给出的实现级参数包括：GAMESS、FMO2、MP2、6-31G*；高通量 SQM 轨道另用 DFTB3/3ob。主文的 70 个化合物覆盖 DNA ligase、MCL-1、MUP-1、JAK-2、HSP90、p38 六个靶点，MCL-1 系列有形式电荷 −1，其余所述系列为中性。主流程结构处理使用 Schrödinger Protein Preparation、保留配体 5 Å 内水并优化氢键网络；Glide SP 对接、与晶体配体核心对齐、OPLS3e 约束最小化；截取配体约 5 Å 内残基，以 ACE/NME 封端。FACIO 按连续主链的 Cα–C 键切割构造片段，并用 bond-detachment hybrid orbitals 处理切口；配体与每个保留水各自为一个片段。不能把每个完整残基直接作为独立片段并省略切口处理。`clogP` 使用 VolSurf。

这些是论文的方法，不代表原作者实现源码已公开。官方 Supplementary Information 已取得并缓存于 `data/references/sophosqm_si.pdf`，其中列出残基、结构、配体、实验值、拟合参数和 pose RMSD 等复现关键材料。Schrödinger/Glide、OPLS3e、VolSurf 是论文工作流组件，不可默认免费或开源；本项目的替代组件均已显式记录，不能声称逐项复现。

### 原文报告口径

SophosQM 正文使用的公式是 `ΔG_sim = α·tPIE + β·clogP + γ`，其中 `α/β/γ` 由一组实验亲和值拟合得到；SI 对每个 series 列出这些系数和对应的 FMO-predicted affinity。正文和 SI 均未出现 leave-one-out、cross-validation、external test set、held-out 或 out-of-sample 等验证表述。因此其报告的实验/预测相关性和误差应按 in-sample fitted performance 解读。原文讨论了拟合参数可靠后可外推到未知化合物的适用条件，但未把 LOO 作为论文中的报告指标。

## 本项目实现的 baseline

主 benchmark 不用 MP2 替代 Boltz-2 论文中的 DFTB3 baseline，而是调用 GAMESS FMO2-DFTB3/3OB-3-1/PIEDA。自行维护的是结构制备与原子映射、FACIO-like fragment/cap/water 规则、GAMESS 输入与运行、真实输出解析、tPIE 聚合、clogP、校准消融和 benchmark。不得把 GAMESS 本身称为我们实现的 FMO 数值算法，也不能把输出 tPIE 直接称作热力学 `ΔG`。

对配体片段 L，建议首个严格定义为 `tPIE(L) = Σ_r PIE(L,r)`，其中 r 为预先声明口袋中的蛋白残基、水片段。只累计唯一无序片段对的 `total` 列；PIEDA 的 ES/EX/CT/DI/solvent 分量另做解释性分析，不能把分量行重复相加。需用 GAMESS 真实输出验证列单位、片段 ID 映射、配体端方向、总量与组件和的一致性。

## Ross2023 适配不是原论文复现

SophosQM 原论文用实验绝对 ΔG 拟合每个靶点/化学系列的三个系数；Ross2023 公共表主要给 pairwise `ΔΔG` edge。对 Ross 做 `ΔΔG_pred = α·ΔtPIE + β·ΔclogP` 是合理的 benchmark 适配，但属于新的 pairwise 训练目标，不能冒称逐步复刻原论文 absolute regression。必须先固定 edge 符号约定，保证 `ΔΔG(i,j)` 与结构化数据列方向一致。训练/测试应按 ligand 节点或 scaffold 划分，测试配体的任何关联标签不得进入拟合。

Boltz 论文明确说其 4-target affinity subset 是 CDK2、TYK2、JNK1、P38 上的 87 个 neutral compounds，来源为 Hahn protein-ligand-benchmark；其较便宜的 FMO baseline 是作者的 in-house DFTB3 FMO code，输入是 Glide docked pose，FMO 平均 target Pearson R=0.55、Kendall τ=0.38。原论文附录未给出可直接复现这套 FMO 实现的源码，因此不能把 Boltz 论文里的 FMO 代码当可复用开源后端。本地 Ross/Chen 四靶点 SDF 有 104 条记录；按 canonical isomeric SMILES 去重后，17 条 JNK1 alternate-flip pose 是重复化合物，保留 non-flip 代表得到精确 16+21+34+16=87。

所以本项目明确区分两项实验：A) 论文式 SophosQM FMO2/MP2/6-31G* 校准思想与官方 SI 系数消融；B) 主 benchmark 的 GAMESS FMO2/DFTB3/3ob raw tPIE，使用同一 87 ligand/pose 并只比较排序指标。A 的拟合回归指标不得冒充 B 的 raw-score `Pearson=0.55 / Kendall=0.38` 对照，也不得在测试集上拟合系数后再报告训练相关性。

## 公开性与复用边界

- **可复用**：GAMESS 的 FMO/MP2/基组计算能力（按 GAMESS 注册与许可条件获取使用）；Ross2023 已公开数据；Boltz 已发布的代码/权重按各自许可使用；常规结构/分子工具需按其许可证检查。
- **不能假设可复用**：SophosQM 原作者 FACIO/GAMESS 生成脚本、输出分析与拟合工作流；Schrödinger/Glide、OPLS3e、VolSurf 许可组件；Boltz 论文的 FEP+/四靶点 benchmark harness、作者侧 prediction manifest/Table 12 生成脚本；Boltz 作者内部 FMO 对照执行代码/全部样本映射。
- **OpenFF 边界**：SophosQM 理论本身不要求 OpenFF；本项目的标准几何替代管线使用 OpenFF 2.2.0/GNN AM1-BCC。RDKit Crippen MolLogP 只是 VolSurf clogP 的开源替代描述符，不能声称数值等价。
- **不需要自写量子化学引擎**：本项目自行实现 baseline orchestration 与 scoring protocol；QM/FMO 数值核交由 GAMESS。

## 当前实现边界

FMO87 主线已完成：精确 87 化合物身份重建、PDB2PQR/PROPKA 制备、OpenMM 标准几何、shifted Cα–C′ fragment、ACE/NME/OXT 修复、`HOP_C`、闭壳层电子校验、GAMESS FMO2-DFTB3/PIEDA、TIE 对账、Boltz2 复现、Table 12 式指标、10000 次 bootstrap 和两个 87/87 显式水消融均已通过。Glide SP pose、MacroModel/OPLS3e 和 VolSurf clogP 不可得，均已由公开替代品和明确偏差声明处理。不能把该项目称为 Boltz 作者私有 FMO scorer 的复现。
