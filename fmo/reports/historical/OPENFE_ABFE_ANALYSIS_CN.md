# Boltz-2 中 OpenFE 与 ABFE 的论文、数据集和复现边界

日期：2026-09-27

## 结论先行

Boltz-2 中的 OpenFE 和 ABFE 不是同一个实验：

| 项目 | Boltz-2 中的用途 | 数据集 | 是否进入 Table 12 四靶点比较 |
|---|---|---|---|
| OpenFE | 公开 FEP+ benchmark 上的相对 FEP/RBFE baseline | 876 个 hit-to-lead / lead-optimization measurements；另有四靶点 87-case 汇总 | 是；Table 11 是 876-case，Table 12 有四靶点 OpenFE 行 |
| ABFE / Boltz-ABFE2 | TYK2 virtual-screening 后的绝对 FEP 验证；另有四靶点方法学 benchmark | TYK2 生成/商业库候选；以及 TYK2、P38、JNK1、CDK2 benchmark | 是；Table 12 有四靶点 ABFE 行，Boltz-2 正文的应用验证主要是 TYK2 |
| FMO | 便宜的结构打分 baseline | 4-target subset，87 个 neutral compounds | 是；对应 Table 12 |

因此，不能把“OpenFE 876-case 结果”“四靶点 87-case OpenFE/ABFE 结果”“Boltz-ABFE TYK2 筛选结果”和“Table 12 87-case FMO 结果”混为一个复现实验。

## 1. Boltz-2 对两个方法的引用和使用

Boltz-2 论文把 OpenFE 描述为公开的 relative FEP 方法，并在完整 OpenFE subset 上比较；论文正文/附录称该集合有 876 个高质量 hit-to-lead measurements。论文 Table 11 的 OpenFE 结果为 target-average Pearson `0.63`、Kendall `0.47`；Boltz-2 为 `0.62`、`0.46`。此外，Table 12 的四靶点 87-case 汇总也包含 OpenFE 行（Pearson `0.66`、Kendall `0.51`），同时比较 FEP+、ABFE、FMO 和 MM/PBSA。

ABFE 在 Boltz-2 中有两种语境：

1. **Table 12 physics baseline**：四靶点 benchmark 上的 ABFE 结果，论文只提供汇总指标。
2. **Boltz-ABFE2 prospective/structure-prediction workflow**：先用 Boltz-2 生成 protein–ligand complex，再用作者的绝对 FEP protocol 计算候选的 ABFE；主要用于 TYK2 virtual screening。论文称该流程在 protein-ligand benchmark 的 TYK2 上 Pearson `0.95`、centered MAE `0.42 kcal/mol`、偏移约 `0.92 kcal/mol`，并据此验证筛选出的候选。

因此，ABFE 论文/方法应按“绝对 FEP protocol”理解；它不是 OpenFE RBFE 的另一个命令行模式。Boltz-2 正文明确将 Wu et al. 的 ABFE protocol 用作四靶点 baseline，并在 TYK2 预测复合物上用于 prospective screening 验证。

## 2. OpenFE：方法和数据集

### 方法

OpenFE 的核心是 relative binding free energy（RBFE）：对同一 congeneric series 的 ligand A/B，分别计算 complex leg 和 solvent leg 的 alchemical transformation，再取 thermodynamic cycle 的差值：

```text
ΔΔG_binding(A→B)
  = ΔG_complex(A→B) - ΔG_solvent(A→B)
```

标准 OpenFE 工作流包括：

1. 读取 protein PDB 与 ligand SDF；
2. 生成 ligand partial charges；
3. 用 kartograf/Lomap 生成 atom mapping 和 ligand network；
4. 为每条边运行 complex 与 solvent 两条 alchemical legs；
5. 用 OpenMM/OpenMMTools 采样；
6. 用 MBAR/BAR 等方法分析；
7. 用 Cinnabar/`openfe gather` 汇总节点自由能、ΔΔG 和不确定度。

公开 OpenFE tutorial 的默认 RBFE 例子使用水 + 0.15 M NaCl、约 1.2 nm buffer、kartograf mapping、LOMAP 网络和 OpenFE 的 OpenMM Hybrid Topology RFE protocol。实际论文 benchmark 的确切参数必须以对应 benchmark/config 为准，不能把当前版本默认值自动当作 Boltz-2 作者运行值。

### Boltz-2 Table 11 的数据集

Boltz-2 论文称 OpenFE subset 来自 Ross/FEP+ benchmark，包含 876 个 hit-to-lead/lead-optimization measurements；这些不是当前项目的 87 个四靶点 compound 的同一数据表。OpenFE 的 876 记录更接近相对 FEP 的 transformations/measurements，而非简单的 876 个独立蛋白–配体绝对 ΔG。

当前仓库包含 Ross public benchmark 的公开结构输入和 FEP+ 结果目录，但没有一份经过 Boltz-2 作者最终筛选、命名和聚合的 876-case OpenFE 输入 manifest；在当前可见的 Boltz 源码/数据树中也没有找到逐点 Table 11 输出包。因此：

- 可以获取公开 Ross/FEP+ 结构和部分公开结果；
- 可以使用 OpenFE 自己搭建 RBFE network；
- 不能仅凭当前仓库确认 Boltz-2 Table 11 的 876 条逐点输入、边、过滤条件和最终结果完全一致；
- 不能把在当前 87-case 上运行 OpenFE 的结果称为论文 Table 11 复现。

## 3. ABFE：两篇论文实际做了什么

### Wu et al.：Optimizing absolute binding free energy calculations for production usage

这篇工作比较了两套 ABFE protocol：

- Aldeghi-like protocol；
- RXRX protocol。

共同的物理循环是：

1. 在水中把 ligand alchemically annihilate 到 non-interacting dummy；
2. 在 complex 中施加 Boresch-style restraints；
3. 在 complex 中 annihilate ligand；
4. 加入 restraint/free-energy standard-state corrections；
5. 取 complex 与 solvent legs 的差得到 ΔGbind。

RXRX 的主要改进不是换成另一种自由能定义，而是：

- 自动选择更稳定的 ligand/protein anchor atoms；
- 优化 Boresch distance/angle/dihedral restraint；
- 调整非键相互作用 annihilation 的 λ-schedule；
- 先处理 dihedral restraints，再处理 distance/angle restraints，以减少数值不稳定；
- 每个 ligand 做三次 repeat，报告均值和重复间不确定度。

论文在 TYK2、P38、JNK1、CDK2 四个 protein-ligand benchmark systems 上比较 RXRX 和 Aldeghi；结果显示 RXRX 的 RMSE、MUE 和 repeat variance 普遍更低。该论文不是 OpenFE 的实现论文，而是 Recursion/BioSimSpace/ArgoWorkflows 生产平台上的 ABFE protocol 优化。

### Boltz-ABFE：Free Energy Perturbation without Crystal Structures

Boltz-ABFE 把问题拆成：

```text
Boltz-2 sequence/SMILES co-folding
    ↓
predicted protein–ligand complex
    ↓
automated structure improvement / restraint setup
    ↓
absolute FEP (ABFE)
```

它的科学问题不是重新比较 Boltz-2 和 FMO，而是检查：没有 experimental crystal structure 时，Boltz-2 生成的 complex 是否足以作为 ABFE 初始结构。

Boltz-2 附录对 TYK2 screening 的实际流程给出了这些关键步骤：

- 先用 Boltz-2 对候选做 affinity/likelihood screening；
- 对入选 ligand 在 pH 7.4 选择 dominant tautomer；
- 再用 Boltz-2 co-fold TYK2–ligand complex；
- 调用 Boltz-ABFE protocol 计算 absolute ΔG；
- 多 stereoisomer 时选择 ABFE 更高者；
- HLL、Kinase Library、SynFlowNet、random 和 public TYK2 binder 各自筛选，再提交候选做 ABFE。

TYK2 virtual screen 的最终候选主要是 10 个一组的筛选结果，不是 87-case 的逐点 Table 12 预测。论文称 HLL top 10 中 8 个、Kinase Library top 10 中 10 个被 ABFE 预测为 binder，random controls 则为 non-binder；这是 prospective ranking/triage 叙事，不是四靶点全量 ABFE benchmark。

## 4. 开源状态

| 组件 | 开源状态 | 可否直接复现 |
|---|---|---|
| OpenFE 核心代码 | MIT，公开 GitHub/文档 | 可以，需准备环境、结构、ligand network 和算力 |
| OpenFE generic RBFE protocol | 公开 | 可以做协议级复现；结果受版本、force field、mapping、sampling 影响 |
| OpenFE generic ABFE protocol | 公开 API/文档，但部分 `gather-abfe` 能力仍标为 experimental | 可以做 exploratory protocol-level ABFE，不等于 RXRX/Boltz-ABFE |
| Ross public FEP+ benchmark structures | 公开 | 可以获得结构、ligand SDF、部分实验/edge 数据和 FEP+结果；许可需逐项遵守 |
| Boltz-2 code/weights | 公开 | 可以复现 Boltz-2 结构/affinity 输入部分 |
| Boltz-2 Table 11 876-case exact manifest/results | 当前未见作者提供完整逐点复现包 | 不能保证 exact reproduction |
| Wu/RXRX production implementation | 论文公开、方法描述公开；实现依赖 Recursion BioSimSpace/Argo production platform | 不能直接复现原生产 pipeline |
| Boltz-ABFE2 end-to-end implementation | 当前没有在本仓库或 Boltz 公开代码中发现等价实现 | 不能直接 exact reproduce |
| OpenMM/BioSimSpace/相关基础组件 | 多数公开，但版本和许可不同 | 可组装近似实现 |

开源边界的关键区别是：**“论文开放访问”不等于“生产代码、输入 manifest、checkpoint、workflow 和结果全开源”。** OpenFE 属于真正可运行的开放工具链；RXRX/Boltz-ABFE 是论文级方法描述加内部生产平台，当前更适合做 protocol-inspired reproduction。

## 5. 对本项目的复现路线

### Stage 0：冻结身份和输入

- 保留当前 87-case cohort，仅作为 focused four-target benchmark；
- 使用 Ross/FEP+ 原始 ligand SDF、protein PDB、实验 `r_exp_dg`；
- 记录 OpenFE 版本、OpenMM/OpenMMTools/OpenFF/RDKit 版本；
- 不把当前 FMO 的 PDB2PQR、FACIO fragment 或 GAMESS deck 传给 OpenFE，除非实验明确研究结构制备影响。

### Stage 1：OpenFE smoke test

先选每个 target 2–3 个 congeneric ligand pair：

- 验证 ligand charge、atom mapping、complex/solvent system 建立；
- 生成两个 legs；
- 跑短 pilot，检查 energy overlap、forward/reverse agreement、cycle consistency；
- 只输出 ΔΔG 和 uncertainty，不与绝对 FMO tPIE 直接比较。

### Stage 2：87-case OpenFE focused benchmark

若 smoke test 通过，再对四个 target 的 87 compounds：

- 按原始 congeneric series 建 network，而不是把 87 个 ligand 随意连成一条链；
- 每条 transformation 同时跑 complex/solvent legs；
- 默认至少 3 repeats；
- 用 cycle closure、MBAR uncertainty、edge overlap 和 network topology 做质量门；
- 用节点 ΔG 或 edge ΔΔG 与实验数据评价；
- 明确这是“87-case OpenFE focused reproduction”，不是论文 876-case Table 11 reproduction。

### Stage 3：完整 876-case OpenFE 追踪

这一步必须先解决数据问题：

1. 获取 Boltz-2 作者最终 876-case manifest；或
2. 从公开 Ross/FEP+ benchmark 重建并记录自己的过滤规则；
3. 对比 ligand count、target/series count、edge count、experimental labels；
4. 只有 identity、edges 和 assay aggregation 都一致，才可以称为 Table 11 reproduction。

### Stage 4：generic OpenFE ABFE pilot

这不是 Boltz-ABFE exact reproduction，而是验证环境和绝对 FEP 基础能力：

- 先选 1 个 TYK2 ligand；
- 用公开 protein/ligand pose；
- 跑 OpenFE `AbsoluteBindingProtocol` 或 SepTop；
- 检查 complex/solvent legs、Boresch restraints、standard-state correction；
- 3 repeats 后评估 variance；
- 再扩展到 2–3 个 TYK2 compounds。

### Stage 5：Boltz-inspired ABFE

若目标是接近 Boltz-ABFE：

- 使用 Boltz-2 对 sequence + canonical SMILES co-fold；
- 从 Boltz CIF 提取 ligand/protein complex；
- 做 protein/ligand topology 和 protonation/tautomer audit；
- 实现自动 restraint anchor selection；
- 采用 RXRX 思路的 λ-schedule 和 staged restraint；
- 每个 compound 运行三次；
- 与公开 TYK2/P38/JNK1/CDK2 benchmark 的实验 ΔG 比较。

但最终报告必须命名为：

```text
Boltz-inspired OpenMM/OpenFE ABFE reproduction
```

而不是：

```text
exact reproduction of Recursion Boltz-ABFE2 / RXRX
```

## 6. 与当前 FMO 项目的关系

OpenFE/RBFE 和 ABFE 都比当前 FMO 更接近热力学 ΔG：

- FMO tPIE 是 ligand–environment interaction/enthalpic ranking score；
- RBFE 直接估计同一 series 的 ΔΔG；
- ABFE 直接估计 standard binding free energy ΔG，但对 pose、restraints、sampling、protonation 和 standard-state correction 极敏感。

因此不能用 OpenFE/ABFE 的结果简单“校准”FMO tPIE，再声称 FMO 精度提升。应将三者作为独立物理 baseline：

```text
Boltz-2 ML affinity
OpenFE RBFE ΔΔG
ABFE ΔGbind
GAMESS FMO tPIE
```

共同使用同一实验标签时，必须保持 target/series identity、符号、单位和 centered/non-centered 口径一致。

## 7. 当前建议

1. **先做 OpenFE smoke test**：先在目标环境冻结 OpenFE、OpenMM、OpenFF、RDKit 版本，不要直接启动 87-case 长跑。
2. **优先做 87-case focused RBFE**：数据和结构在本地，科学问题清晰，且 OpenFE 是真正开源工具链。
3. **ABFE 先做 1–3 个 TYK2 pilot**：验证 topology、restraints、standard-state 和 repeat variance。
4. **不要宣称复现 Boltz-2 的 876-case OpenFE 数字**，除非拿到 exact 876 manifest。
5. **不要宣称复现 RXRX/Boltz-ABFE2**，除非获得其生产 workflow、参数和输入；当前只能做公开方法启发的近似复现。
