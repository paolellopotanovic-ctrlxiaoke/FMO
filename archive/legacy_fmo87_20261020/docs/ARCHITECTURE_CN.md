# 项目架构

## 目标

实现一套独立、可审计的蛋白–配体 FMO scoring 工具链：底层电子结构由 GAMESS 求解；本项目负责队列冻结、蛋白制备编排、标准几何、FMO 分片与封端、输入生成、作业调度、PIE 解析、tPIE 对账、Boltz-2 复现实验和指标评估。

## 数据流

```text
public FEP benchmark structures
    ↓
canonical-isomeric-SMILES identity and duplicate collapse
    ↓ frozen fmo87 cohort
    ├── FMO branch
    │   ├── PDB2PQR/AMBER/PROPKA pH 7.4
    │   ├── OpenMM AMBER14/OpenFF restrained geometry
    │   ├── 5 Å pocket and shifted Cα–C′ fragments
    │   ├── ACE/NME ends + internal FMOBND HOP_C
    │   ├── GAMESS FMO2-DFTB3/3OB-3-1/PIEDA
    │   └── parser + selected-TIE reconciliation
    └── Boltz-2 branch
        ├── canonical SMILES + shared target MSA
        ├── Boltz 2.2.1 affinity prediction
        └── checkpoint/input/JSON provenance audit
    ↓
exact identity merge
    ↓
per-target metrics, compound-count weighting
    ↓
bootstrap + SophosQM calibration ablation
    ↓
explicit-water scorer/preparation ablations
```

## 模块职责

- `cohort.py` / `build_fmo87_cohort.py`：重建并冻结 87 化合物身份。
- `preparation.py`：过滤配体、保留或制备前移除水、调用 PDB2PQR/PROPKA，并保存 provenance。
- `minimization.py`：OpenMM AMBER14/OpenFF flexible complex 与 ligand positional restraint。
- `structure.py`：PDB/SDF 坐标和残基/配体结构表示。
- `fragmentation.py`：5 Å pocket、shifted Cα–C′ fragment、完整 ACE/NME、水 fragment、HOP 边界和电子守恒审计。
- `fmo_input.py` / `gamess.py` / `jobs.py`：生成 GAMESS FMO2-DFTB3 deck、manifest 和可执行环境。
- `runner.py` / `parser.py`：运行 GAMESS、解析真实 PIEDA pair 表、聚合 tPIE，并强制与 selected final TIE 对账。
- `boltz_inputs.py`：生成 Boltz YAML、共享 MSA manifest 和输入哈希。
- `exact_evaluation.py`：严格身份合并、Table 12 式逐靶点加权指标和分层 bootstrap。
- `scripts/summarize_water_ablation.py`：汇总 retained/scorer/preparation 三种水策略的全队列指标和 P38 成对位移。

## FMO 计算边界

本项目不重写 SCF、DFTB、基组积分或 GAMESS 本体。复用 GAMESS 官方求解器；自行实现的研究部分是：

1. 可审计的队列、结构、质子化和几何处理；
2. FACIO-like fragment、ACE/NME、HOP 和电子计数协议；
3. GAMESS FMO 输入和作业调度；
4. 原始 PIEDA pair table 解析与 tPIE 聚合；
5. 与 Boltz-2 论文一致的 exact-cohort 评估；
6. bootstrap、拟合消融、显式水消融和失败关闭审计。

## 不允许的隐式行为

- GAMESS 不可执行时不能自动换成点电荷或 LJ 能量并标记为 FMO。
- 解析不到 PIE 时不能用总能量差伪造 PIE。
- GAMESS 正常退出但存在未收敛 monomer/dimer、geometry error 或 TIE 对账失败时，不能标记成功。
- 测试标签参与拟合时不能称为 held-out 结果。
- Boltz `affinity_pred_value`、binary probability 和 FMO tPIE 必须分列保存。

## 当前验证边界

- `aggregate_tpie` 将单条片段对的 `total` 列计一次，只保留恰好一个端点为配体的片段对；不会重复相加 PIEDA 分量。
- 87/87 生产 GAMESS 输出通过 normal termination、无 SCF/geometry/electron failure、pair table 存在、TIE reconciliation、输入哈希和化学协议审计。
- tPIE 是构象、质子化、水、片段与溶剂协议依赖的 enthalpic score，不等于严格结合自由能。
- 主比较是 raw tPIE 对 Boltz-2 论文内部 FMO；拟合 SophosQM-style 分数只作为次要分析。
