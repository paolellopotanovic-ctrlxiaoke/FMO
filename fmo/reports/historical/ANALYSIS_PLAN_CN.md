# 完整分析与改进计划

## 目标与约束

目标不是把任意后处理分数调到论文数字，而是判断差距来自哪一层，并在每一层给出可复现、可失败、可停止的证据：

1. 先确认数据集、指标和 baseline 没有实现错误；
2. 再隔离 QM deck、溶剂/水、片段、几何、pose 和 scorer 定义；
3. 只把完整 87/87 且审计通过的协议与论文 Table 12 做主比较；
4. 对不完整协议只报告运行率和同成功子集 paired ranking，不把选择性成功当作性能提升。

## 阶段 A：基线与论文合同

**输入**：`data/fmo87_cohort.csv`、Boltz2 manifest、GAMESS result、SophosQM 主文和 SI。

**检查**：

- 87 个身份、四靶点计数、neutral formal charge、SMILES/hash；
- target-average Pearson/Kendall 的实现；
- tPIE 只累计 ligand-environment pair total；
- tPIE 与 selected final TIE reconciliation；
- Boltz2 的输入 identity、affinity 转换和目标权重。

**停止标准**：若 baseline 审计失败，先修审计，禁止做物理协议实验。

## 阶段 B：已知 deck 变量

按单变量/小组合顺序测试：

1. PCM：开启/关闭；
2. DFTB dispersion：NONE/UFF；
3. GAMESS FMO DFTB 模式：默认值与 `MODESD=2/MODGAM=1`；
4. 结构水：保留/去除；
5. pocket cutoff：4/5/6 Å。

**实验设计**：每个协议先跑每靶点 3 个均匀抽样 case；只有成功率接近 baseline、paired ranking 不恶化，才进入全量。

**已完成结论**：

- PCM + UFF 默认模式全量仅 79/87，成功子集 0.2068/0.3324，不采用；
- gas 协议大面积 SCF failure；
- `MODESD=2/MODGAM=1` 触发当前 build 的 `lumpesd_` segmentation fault；
- 去水、4 Å、6 Å 在 12-case 上都不优于同成功子集 baseline。

## 阶段 C：几何协议

比较可获得的公开替代几何：

- 生产 flexible receptor/water + 100 kJ/mol/Å² ligand restraint；
- fixed prepared receptor/water；
- 10 kJ/mol/Å² 弱 restraint；
- 1000 kJ/mol/Å² 强 restraint。

记录 ligand RMSD、最大位移、能量下降、成功率和 paired ranking。

**已完成部分结论**：fixed receptor 12/12 成功但 paired Pearson/Kendall 从 0.4349/0.5000 降至 -0.0084/0；weak/strong restraint 仅 9/12，不能升级。

## 阶段 D：片段与边界规则

这是当前最高优先级的下一阶段，因为它直接对应 SophosQM/FACIO 描述且不依赖私有软件：

- shifted Cα–C′ cut vs 未 shift 的严格对照；
- HOP_C 开启/关闭或不同 HOP 类型；
- 端点 ACE/NME 与简化 cap；
- pocket residue 连通性 bridging 开启/关闭；
- 水作为独立 fragment 与水并入环境的对照；
- fragment charge/multiplicity 及 protein protonation audit。

每个变体必须保持相同坐标，先用单个 CDK2/JNK1/P38/TYK2 case 验证 deck 物理性，再跑 12-case。

**已完成边界检查**：关闭 HOP 导致 GLU 片段 57 个电子的奇电子 singlet，12/12 在 deck 电子校验失败；关闭 sequence-gap bridging 的初次成功样本出现 126.777 kcal/mol ligand pair。pair energy 质量门已从 1000 收紧到 100 kcal/mol，并以 fresh rerun 验证该协议 0/12；生产 baseline 最大 pair 为 28.019 kcal/mol。

## 阶段 E：pose 与局部结构

公开输入无法证明等价于作者内部 Glide pose，因此下一优先级是 pose/局部几何，而非继续调 PCM：

- 比较 source/reference/fmp 坐标是否一致；
- 对存在多 pose 的记录计算 ligand heavy-atom RMSD、接触残基变化和 tPIE 差异；
- 对每个 target 选择少量高/低 tPIE 与大残差 case 做局部 pose 扰动；
- 只在预先声明的 pose 集合中选择最低 FMO energy，不能用实验标签选 pose。

## 阶段 F：PIEDA 与校准诊断

已完成全量分量分析：ES、EX、CT、DI、solvation 的简单替换都没有达到论文 FMO 参考；`tPIE - EX` 仅为 0.4932/0.3236。因此不再把“遗漏某一个 PIEDA 分量”作为主假设。

SophosQM `α·tPIE + β·clogP + γ` 只作为 secondary model：in-sample 0.6088/0.4996，LOO 0.4020/0.3309。除非有独立训练/验证 series，不把系数拟合提升称作私有 FMO 复现。

## 阶段 G：全量决策

候选进入全量的条件：

- 12-case 成功率不低于 baseline 的可接受范围；
- 四靶点没有明显单靶点崩溃；
- paired ranking 至少不比 baseline 明显恶化；
- manifest、deck hash、GAMESS hash、失败日志齐全。

候选升级为生产协议的条件：

- 严格 87/87 成功；
- 统一 evaluator 输出完整 target-average 指标；
- 与 baseline 的提升不由失败样本选择造成；
- 对关键靶点的改善有 bootstrap/paired uncertainty 支持，或至少有清晰机制解释。

## 当前下一步

1. 不再运行已证伪的 gas、MODESD=2 和 PCM+UFF 全量变体；
2. 补齐片段/边界规则对照；
3. 检查公开结构的 pose provenance 和多 pose 情况；
4. 若片段/pose 仍无法缩小差距，则结论应是“公开 GAMESS baseline 与 Boltz 私有 scorer 的不可辨识协议差异”，而不是继续无界搜索。

## Pose 实验诊断（2026-09-26）

公开 JNK1 SDF 确实包含 `flip` alternate pose；对 18636-1、18629-1、18660-1 的初步运行显示，直接用翻转坐标重新选 pocket 会同时改变 pocket residue 和 `NFRAG`（例如 18636-1 从 25 变为 26），因此不能把得到的 tPIE 变化解释为纯 pose 效应。进一步尝试固定 baseline pocket 后，现有 `prepare_job` 的自动 fragment-bond 一致性校验拒绝了该输入。下一步若继续研究 pose，必须先增加“固定 protein atom/fragment membership、仅替换 ligand coordinates”的显式 API；在此之前不进行 21 个 JNK1 pose 的全量扩展。

Boltz-2 的 JNK1 输入 YAML 只有 protein sequence 和 canonical ligand SMILES，没有 non-flip/flip 坐标；它自行生成 pose。将 Boltz 输出与公开 protein 坐标对齐后，18636-1、18629-1、18660-1 均更接近 non-flip（分别为 `0.714 vs 2.032 Å`、`0.785 vs 1.767 Å`、`1.248 vs 2.673 Å` 的 ligand heavy-atom RMSD）。因此当前最准确的表述是“Boltz 输出 pose 更接近 non-flip”，而不是“Boltz 明确选择了 non-flip 输入 pose”。
