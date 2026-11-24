# 亲和力论文的参数拟合与“新几何”验证：文献对照与 FMO87 方案

日期：2026-09-27。这里的“预测新化合物”指 **测试化合物的实验标签从未参与参数、制备规则或姿势选择**；仅仅对已有样本换一套几何重算，不等于测试新化合物。本文是实验设计，不声称已运行新的 QM 实验。

## 1. 文献如何得到最终结果

| 文章/类别 | 参数如何得到 | 测试实际隔离了什么 | 本项目可借鉴与限制 |
|---|---|---|---|
| SophosQM，Guareschi 等（2023，DOI `10.1021/acsomega.2c08132`） | 六个靶点各用同系列**全部已知实验 ΔG** 拟合 `ΔG=α×tPIE+β×clogP+γ`，共六套系数 | 主文 `R²≈0.88` / MAE `≈0.4` 为拟合内值；没有冻结系数后独立新化合物的主指标 | 证明了同系列校准可拟合，不证明零标签新靶点或新 pose 性能；见 `SOPHOSQM_CLAIMS_AUDIT_CN.md` |
| 自动化 LIE，van Dijk 等（2017，DOI `10.1021/acs.jcim.7b00222`） | 对一批 aromatase 抑制剂，用实验亲和力拟合 MD 中静电、范德华特征的 `α,β`（可加 `γ`）；通过模型选择得到多套子模型 | 报告留一和 bootstrap 交叉验证；选择可解释数据子集仍是回顾性，不能把覆盖子集当全部化合物外部预测 | 同样要把“拟合哪些系列、覆盖率、剔除/分簇规则”写清；不能只报最佳子模型 |
| FXR/LIE D3R 盲测，Rifai 等（2017，DOI `10.1007/s10822-017-0055-0`） | 用先前公开的同靶点已测 ligand 训练不同化学类的 LIE 局部系数，再对 D3R 提交计算结果 | D3R 实验结果公布后才对提交的**新 ligand**评价，同时按训练化学/相互作用适用域分层；最高信心仅 9 个 ligand、标准预测误差约 `5 kJ/mol`，其他类有明显失败 | 最接近所需“锚点校准→冻结→新 ligand 新 pose→盲测”；必须报告全部提交/拒绝数，而不只报 9 个好结果 |
| XLPFE，*ACS Omega*（2022，DOI `10.1021/acsomega.2c01723`） | 用 PDBbind 2018 年前的 refined set 训练线性/非线性结构打分模型，并比较特征/模型 | 另在 2018 年后的 refined set 检验；也讨论相似度降低对性能的影响（CASF-2016 另列比较） | 最接近“冻结模型后打分新复合物”；时间划分和相似度审计比随机按行拆分可信 |
| Moore 等，自动化 RBFE from SMILES（2023，DOI `10.1038/s42004-023-00859-9`） | FEP 主要用固定力场和模拟协议，不靠每个 ligand 用实验 ΔG 拟合回归系数；相对能量需网络及绝对零点 | 对比人工建模的 reference pose 与不同自动 docking pose；明确指出人工 pose 基准不能直接推断自动前瞻应用效果 | **结构生成本身也要盲测**：用 test ligand 的共晶 pose 与从 SMILES 预测 pose 是两种不同任务 |
| Boltz-2（2025，DOI `10.1101/2025.06.14.659707`） | 神经网络参数由训练集拟合；测试时固定权重，不对 FMO87 每个靶点重训 `α,β,γ` | FEP+ 相关 held-out affinity 子集，论文说明对大多数测试集做蛋白序列层面过滤；四靶点 Table 12 的 FMO 行是另一物理 baseline | 区分“冻结大模型直接预测”与“用测试集本身做 per-target 校准”；两者数字不能等价比较 |
| Leak-Proof PDBBind 不确定性研究（2025，DOI `10.1038/s41598-025-27167-7`） | 在 LP-Train 训练模型、LP-Val 选择超参/停训，再在 LP-Test 及其他外部集评价 | 显式警惕蛋白与 ligand 的近重复产生 soft overlap，并检验测试域和不确定度 | 外部拆分仍要审计目标/骨架相似性；单一 MAE 或总体相关无法证明新系列可靠 |

**两条路线不能混用：** 物理 RBFE/ABFE 在冻结力场下计算能量（RBFE 以一个已知实验 anchor 定零点时，仅改变绝对偏移，不改变系列内排序）；经验能量模型以部分 ligand 的已知实验标签拟合 `α,β,γ`。后者必须清楚公布训练 ligand 和测试 ligand 的边界、适用域及校准所花标签成本。注意 Moore 等的 RBFE 论文在其**回顾性**对实验的绝对 ΔG 比较中，用同系列**全部实验值的平均数**给相对自由能定偏移：这种绝对误差不能当作“只用一个 anchor”的盲测 MAE；排名/ΔΔG 指标不受常数偏移影响。上述文章中也有回顾性或筛选偏差，并不是都达到了前瞻性盲测标准。

## 2. “固定参数＋新几何”至少是三个独立问题

1. **新化合物 / 同靶点、同系列**：部分实验化合物作 anchors，`α,β,γ` 只在 anchors 上拟合并冻结，测试从未参与拟合的新 ligand；这是检验 SophosQM 在先导优化中是否有用的首选任务。
2. **同一化合物 / 不同姿势或蛋白构象**：固定 ligand 和系数，改变结构生成（晶体、Boltz、对接、JNK1 flip/noflip）；衡量**几何鲁棒性**，但不是新 ligand 的预测泛化。
3. **新靶点 / 新系列**：不能移植为旧靶点专属的 `α,β,γ` 后宣称是同一模型；应使用多靶点预训练后整体冻结，或明确新增靶点的实验锚点成本。四个目标不足以强力证明跨靶点通用。

“几何是新的”也不自动等于盲测：若新几何是用**测试 ligand 的晶体复合物**做模板、根据测试集 ΔG 挑选 poses/水/质子化、或看到测试结果后再选较优协议，都泄漏了测试任务的信息。用晶体结果作**单独的 oracle-pose 对照**是可以的，前提是不要把它冒充只输入序列和 SMILES 的生产预测。

## 3. FMO87 的最小可执行、预注册设计

### 3.1 冻结数据与划分

- 输入 `data/fmo87_cohort.csv`：CDK2 `16`、JNK1 `21`、P38 `34`、TYK2 `16`，每个 case 都有公开的 `experimental_dg_kcal_mol`；不存在“只有四条有 G，另外 83 条没标签”的问题。
- 示范性的**回顾性保留集**：仅根据目标、SMILES/骨架多样性、唯一身份、固定 seed 选训练锚点 `8/8/10/8`（合计 34）；其余 `8/13/24/8`（合计 53）为测试。不得用实验 ΔG、FMO 得分或 Boltz 亲和力挑锚点。实验前冻结名单、标签遮蔽规则、hash、失败处理与统计口径。
- 因本项目已经用**全 87 标签**研究协议/指标，所以这一拆分无论现在怎么固定，都只能称为 **retrospective held-out-style audit**，不能称真正前瞻性盲测。真正盲测要找尚未查看标签的新增同系列实验，或未来合成测定的 ligand；确定 pose 处理时也不能利用这些标签。
- 同一靶点内部若跨多个 scaffold/实验批次，训练和测试必须记录 scaffold/series 和 assay 一致性；不存在可信时间戳时不要声称“时间外推”。先做同系列新 ligand，再单列严格的 scaffold-out 挑战。

### 3.2 冻结结构和特征

- 固定 QM 版本、fragment、水/质子化、5 Å 截取、FMO 能量汇总、clogP 软件/版本，以及缺失数据处理；**拟合和测试必须同一特征定义**。SophosQM 官方 MP2/VolSurf 系数不能直接套进本项目 DFTB3/RDKit 特征。
- 几何条件 A（上限对照）：用公开静态实验/建模 SDF 姿势，训练和测试均如此；说明测试姿势有结构先验，不能称 pose-blind。
- 几何条件 B（部署测试）：只用 anchor 的蛋白参考结构/指定序列和测试 ligand 的 SMILES，生成测试 poses（如同一 Boltz 版本或固定 docking 工作流）；固定种子、结构选择和 QA 阈值，不用 test 实验 ΔG、test ligand 共晶 pose 或 FMO 结果选“最优”。anchor 和 test 采用一致的 pose 生成方法训练模型；另外可单列“用 anchor 晶体特征拟合、在预测 pose 上打分”的几何迁移实验。
- 几何条件 C（敏感性）：同 ligand 多候选姿势/flip-noflip，在**完全冻结**参数和无标签选择规则下同时报告 score 的离散度、选错姿势率、失败率。不能为每个 ligand 看实验挑最接近 ΔG 的 pose。

### 3.3 两个明确分开的校准问题

**目标内监督校准：** 在每个 target 的 8–10 个 anchor 上，用训练数据确定 `α,β,γ`（或在 anchors 内做标准化/岭回归、选正则强度）；冻结系数、变换、姿势规则后，对本目标剩余所有 test ligand 预测。小样本两特征可能共线：报告条件数/置信区间；正则超参数只能在 anchors 内用嵌套验证选择。推荐先并列三种模型：`clogP-only`、`tPIE-only`、`tPIE+clogP`，另比固定符号 raw tPIE 排名。

**无标签泛化：** 用完全不含该靶点/该系列的外部训练数据学得共用系数，再冻结后在新靶点预测；若允许只用新靶点 4 个 anchor 定一个常数偏移 `γ`，要明确标注 **4-shot intercept calibration**，不要当无标签评估。偏移只改 kcal/mol 绝对误差，不改变同靶点 Pearson/Kendall；若用 4 点拟合 `α,β,γ`，剩余自由度太小，不建议报告为稳定的 SophosQM 复现。

### 3.4 只在测试集计算成绩

- 先写出 `train_ids.json`、`test_ids.json`、`pose_protocol.json`、模型系数和冻结时间/哈希；生成一份只含 test ID、结构 hash、预测 ΔG 的文件后，再关联测试实验标签。
- 每靶点测试 MAE、Pearson、Kendall、pairwise ΔΔG 误差和失败率；总结 target-average 时明确按**测试集**数目加权，并逐靶点报告样本数、置信区间，不要把训练点混进 Figure 19 式全局散点。
- 主要科学问题是新 ligand 的**相对排序是否超过基线**，次要问题才是 kcal/mol 的绝对误差；每靶点内算排名，避免四靶点间的平均 ΔG 差异抬高 pooled 相关。绝对 MAE 在四个靶点分别评价；不允许用 test 平均 residual 做 centered MAE 来冒充部署绝对误差。
- 对结构产生失败、FMO 不收敛、带电状态变化和变更姿势的 case 记录原始 ID 和原因，报告全体覆盖率及成功子集指标；若多协议/多个训练比例，必须在训练内选择主结果，其他标为探索性。

**成功标准应先写后看**：在目标内 held-out ligand 上，联合模型同时优于 `clogP-only` 与 raw tPIE，具有稳健的每靶点 Kendall/pairwise 优势、合理的绝对误差和高覆盖率；若仅拟合内提升而测试降级，应明确否定“校准改善泛化”。四靶点最小测试样本数仅 8–24 个，结论须附不确定性，最终要靠新同系列真实实验复核。

## 4. 与现有结果的关系

`docs/FINAL_REPORT_CN.md` 已经记录：本项目原始 FMO 的 target-average Pearson/Kendall `0.479/0.320`，全 87 个标签进行同靶点拟合后为 `0.609/0.500`，逐个留一为 `0.402/0.331`；TYK2 的留一 Pearson 为负。这是**已经观察到的校准不稳定性**，不是新协议的测试结果。应保持现有实验数据与模型不变作为探索性历史记录，把上面的冻结子集实验独立命名，禁止反复调参追论文数字。

## 主要原始来源

- SophosQM，Guareschi 等，*ACS Omega* 2023，DOI `10.1021/acsomega.2c08132`；本地补充表 `data/references/sophosqm_si.txt`。
- Automated LIE for aromatase，van Dijk 等，*Journal of Chemical Information and Modeling* 2017，DOI `10.1021/acs.jcim.7b00222`，完整正文 PMC5615371。
- FXR LIE / D3R Grand Challenge 2，Rifai 等，*Journal of Computer-Aided Molecular Design* 2017，DOI `10.1007/s10822-017-0055-0`，完整正文 PMC5767202。
- XLPFE，*ACS Omega* 2022，DOI `10.1021/acsomega.2c01723`，完整正文 PMC9245135。
- Automated RBFE from SMILES，Moore 等，*Communications Chemistry* 2023，DOI `10.1038/s42004-023-00859-9`，完整正文 PMC10140266。
- Boltz-2，2025 预印本 DOI `10.1101/2025.06.14.659707`；本地 `data/references/boltz2.txt`。
- Uncertainty quantification / LP-PDBBind，*Scientific Reports* 2025，DOI `10.1038/s41598-025-27167-7`，完整正文 PMC12678490。
