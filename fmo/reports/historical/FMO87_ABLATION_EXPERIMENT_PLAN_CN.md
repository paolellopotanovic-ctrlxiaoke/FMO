# FMO87 全因素消融工程计划与运行台账

日期：2026-09-30。目标不是从 87 个实验标签中挑出最好看的组合，而是区分：(i) 单一物理/计算因素的边际影响；(ii) 计算协议替换；(iii) 可解释的协同/抵消；(iv) 当前数据无法识别的因素。所有已尝试协议都保留，不因结果差而删掉。实验标签在首轮因素和协议冻结前不得用于选参。

## 0. 队列与预注册分析约定

- 固定队列：CDK2 16、JNK1 21、P38 34、TYK2 16，共 87 个 target-ligand 唯一键。每一个进入主表的分数都须 87/87；失败只能报缺失/失败，不准按成功子集与全队列比较。
- 主排序指标：每靶点分别算 Pearson r 与 Kendall tau，再按靶点样本数加权。实验 `dg` 越负表示结合越强；所有预测分数统一定义为越低越强，因此正相关为方向正确。另报告每靶点数值、靶点内成对排序准确率、分层 bootstrap CI 和相对 tPIE 的配对 CI。不得把跨靶点原始 QM 绝对能量池化排序。
- 预先指定基线：旧 FMO2-DFTB3/PCM tPIE。任何新增项先报告“单独项”和“基线 + 单项”，再报告预定义组合；未经外部验证或嵌套交叉验证的全 87 标签最优组合只叫探索结果。
- 物理项相加只在量纲/状态循环/几何相容时称为近似热力学组合。不同 Hamiltonian、不同 pocket、PCM 与 PBSA 的交叉组合单列为经验消融，不能宣称严格自由能分解。
- 几何刚体对齐只允许对蛋白共同原子求 Kabsch 变换，再把配体整体带入同一受体坐标框架；禁止为追求更好的亲和力分数而独立旋转/平移配体。刚体对齐不改变分子内构象，也不等于 pose docking。
- 每项需记录 protocol JSON、原始输入与日志哈希、程序/参数版本、成功数、评分公式、全量逐例值、靶点内指标和审计状态。

## 1. 因素及实验矩阵

| ID | 科学问题/对照 | 计算及输出 | 当前状态与下一步 |
|---|---|---|---|
| E00 | 原始 tPIE 是否可复现？ | AB FMO2/PIEDA，筛选蛋白-配体 pair 的总值；另保留 ES、EX、CT+mix、DI、SOLV 分项 | 已有 87 例；重审 pair 汇总与日志哈希、方向和队列身份 |
| E01 | PIEDA 哪个分项携带/破坏排序？ | tPIE 基线；各 PIEDA 子项独立；tPIE+单项；全部 2^k 无拟合子集组合；禁止把子项都称独立物理能 | 有 87 例 `component_scores.csv`；做完整目标分层与不确定性分析 |
| E02 | PCM 逐层响应的边际作用？ | 同一 DFTB3/PCM 三态日志从 Eunco+D → EunD_es → EunD+so；单独增量、累计项、tPIE+各增量仅作经验比较 | 有 87 例日志/诊断；复核定义并报告增量，而非把三个累计值当三项相加 |
| E03 | native subsystem 分解中的 interaction 项是否等同 tPIE？ | 对每个 AB 用同一 FMO2 Hamiltonian 输出 `DEint/DEsolv/DEmonomer/DE`；另算 AB PIEDA 跨 A/B pair 和；对比 `DEint`、tPIE、手动三态 | 仅少数 native validation；需先 1-case 试跑核实 GAMESS 组分定义，再决定 87-case 批量运行 |
| E04 | 单体/变形项在三态中贡献多大？ | 从 native output 独立得到 monomer destabilization 与 interaction，不用总值反推；与固定结合几何 B 的手动 A/B/AB 闭合检查 | 单例已核对；全量 native 若 E03 可稳定解析则补跑 |
| E05 | B bound geometry 与游离态几何的影响？ | 同一 AB/A，比较 B_bound 与 B_relaxed；再用多构象游离 B 的低能态/系综，不把单个局部极小值叫自由能 | 已有 DFTB3/PCM B-relax 87 例；审计优化收敛、几何、relax penalty、分层效应；HF 气相版需另行计算并标协议 |
| E06 | PCM 与 PBSA 的极性/非极性贡献？ | PCM 项、PB 极性、SA 非极性分别单项、累计项；tPIE+PB/SA仅经验模型；气相 FMO+PBSA 是独立热力学启发模型 | 已有气相 HF/PBSA 87 例及 PCM/PBSA hybrid diagnostic；修正总分字段审计，拆分每项且不重复溶剂响应 |
| E07 | PBSA 与 QM 的结构/电荷/原子映射不一致是否造成偏差？ | 逐态比较 A/AB/B 原子身份、口袋残基、cap、原子电荷/半径、三态 PQR 闭合；重建 same-residue-pocket PBSA 与备选共同原子策略 | 现有 PBSA 三态内部 A+ligand=AB 已审；QM cap 与未封端 PQR 有跨模型差异。先出逐例差异审计，再生成可辩护的 matched-pocket variant；cap 电荷不可臆造 |
| E08 | 蛋白/配体结构处理影响？ | protein frame RMSD/Kabsch、共有原子/链/水/altloc/质子化审计；比较原始 pose vs 标准化/最小化 pose、固定 vs 柔性蛋白、水保留 vs 删除、共同 pocket cutoff | 有多个几何与部分 FMO sensitivity 运行但不全是匹配 87 例；先完成结构审计，再对可比协议做 87 例 FMO/PBSA 重跑 |
| E09 | 需要何种几何对齐？ | 每个复合物蛋白对 target reference 做共同重原子 Kabsch；记录蛋白 RMSD 和变换后配体坐标；比较对齐前后原子映射/重原子接触。整体刚体不应改变内禀能量，坐标框架检查为主 | 实际 pose 已基于共同 prepared protein 生成；仍需 87 例数值审计。不能用独立 ligand fit 作为能量“优化” |
| E10 | relax 是否只是局部极小值？ | 每配体多起点 MMFF/OpenFF conformer search；记录 bound-like basin、最低能构象、能量分布、RMSD、重复性。和单次 B relax 区分 | 现有是单一 DFTB3/PCM 优化；多构象采样尚未做 |
| E11 | 熵项是否贡献独立信息？ | 先做可复现 ligand conformer ensemble 熵/构象选择代价 proxy；更强方案是 complex/free explicit-solvent MD 的重采样与 quasi-harmonic/dihedral entropy，至少独立重复；熵以 kcal/mol 的 `-TΔS` 进入组合 | rotor count 不是熵。ensemble proxy 可跑；MD 预算/力场有效性需先做 pilot，不将有限短轨迹包装成收敛自由能 |
| E12 | receptor cutoff/fragmentation/capping 效应？ | 同 Hamiltonian、geometry、PCM 的 3/5 Å pocket；残基并集固定、 cap 模式一致；比较 FMO 得分、收敛率、pair closure、相对排序 | HF 3 Å 已完成，DFTB PCM 5 Å 已有；cross-method 不能当 cutoff ablation。需同方法配对重跑至少 pilot，再评估全量成本 |
| E13 | 模型整合能否增益？ | 逐项 tPIE+PCM/PB/SA、native 分项已审计；8 组渐进特征做四折外层留一靶点、内层留一靶点选 Ridge 正则及特征组，严格只在训练目标拟合 | 已完成 87/87 匹配特征嵌套预测；`runs/fmo87/ablation/E13_nested_target_fit/`。内层选模型外测 `0.589/0.459`，固定组 `ES+SOLV` 外测 `0.608/0.465`，更多特征没有稳定提高；外部独立队列仍缺 |
| E14 | 晶水是否破坏固定姿态的 tPIE 排序？ | 在仅有晶水的 P38 34 例上固定全部其它输入，单独去水、重新算 FMO2/PIEDA；另外 53 例不改；核对逐例 pair/TIE 与分项变化 | 已完成全量，见下方台账。已写入的审计目录历史名称为 `E13_P38_no_waters`，**实验矩阵编号以 E14 为准**，不与上行 E13 模型整合混淆 |

## 2. 明确的输出分层

### 已有数据可完成，不需 QM 重算

E00/E01 PIEDA 分项；E02 PCM 三层能量与增量；E05 既有 B-relax；E06 已有 HF gas、PB polar、SA、PCM/PB diagnostic；已有几何最小化 RMSD；这些均重新进行 87-case identity audit、靶点分层指标、pairwise concordance、bootstrap uncertainty 和误差/波动分解。

### 必须重算的量子/结构实验

- E03/E04 全队列 native decomposition 与 AB PIEDA interaction-vs-tPIE closure：先单例试算、判据通过后再排队。不能拿 native 少数样本外推 87 例。
- E05 HF/STO-3G gas 版 B-relax 或明确保持 PCM relaxed 只作旧协议实验。优化后需 AB/A/B 采用定义一致的几何与 Hamiltonian。
- E07 matched PQR：只用有可追踪原子电荷/半径的同一残基/原子映射；QM link-cap 电荷若无可靠模型，须作为边界敏感性而非“修正成真值”。

**E25 QM/PBSA 受体封端边界：同力场 PQR 参数敏感性已完成 87/87，跨 Hamiltonian 等同与严格数值收敛仍未证明。** E04/E17/E27 已验证同靶点 QM A/AB 的真实受体和 ACE/NME 封端。新增 E25 从现用 PDB2PQR 自带 `AMBER.DAT`（SHA256 固定）读取原版 ACE/NME 电荷与 PB 半径，原 PQR 蛋白真实原子、电荷/半径与 B 逐字节不变，仅为 PBSA A/AB 添加 E04 冻结 QM cap 原子；四靶点 A 净电荷与 QM 一致，`PQR A+B=AB`、全部 87 例 PBSA 输出及三态差值独立复核完成：`runs/fmo87/ablation/E25_amber_capped_pbsa/independent_audit.json`。CDK2/17、P38/2ee 的 `0.5/0.375/0.25 Å` 网格先导保持 cap 效应方向但未达到严格收敛；没有事先为新模型选权重。封端后 PBSA polar 加 `ES+SOLV` 仍强烈不利，说明缺 cap 不是唯一根因。**后续如要宣称 PB 极性项的 0.1 kcal/mol 精度或优化 PB 数值设置，必须做预注册的全目标网格/半径敏感性；QM DFTB 与 PBSA/AMBER 电荷及溶剂响应仍不是同一个 Hamiltonian，E25 不等于实验结合自由能。**
- E08/E12 结构、柔性蛋白、explicit waters、3/5 Å 收敛诊断：冻结各统一协议，先 pilot 再全量；禁止混合不同成功子集。

### 熵的证据等级

1. RDKit/MMFF 多构象 Boltzmann/构象选择代价是**配体内构象熵代理**，不等于 binding entropy；需报 conformer coverage、温度、能量模型、RMSD basin 定义。
2. rotor count 是描述符，不是物理熵，作为 benchmark only。
3. 结合/游离体系的构象熵需采样两态；短 MD 或单构象 harmonic 估计只作敏感性检查。要称收敛的 ΔS，需独立 replicas、平衡/重采样诊断与不确定性；若计算预算/力场不允许，明确未完成，不能制造一个“熵值”。

**E11 真正结合熵的工程验收（与已跑的扭转角代理分开）**：输入至少需要同一结构/质子化/离子浓度与电荷方案下的游离受体 A、游离配体 B、复合物 AB 各自显式溶剂状态，预定义 1 mol/L 标准态；AB→A+B 的平动/转动自由度、蛋白/配体构象和溶剂相关熵不能用仅 B 的扭转角直方图代替。先写冻结的自由能/熵估计器及热力学路径（如温度导数或可重复的束缚/解束缚热力学循环）和重叠/闭合单元测试，再做四靶点各三独立重复的长程 pilot。验收至少包括采样分块稳定、跨重复 CI、状态重叠、能量/标准态与电荷守恒、在不挑成功配体子集下 87/87 全部纳入；任一门槛失败时**不能进入 87 例熵排序**。现有 2/50/500 ps 和正在运行的 2 ns 实验仅探索 B 扭转角边际分布的采样时间，不具备上述自由能闭合或物理意义上的 `−TΔS_binding`。

## 3. 统计与成功标准

- 所有分数方向、单位及 87 个键先通过自动校验；每个主实验必须 87/87 可评分，或报告失败矩阵且不报全队列指标。
- 主表含 `n`、每靶点 Pearson/Kendall、样本加权值、相对基线 Δr/Δtau、按靶点分层 paired bootstrap 95% CI。单个组合改善须差值 CI 支持且跨靶点方向不由一个 target 独占；所有组合需提示多重比较和探索性。
- 对组合开展 exhaustive subset screening 时，主报告只按预定义公式汇报；标签驱动的最佳子集必须嵌套交叉验证，并且不可宣称外部验证。
- “提升”须定义为同一 87 队列上的相关/成对顺序指标改善；不把绝对能量更负、单靶点变好或回归拟合内 R² 叫作泛化改善。

## 4. 当前运行台账

| 阶段 | 状态 | 产物/验收 |
|---|---|---|
| 旧数据的快速、无拟合 89-score 组合探索 | 已完成，属于初筛不是完整工程 | `runs/fmo87/ablation/`；本计划要求继续扩展并标注其边界 |
| 所有输入/输出源身份审计及新汇总器 | 进行中 | 已修正 HF gas+PBSA 总分列误用；新旧三态几何版本不同，固定统一几何后重算 |
| E03 native 87-case + PIEDA closure | 已完成 87/87；**仅验证 E08 逐配体 AB 口袋** | `runs/fmo87/ablation/E03_native_fixed_protein/`；DEPLd/DEsolv/DEint/DE 内部闭合、DEint 与 E08 tPIE 最大差 0.005 kcal/mol；不能据此验证 E04 并集口袋三态 |
| E07 matched PBSA rebuild | 已完成 87/87 并审计；原跨 QM 组合**口袋未对齐，现已以 E17 重做消融** | 纠正 `/usr/bin/apbs` 误用后 AmberTools `pbsa` 87/87 成功；PBSA 用同靶点并集口袋，而 E08 QM 用逐配体口袋；E17 并集选择性 PIEDA+PB/SA 87/87 已重新分析；`runs/fmo87/ablation/E17_shared_pocket_ablation/`。PBSA/QM 不同 Hamiltonian/原子帽，不能称严格热力学闭合 |
| E15 口袋协议审计/并集 native 闭合 | **87/87 完成并独立审计** | `runs/fmo87/ablation/E15_pocket_protocol_audit.json`：82/87 受体片段数不同，E03(native 旧口袋)与 E04(独立三态并集)绝对偏差中位 0.433、最大 1.579 kcal/mol。基于 E04 AB 另做并集 `MODMOL=24` native 87/87，独立审计 `runs/fmo87/ablation/E15_shared_native_audit.json`：对 E17 tPIE 最大差 0.004 kcal/mol，对 E04 A/B/AB 最大差 0.00050 kcal/mol；DEsolv 单项 `−0.486/−0.403`，tPIE+DEsolv `0.054/0.071`，DEPLd 单项 `−0.244/−0.158`。 |
| E15/E16/E17 全量 AB 与选择性 PIEDA 分项辨识 | E15 87/87、E16 P38 34/34、E17 选择性 87+34 **全部完成并审计** | E17 同靶点并集口袋 `MODMOL=1` 后 tPIE `0.308/0.246`、ES `0.555/0.421`、ES+SOLV `0.624/0.487`、ES+CT+SOLV `0.636/0.495`；19 水去除后完整 tPIE `0.521/0.368`。E15 全量 AB 的 EX/CT 在 56/87 例 >1000 kcal/mol 大反向对消；**最终 PIEDA 分项只采用 E17**；`runs/fmo87/ablation/E17_shared_pocket_ablation/` 提供 87 行、31 子集和 CIs；31 子集增益置换校正 p≈0.046，仅限同队列探索。 |
| E07 原子级跨模型对齐 | 已完成 87/87；cap 物理敏感性未重算 | 87 例 B 单体的 QM/PBSA 原子元素顺序、坐标完全一致；QM AB 前 A 后 B 的坐标也与独立两态严格闭合。4 个靶点 PBSA A 原子均与 QM A 对应，但是 QM A 因 ACE/NME/HOP 比 PQR A 分别多 96/84/120/96 个原子；`runs/fmo87/ablation/E07_cross_model_atom_audit.json`。因此“配体姿态没对齐”已排除，**溶剂/受体边界差异仍在**，不得视作已修正跨模型物理 |
| E08/E09 common-frame geometry audit | 已完成并已修复主 FMO | 87/87 都是自由蛋白 minimization；prepared-vs-minimized protein RMSD 中位数 1.81 Å，Kabsch 后 1.78 Å；固定蛋白 tPIE/三态/native 均完成 87/87 |
| E14 结晶水是否造成排序噪声 | P38 34/34 已完成并审计；另 3 靶点均无晶水 | P38 每例去掉原始口袋中的 10–15 个 water fragments，原有 5 Å 固定姿态、FMO2-DFTB3/PCM/PIEDA 与全体输入保持不变。所有 34 例日志、pair/TIE、配体坐标对齐通过；未去水的另外 53 例保持原分数。87 例完整 tPIE 从 `0.356/0.255` 变 `0.556/0.395`，Δr 配对 CI `[+0.033,+0.364]`。仅将 P38 的 EX 换成去水计算值可得约 `0.575/0.418`，主要效应在 EX；历史目录 `runs/fmo87/ablation/E13_P38_no_waters/audit.json`；不可宣称去水在外部靶点有益 |
| E05 固定几何 B 单体 QM relax | 已完成 87/87 并审计 | 固定几何 A/AB 下重新独立优化 87 个 B；`runs/fmo87/ablation/E05_fixed_protein_ligand_relaxed/`。三态分数约 `−0.252/−0.162`，比同几何未 relax 的 `−0.018/0.043` 更差；单起点不能替代多起点 |
| E10 双起点 B relax | 已完成 87/87 并独立审计；非复合物多起点 | 每例另选 32 个 ETKDG/MMFF 构象中最低能者，以原子顺序不变的第二起点重跑 B 的 FMO1-DFTB3/PCM 优化。1 例 160 步未收敛，改 500 步定点续跑收敛；87 个起点均按种子/原子/能量重建，并核对 deck、QM 参数和 GAMESS 日志。第二起点有 48/87 例更低，选两者中更低能后，三态 Pearson/Kendall `−0.217/−0.132`；`runs/fmo87/ablation/E10_second_start_qm_relax/audit.json`。两起点仍不代表完成构象自由能采样 |
| E11 conformer entropy proxy / 两态采样先导 | 两次 conformer 87/87；四靶点双态各三重复 24/24；CDK2 50 ps 六条重复通过身份审计、**扭转角分布未收敛** | 原 2 ps×24 条各20帧未达到门槛；CDK2 50 ps×6 条各500帧已过帧数/有效样本门槛，但 free 前后分布 JS=0.265>0.15；`runs/fmo87/ablation/E11_explicit_two_state_50ps_CDK2/entropy_audit.json`。均只含配体扭转角边际熵，仍非完整 87-case binding entropy。 |
| E11 采样长度敏感性扩展 | **CDK2 500 ps、2 ns 两态各三重复均 6/6 成功，质量门槛都未通过** | 500 ps×6 各500帧：bound JS=0.241、free ESS105.1、JS=0.187；2 ns×6 各1000帧：bound JS0.213、free ESS64.7、JS0.412，仍不满足预注册 JS≤0.15/ESS≥150。2 ns 轨迹 SHA/身份和蛋白接触审计通过，但 `converged_torsional_proxy=false`、`full_binding_entropy_available=false`、`ranking_score_generated=false`；`runs/fmo87/ablation/E11_explicit_two_state_2ns_CDK2/entropy_audit.json`。这是 CDK2 的配体扭转角诊断，不是 87 例结合熵。 |
| E18 固定受体弱配体约束的复合物几何 ablation | **87/87 选择性 QM 同口袋重算、闭合及成对分析完成** | 87 例 k25 与原 k100 都是固定受体配体局部最小化；受体坐标、cap/片段、QM 方法相同，只替换 B 坐标。重原子 RMSD 中位 `0.070 Å`；完整 tPIE `0.308/0.246 → 0.462/0.371`、ΔPearson CI `[+0.072,+0.247]`，ES+SOLV `0.624/0.487 → 0.618/0.487`、CI 跨零。最初 1 例因打印精度误报闭合失败，保留旧日志、修复容差并重跑。逐例数据和靶点 CI：`runs/fmo87/ablation/E18_geometry_ablation_analysis/`；不得将 k25 新 AB 与旧 PBSA/三态当同姿态相加。 |
| E19 同输入气相 FMO2-DFTB3/PIEDA | **87 例首版尝试仅 12 例收敛，且意外同时删除经验色散 DC；不得生成 87 例指标** | 首版从 E17 AB 删去 PCM 三块和 `$DFT DC` 一块，75 例片段 SCF 振荡，且**并非纯 PCM 单变量**：CDK2/17 首版 tPIE `−20.335`，保留 DC 仅去 PCM 是 `−72.006`，水 PCM `−69.034` kcal/mol。正确纯 PCM 消融先导 CDK2/30 仍不收敛；`ITYPMX=2` 电荷 DIIS、`ITYPMX=-1` 轨道 SOSCF 未修复原失败例。完整日志/审计 `runs/fmo87/ablation/E19_shared_pocket_gas_pieda/`、`E19_gas_keep_dc_pilot/`、`E19_gas_diis_pilot2/`、`E19_gas_soscf_pilot/`；不能用幸存者外推全 87。 |
| E20 只改变 PCM 静态介电常数 | **87/87 成功且逐例 deck/日志/Pie、分项及 CI 审计完成** | E17 同一个 AB/GAMESS FMO2-DFTB3/PIEDA、PCM 水腔及非静电参数，唯一增加 `$PCM EPS=4.000`，日志 active dielectric 逐例确认；完整 tPIE `0.308/0.246 → 0.295/0.250`，ES+SOLV `0.624/0.487 → 0.619/0.499`，后者 ΔPearson CI `[−0.019,+0.014]`；SOLV 单项 `0.220 → 0.241` 但组合无独立稳定增益。逐例和按靶点 CI：`runs/fmo87/ablation/E20_pcm_dielectric_analysis/`。**这是连续介质响应敏感性，不是真气相**。 |
| E21 同几何 PBSA 极性/非极性消融 | **87/87 PBSA 完成，四靶点受体 PQR 哈希相同** | 用 E18 k25 ligand 姿态重算 PBSA，但用 E08 k100 姿态冻结口袋选择；配体 PQR 原子身份/电荷/半径及坐标逐一核对，`runs/fmo87/ablation/E21_pbsa_k25_fixed_pocket/full/audit.json` 成功。k25 `ES+SOLV` `0.618/0.487`，加相同姿态 PB polar `0.182/0.152`，配对 ΔPearson CI `[−0.659,−0.255]`；SA 增益 CI 跨零。逐例和双条件 CI 在 `runs/fmo87/ablation/E21_geometry_pbsa_analysis/`；跨 QM cap/电荷仍不一致。 |
| E22 P38 水×几何 2×2 | **去水 k25 新 QM 34/34 成功；四格共 87 行的成对分析完成** | P38 固定 5 Å 并集口袋，仅在 QM 评分阶段移除统一的 19 个晶水片段；配体姿态分别来自保水最小化 k100/k25，其余 53 例各自保持不变。tPIE 四格 `0.308/0.246`、`0.462/0.371`、`0.521/0.368`、**`0.584/0.424`**；水效应×几何 87 例加权 Pearson 交互 `−0.091`、CI `[−0.161,−0.032]`。对 ES+CT+SOLV 去水仍不提升，最好保水 k100 `0.636/0.495`；`runs/fmo87/ablation/E22_water_geometry_factorial_analysis/`。队列内组合选择，不是盲测。 |
| E23 更弱姿态约束 k=10 | **固定蛋白 87/87 最小化、新 QM 87/87、三档成对 CI 全部完成** | `runs/fmo87/ablation/E23_fixed_protein_k10/summary.json` 最小化与 `runs/fmo87/ablation/E23_fixed_geometry_k10_pieda/audit.json` QM 是独立步骤；新 AB 只改变 ligand 姿态，重原子 RMSD 中位0.139 Å，无重原子移动>1 Å。tPIE k100/k25/k10 Pearson/Kendall `0.308/0.246`、`0.462/0.371`、`0.533/0.412`；k25→10 ΔPearson CI `[+0.010,+0.138]`。ES+CT+SOLV `0.636/0.495`→k10 `0.644/0.504`，CI跨零。`runs/fmo87/ablation/E23_restraint_dose_analysis/` 提供87行及三档配对区间；未在k10重新算PBSA/三态，不得拼接旧几何项。 |
| E24 原有配体熵代理与并集 PIEDA | **两批已有 87 例代理完成 E17 同输入配对消融** | 已有游离配体 ETKDG/MMFF94s 64 构象(seed1)/128 构象(seed2) 的 `T×S_free`，按固定物理单位 `kcal/mol` 分别加 tPIE/ES/ES+SOLV/ES+CT+SOLV；ES+SOLV `0.624/0.487` 原始，对应 `0.623/0.498` 和 `0.616/0.495`，ΔPearson CI 均含0；逐例 `runs/fmo87/ablation/E24_shared_entropy_proxy_analysis/`。这不是完整结合熵，两个采样量与种子同时变化；独立的显式双态质量门槛尚待完成。 |
| E25 同源 AMBER 封端 PBSA A/B/AB | **先导 2/2、四靶点 87/87 全计算和独立审计、配对消融完成；网格仅先导** | PDB2PQR AMBER.DAT ACE/NME 电荷/半径参数 + QM 冻结 A/AB 中 cap 坐标；真实受体与 B PQR 原封不动，独立核对逐原子 87/87 AB 和 4 个 A、三态能量及 SHA。PB polar 单项 Pearson/Kendall `−0.479/−0.352 → −0.465/−0.344`，`ES+SOLV+PB` `0.240/0.173 → 0.238/0.155`，PB 对无 PB 的 `ES+SOLV` 增量 Pearson `−0.386`、CI `[−0.561,−0.242]`；两个病例 0.5/0.375/0.25 Å 网格 cap 差值方向一致但未严格收敛。结果与逐例：`runs/fmo87/ablation/E25_capped_pbsa_analysis/`，独立审计 `runs/fmo87/ablation/E25_amber_capped_pbsa/independent_audit.json`。QM DFTB 与 PBSA 电荷/响应仍不同，不可称亲和力真值。 |
| E26 P38 k=10 去水 3×2 矩阵 | **P38 34/34 新 QM 成功，完整六格 87 例、分项与交互 CI 全完成** | E23 k10 pose 与 E16/E17 同靶点去水受体/封端仅替换 ligand 坐标；**只在 QM 评分阶段删水**。k100/k25/k10 去水 tPIE `0.521/0.368`、`0.584/0.424`、`0.606/0.442`；对应保水 `0.308/0.246`、`0.462/0.371`、`0.533/0.412`。去水边际 Pearson增益 `+0.213/+0.122/+0.072`，k25→k10 加权交互 `−0.050`、CI `[−0.109,−0.003]`；ES+CT+SOLV 在三档仍保水优。`runs/fmo87/ablation/E26_six_cell_water_restraint_analysis/` 87 行完整配对；同队列探索非盲测。 |
| E27 共享口袋 QM 封端所属片段敏感性 | **87/87 冻结受体/封端映射、选定 PIEDA pair 和哈希审计完成；不等同 E25 真正重跑 PBSA** | E04/E17 冻结并集受体逐原子 A/AB 同一，封端片段 15/13/18/14 个；排除整个含 cap 片段 pair 的 Pearson/Kendall 从 `0.308/0.246` 到 `0.324/0.264`，配对 CI 跨零；被排除片段还含真实蛋白，**不能**称为 cap-only 更正。老姿态文件与冻结输入的差异逐例记录，最终以冻结 AB deck 为准；`runs/fmo87/ablation/E27_shared_cap_fragment_diagnostic/`。 |
| E28 两档网格全队列原 PQR 与封端 PQR 消融 | **新算 PBSA 269 状态，87/87 配对分项/排序与独立能量审计完成** | 同一 PQR A/B/AB、PB 参数/电荷/半径只将 `space=0.5→0.25 Å`；细网格封端 PB 单项 `−0.463/−0.340`，`ES+SOLV+PB` `0.239/0.155`，相比 `ES+SOLV` Pearson −0.386、CI `[−0.553,−0.244]`。排序负面结论跨两档网格成立，但 cap PB 绝对值按靶点中位可能大幅变化，TYK2 改变量 `+0.044→−0.198`（符号反转），**不等于 PB 绝对能量严格收敛**；`runs/fmo87/ablation/E28_full_grid_capped_pbsa_analysis/`，独立审计 `runs/fmo87/ablation/E25_capped_pbsa_grid_pilot/full_grid_0.25_independent_audit.json`。 |
| E29/E30 配体姿态 k25/k10 × PBSA 封端交互 | **两档几何各 87/87 无 cap 与同源 AMBER-cap PBSA；各几何 A/B/AB、QM 配体输入逐例对齐审计，三档 87 行配对交互完成** | k25/k10 不再拼接 k100 的 PB。`ES+SOLV` 三档 `0.624/0.487`、`0.618/0.487`、`0.634/0.508`；各自加同姿态 cap-PB 降到 `0.238/0.155`、`0.186/0.163`、`0.245/0.181`，三档 PB 增量 Pearson CI 都严格为负；封端×几何交互 CI 跨零。SA 带来的独立增益很小。完整三档每例 tPIE/ES/CT/SOLV/PB/SA 及交互见 `runs/fmo87/ablation/E30_geometry_cap_factorial_analysis/summary.json`；新 PBSA 输入 `runs/fmo87/ablation/E29_amber_capped_pbsa_k25/`、`runs/fmo87/ablation/E30_pbsa_k10_fixed_pocket/`、`runs/fmo87/ablation/E30_amber_capped_pbsa_k10/`。 |
| E31 MM/GB 正常模三态可行性 | **CDK2/17 三态拓扑与 B 单项试点完成；AB/A 未放行；87 例结合熵未运行** | OpenMM ff14SB/OpenFF→ParmEd AB/A/B 电荷与原子数闭合；必须用 mbondi2 补齐转换后全零的 GB 半径。修复后 B 的正常模能量最小化收敛、零虚频；AB/A 全蛋白分别 4847/4802 原子，一个双精度稠密 Hessian 至少 1.575/1.546 GiB，超过**本试点设定**的 0.2 GiB 单矩阵预算（真实工作内存更高）。不得将 B 的绝对谐振熵当结合熵，保留 MD 双态与标准态验收门槛。`scripts/pilot_fmo87_normal_modes.py` 与 `runs/fmo87/ablation/E31_normal_mode_CDK2_17_mbondi2/`。 |
| E12 controlled same-method cutoff/pocket | 已完成 87/87，PIEDA 子集已配对评估 | 固定蛋白、DFTB3/PCM 的 3 Å 同例对照，`runs/fmo87/ablation/E12_cutoff_subset_robustness.json`；3 Å 的 ES+CT+SOLV `0.626/0.467`，5 Å `0.676/0.542`，二者均比各自完整 tPIE 更强；不能声称 3 Å ES 比 5 Å ES 更强 |
| E13 训练内嵌套目标外测 | 已完成 87/87、输入哈希及训练系数留档 | 4 外层靶点，3 靶点内折选 8 种特征组与 Ridge alpha；外层标签只用于最终评价。内选 `0.589/0.459`，组固定 `ES+SOLV` `0.608/0.465`，含五个 PIEDA 分项 + PB + SA 的七特征组 `0.417/0.292`，未见“越复杂越好”；`runs/fmo87/ablation/E13_nested_target_fit/summary.json` |

| E32 游离配体构象数 × 采样种子 | **64/128 × 两个固定数值随机种子，四格均 87/87 完成；已分离 E24 种子/构象数混杂** | 固定 E17 QM、姿态、MMFF、温度与 RMSD；最小化全部收敛。`ES+SOLV` 原 `0.624/0.487`，四格 `0.617/0.501`、`0.616/0.495`、`0.619/0.494`、`0.623/0.505`；加 `+T×S_free` 的 Pearson 增量 CI 都跨零。同种子增加候选数对排序的结果随种子改变；bound-like 零概率不支持 87 例 `−RT ln p`。种子为 RNG 整数，不是运行日期；详见 `runs/fmo87/ablation/E32_free_ligand_factorial/`。 |

| E33 冻结复合物姿态→游离 B 单体应变 | **两套 128 构象数值随机种子各 87/87，E04 A+B=AB=E17 受体/配体几何逐例匹配，与 E05/E10 两起点 QM B 同例对照完成** | MMFF94s 同一配体重原子冻结、只优化 H 的 bound 能量减去独立 B 全自由最小值；按同 Hamiltonian 把局部位移与更低自由构象盆地拆开。与 QM 双起点 B relax 逐例 Pearson `0.957`、Spearman `0.860`，但直接将应变加 `ES+SOLV` 使 `0.624/0.487→0.058/0.121`，配对 CI 严格为负。只代表当前几何/MMFF 的蛋白结合前配体内应变，非结合熵或标准态自由能；`runs/fmo87/ablation/E33_bound_strain_factor_analysis/`。种子数字 `20261001/20261002` 是 RNG 整数，**并非实验日期**。 |
| E34 应变的七项 MMFF 几何/非键能量分解 | **87/87 同一冻结几何，七项每例和总局部应变闭合；重算 H 极小化数值漂移单列审计** | `ES+SOLV` 原 `0.624/0.487`，加键伸缩 `0.376/0.403`、加角弯曲 `0.437/0.393`、加两项 `0.262/0.301`；分别单项配对 ΔPearson CI 键 `[−0.431,−0.039]`、角 `[−0.298,−0.034]`。加扭转项 `0.648/0.496` 的增量 CI 跨零，靶点异质性显著，**不列作已证实提高**。高应变与重原子键长偏移相关但删异常样本仍伤排名；`runs/fmo87/ablation/E34_bound_strain_term_analysis/`。 |
| E35 弱化配体约束 × 游离应变 | **k100/k25/k10 三档 QM/PIEDA 87/87 既有冻结 AB 输入，新增 k25/k10 MMFF 2×87 全部几何/分项审计；交互及区间完成** | 三档 `ES+SOLV 0.624/0.487、0.618/0.487、0.634/0.508`；匹配各自几何再加 B 局部应变 `0.132/0.198、0.238/0.270、0.313/0.270`。k100→k10 应变增量的 Pearson 交互 `+0.172`、CI `[+0.057,+0.299]`，但 k10 单独加应变仍显著负向（Δ`−0.321`，CI `[−0.473,−0.141]`）。**几何改善噪声未能挽救加应变排序**；`runs/fmo87/ablation/E35_geometry_strain_factor_analysis/`。 |
| E36 独立来源/靶点验证 | **预先冻结，待执行：Galectin3 26 个中性配体；从未进入四靶点 87 例的规则选择** | 只按结构、配体中性/有氢和口袋大小选 `vendor/public_binding_free_energy_benchmark/fep_benchmark_inputs/structure_inputs/misc/galectin3_extra_{protein.pdb,ligands.sdf}`；未读取实验标签数值。PDB2PQR pH 7.4 AMBER，移除游离氯离子、保留原晶水；固定蛋白/水、配体 OpenMM/OpenFF k=10 kJ/mol/Å²，26 例并集 5 Å 同一受体/封端，FMO2-DFTB3/PCM 选择性 PIEDA。**固定无拟合评分四条**：完整 tPIE，ES，ES+SOLV（主要比较），ES+CT+SOLV（次要比较）；全部采用 GAMESS pair 表中的原始 kcal/mol 直接相加，不调权、不看标签排除样本。26/26 成功后才解盲 `r_exp_dg`，每条报告 Pearson/Kendall 与 ES+SOLV 相对 tPIE 的配对重采样区间；失败病例单列、不得挑成功案例声称独立验证。此验证不是原四靶点总体 87 例指标，且 PCM/几何协议的外推可能失败。 |
| E37 外部靶点几何强约束敏感性 | **E36 解盲后提出；必标注事后分析，不可当独立盲测** | 对同一 26 个 Galectin3 结构复用 E36 已冻结的 PDB2PQR 蛋白和标签脱敏 SDF，唯一改变固定受体最小化的配体约束 `k=10→100 kJ/mol/Å²`；为全部配体重建单独 5 Å 并集口袋/封端并重新计算 FMO2-DFTB3/PCM PIEDA，重新算 tPIE、ES、ES+SOLV、ES+CT+SOLV 与 E36 逐例配对差和几何位移。不得把 k=100 挑选为新的未看标签独立验证协议；E36 的反向表现保持原样报告。 |
| E38 外部靶点同受体纯几何臂 | **E37 运行中发现口袋已变，补充 E38 匹配口袋；解盲后分析** | 复用 E37 k100 同 26 个最小化姿态，但选取口袋仅用 E36 k10 的完整 26 例并集，冻结 E36 蛋白/水/封端/HOP；QM 只改变 AB ligand 坐标，其余 FMO2-DFTB3/PCM PIEDA 输入设置不变。必须逐例验证 receptor + cap 与 E36 完全相同；与 E36 对比才可归因纯配体几何变化，E37 自身是几何+口袋混合差，不能代替 E38。 |

**E37 运行后质量控制修正（2026-09-30）**：两个几何各自独立选择 5 Å 并集，使口袋原子数 418→446、ARG183/GLY182 等边缘片段/水归属发生变化。因此原 E37 行中的“唯一改变约束”仅适用于**最小化阶段**，不能用于归因最终 QM 分数差异；不改写 E36 预注册或挑选成功样本。E38 专门冻结 E36 口袋以正交隔离。

**历史计划状态勘误（2026-10-01）**：上表 E36 一行写的「待执行」是当时冻结分析方案的历史状态，**不是当前执行状态**。E36 的 26/26 QM 成功且已审计解盲，E37/E38 各 26/26 且属于 E36 解盲后的因素消融；E36 的 ES+SOLV `−0.538/−0.362`、E38 `−0.505/−0.424`。E39 FAAH 由于 SCF 未全部收敛而**未解盲**；E40 FXa 27/27 审计解盲，ES+SOLV `+0.529/+0.377`。冻结规则和原始失败全部保留，不能把 E37/E38 再称为新外部盲测。

| E39 FAAH 第二外部靶点盲测 | **E36 解盲后、FAAH 标签读取前冻结，24 个中性有氢分子** | 在同一个公开结构源 `misc/hfaah_{protein.pdb,ligands.sdf}`，按结构全部 24 例纳入（不按标签筛选）。与 E36 完全同样的 PDB2PQR pH7.4/AMBER 保水、固定受体/水、OpenMM/OpenFF k10、全 24 例共用 5 Å 并集、FMO2-DFTB3/PCM PIEDA 和无权重 tPIE、ES、ES+SOLV、ES+CT+SOLV 四分数；仍以 ES+SOLV vs tPIE 为主要配对比较。读取 `r_exp_dg` 须在 24/24 QM 分项、逐例 SHA 和 GAMESS 成功审计之后，不允许按效果调整水、约束或靶点后声称盲测。FAAH 与 Galectin3 同在 misc 目录但来自不同研究/蛋白/化学系列，仅可称靶点/研究分离，不可称来源目录分离。 |

**E39 盲态下技术性修订**：原始 `$SCF DIIS=.FALSE. SOSCF=.TRUE.` 前 7 例成功，第 8 例配体 `12` 的一个 dimer 最后 200 次迭代仍在约 1 kcal/mol 幅度振荡，因此初版 **7/24 不予解盲/排序**。在完全不读取 FAAH 的 `r_exp_dg` 前，试 `DIIS=.TRUE. SOSCF=.FALSE.` 数值求解器先导；若先导收敛，**新目录下 24/24 全部相同 solver 重跑**后才解盲，不混用两 solver 幸存样本。单改最大迭代数不足以证明收敛，因为末尾仍持续振荡。必须报告原始失败和新 solver 改动，严格检验两 solver 已成功病例的分数一致性；若无统一收敛则保留失败，绝不偷读标签或声称 24 例独立验证。

**E39 第二次盲态定位**：`DIIS=.TRUE.` 先导与原失败 dimer 迭代序列几乎一致，也失败。核对安装目录的 GAMESS `docs-input.txt` 中 `$DFTB ITYPMX`：默认 `0=Broyden` 做 SCC 原子电荷混合，`2=电荷 DIIS`，`−1` 才采用普通 GAMESS SOSCF/DIIS；所以仅改变 `$SCF` 根本不针对当前默认混合器。继续**不读标签**用固定同一结构 `ITYPMX=2` 对失败病例做先导；若有效，再在新输出根目录对全部 24 例统一重算、审计并报告与 7 例旧收敛输入的能量差。数据尚未解盲。

**E39 第三次盲态定位**：`ITYPMX=2/−1` 两个先导仍有不收敛 dimer；原先失败对涉及阴性 ASP403，距最近配体重原子约 `4.93 Å`，未见 <2 Å 严重重原子碰撞。新开独立 `E39_faah_cutoff3` 臂，预先指定 **3 Å 完整 24 例并集**、原始 k10 姿态、DFTB3/PCM 同设置、统一 24 例；先对失败病例技术先导、若收敛才完整 24/24 重算后解盲。该臂是**技术修订后的另一个 3 Å 物理口袋协议**，绝不当作原 5 Å FAAH 或 E36 5 Å 的同协议重复，不用 7/24 原 5 Å 幸存病例计算排序。3 Å 是此前 E12 在原 87 例已完整测试的口袋半径，不依据 FAAH 亲和力标签挑选。

**E39 最终先导门槛：失败、绝不解盲**。3 Å 配体 `12` 依旧有 SCF 未收敛 dimer；5 Å 初版 7/24 已完成、其余未算，三个数值迭代器先导及 3 Å 先导均失败。上述试验全程未读取 FAAH 的任何亲和力数值，也没有 24/24 完整一致协议，因此**没有 FAAH Pearson/Kendall，更不能基于 7 个幸存者报独立验证**。保存全部计算日志以便将来专门处理蛋白/配体电荷状态或替换 Hamiltonian，但那是新协议及新实验，不影响 E36 外部负结果。

| E40 FXa 第二次独立外部验证 | **E39 技术失败后、查看 FXa 实验标签前冻结：27 个中性显式氢配体** | 选另一结构来源 `opls_stress/fxa_yoshikawa_set_{protein.pdb,ligands.sdf}`，27/27 SDF 均有独立坐标与中性显式氢；原 PDB 两个 Ca 离子距最近配体重原子 ≥8.4 Å，**不在 5 Å 口袋内**，和原协议一致由 PDB2PQR 过滤（仍须承认远处离子静电被省略）。固定 PDB2PQR pH7.4/AMBER/保水、固定蛋白/水、OpenMM/OpenFF k10、全 27 例共享 5 Å 口袋、FMO2-DFTB3/PCM，GAMESS PIEDA 无拟合 tPIE/ES/ES+SOLV/ES+CT+SOLV；主要比较 ES+SOLV vs tPIE。必须 27/27 SCF 成功、日志/分项/SHA 对齐，**先归档盲预测再读 `r_exp_dg`**；如果失败则如实汇报无统计。由于已经看过 Galectin3，但没读 FXa 标签，FXa 可以检验预定四规则在另一个新来源/靶点的迁移，不能假装 E36 结果不存在。 |

**E40 输入阶段蛋白结构修订（标签未读）**：初版 27/27 k10 姿态已有，但 QM deck 的片段电子数检查在配体 `5f`、`A:CYS:191` 报奇电子并拒绝运行，非 GAMESS 计算失败。追踪到 PDB2PQR 生成的口袋含 CYS191–CYS220 约 1.996 Å 二硫键，分离为两个 CYS 片段时每片段是单电子奇数；已单独实现**只对完整 SG–SG 两端在口袋内的共价二硫键合并为同一闭壳层片段**，不增加虚构氢、不断开天然键，保持全队列同一受体/口袋并可对每例受体片段做哈希审核。新目录 `E40_fxa_disulfide` 复用完全相同蛋白/k10 配体姿态，但使用该结构修订后的 FMO 碎片；修订在读取任何 FXa 亲和力数值前完成，不能与旧未运行的错误输入混用。如新 27/27 无法通过奇偶性或收敛门槛，则 FXa 不解盲。

**E40 第二次蛋白修订（继续完全盲态）**：单独合并二硫键后首例 GAMESS 在距离门槛检测到跨片段重原子仅 `0.073 Å`，原子对应人为 NME 终端封端和真实 `CYS220/N` 重合。原因是 PDB 编号 GLY218→CYS220 缺 219，但实测 `GLY218/C—CYS220/N = 1.325 Å`，是完整肽键；旧代码只看残基编号 +1，误当链断裂补了 ACE/NME。新增可选 **按 C–N 1.15–1.65 Å 几何识别编号缺口处连续肽键**，与二硫键合并一起应用于新目录 `E40_fxa_repaired`。准备检查显示片段闭壳层、原始口袋没有重原子跨片段 <1.1 Å 非法碰撞；均对 27 例同一受体生效，且不得修改 E36 旧冻协议或读 FXa 标签直到全部成功。

**E40 执行结果（2026-10-01）**：按盲态冻结的四个无拟合 PIEDA 得分在修复输入后统一计算，27/27 SCF 成功，逐例审计输入/日志 SHA、冻结协议和蛋白/来源 SHA、配体对分项闭合；盲预测写入 `E40_fxa_repaired/blind_predictions.json` 后才运行解盲脚本。独立靶点 FXa 的 tPIE `0.473/0.382`、ES `0.528/0.394`、ES+SOLV `0.529/0.377`、ES+CT+SOLV `0.509/0.382`（Pearson/Kendall）。主要比较 ES+SOLV−tPIE 的 Pearson `+0.056`、Kendall `−0.006`，两项成对 bootstrap 区间均跨零；不能把 E40 当作四靶点选定组合的成功外推证明，也不能因 Galectin3 失败而抹去其正向结果。详情与原始审计见结果报告 E40。

## E76 全构象 ensemble（无 cherry-picking，预注册待执行）

每个 target/ligand 预先生成固定数量的独立 3D 构象（建议统一 4 或 8 个；相同生成器、
种子、受体和 QM fingerprint），全部完成 FMO2/PIEDA 四分数。使用
`scripts/aggregate_fmo_conformer_ensemble.py` 的固定 weighted-mean 作为主分数，同时
报告 median、全构象 Boltzmann 诊断和逐例 min/max；构象缺失、QM 失败、混合 fingerprint
或 ligand 数量不一致则整组 `FMO_REVIEW`，不能用幸存构象补齐。首轮比较单构象
`ES+SOLV`、全构象 `ES+SOLV`、全构象 `tPIE/ES/ES+CT+SOLV`；只报告预先冻结的全队列、
全 target 指标和跨构象稳定性，不能按 ligand 选择构象或聚合公式。先在原 87 例完成
内部盲态方法审计，再用 Galectin3/FXa/TNKS2 做独立迁移验证；当前尚无结果，不能
宣称会提高指标。

## 5. 运行顺序

1. 固定 cohort/hash，审计现有 feature sources 与字段语义；旧错误 score 列不得作为 gas ΔE 使用。
2. 完成并冻结本文件中原始数据 reanalysis、置信区间和逐项贡献；将每个结果标记为复算/新算。
3. 做结构映射与 native/PBSA 单例试验；失败先停下来修输入，不要启动 87 例。
4. 按 E03/E05/E07/E08/E10/E11/E12 分开输出目录和 manifests；每队列完整后才能进合并评分。
5. 生成总矩阵、失败审计、源数据哈希、完整逐例表和学长汇报版总结。

## 6. 未完成项目的明确边界及验收条件

**E41/E42 更新（2026-10-01）**：E41 复用同一 AB、同一姿态、同一受体的 E23 三档完整 87 例 PIEDA，在固定单位系数下逐项审计 `ES`、`SOLV`、`ES+EX`、`ES+CT`、`ES+DI`、`ES+SOLV`、`ES+CT+SOLV`、无 PCM pair 的 `ES+EX+CT+DI` 与原 tPIE；每档 7 条预设边际对比均按靶点分层配对重采样，独立验证已有 87 例汇总分数，逐例最大五项小数舍入闭合残差 `0.014 kcal/mol`。新脚本 `scripts/analyze_fmo87_pieda_marginals.py`，结果 `runs/fmo87/ablation/E41_pieda_marginal_analysis/`。这属**原队列事后归因**，不是 21 次比较里挑最大值替换无拟合主评分；详见结果报告 E41。

**E42 固定可迁移政策**：`scripts/freeze_fmo_affinity_protocol.py --cohort data/fmo87_cohort.csv` 只保留结构输入的白名单列，不使用实验 `ΔG` 或原 SDF 的亲和力属性；在打分前按 `target+ligand` 去重、检查 SDF 质子化/电荷/显式 H/3D 构型、蛋白链及蛋白–配体重原子最近距，所有靶点按全队列共用 5 Å 口袋和 `k=10` 冻结唯一 ES+SOLV 主分数。87/87 结构 gate 通过，冻结数据 `runs/fmo87/ablation/E42_frozen_production_policy.json`；**该文件仅是上游分析/准入政策，不会替代实际 GAMESS 的输入/SCF/PIEDA 审计，也不保证新靶点算例收敛**。任一靶点遇结构/碎片/SCF 问题时停并在盲态统一修订整组，不能仅对成功分子挑阈值或方法。若将来根据配体化学类别动态切换计算方法，分类规则和阈值必须在新数据标签读取前冻结；目前没有证据支持如此动态分配会改善真实亲和力排序，暂不自动切换分数/删水/重加权。

**Boltz-2 公平对照（后续阶段）**：原 87 例 Boltz-2 已计算并通过 `runs/fmo87/final/boltz_production_audit.json` 87/87 严格审计，无须假称未跑或重复消费；待因素与主评分冻结后，Galectin3 E36 与 FXa E40 的**各全部** 26/27 例需各自构建同靶点统一序列/MSA + ligand canonical SMILES 的无标签 Boltz-2 YAML，使用原 87 Boltz 对照相同 checkpoint、seed 0、recycling 3、affinity sampling 200、5 diffusion samples、`affinity_pred_value`、单靶点 Pearson/Kendall 及相同完整集合同例置信区间。Boltz-2 的输出是预测的 `log10(IC50 in µM)`，不是 FMO kcal/mol；只比较各靶点相关/排序，不直接相加或把输出当实验 ΔG；必须先核实 53 例的 MSA、蛋白序列、SMILES 标准化及原研究与模型训练数据潜在重叠。若 MSA 不可得不私自改用单序列模式而仍声称同协议。要区分 Boltz-2 推理结果与原 87 例论文已有基线，也不能按 FXa/Galectin3 的标签选择 Boltz 超参。

**E43 已做的迁移诊断**：读取 E36/E40 已审计的 GAMESS 原始 pair 日志与已公开的各 26/27 实验标签，事后复核 E41 的 ES+EX 等七条增量；不改变两项历史盲预测，不能冒充 E36/E40 预注册分项。`ES→ES+EX` 在内部 k10 Pearson `−0.276`，Galectin3 `+0.481`、FXa `−0.197`，方向相反，证实当前**不能**据 87 例结果自动按靶点添加/去除 EX。输出 `runs/fmo87/ablation/E43_external_pieda_marginals.json`，新脚本 `scripts/analyze_fmo87_external_pieda_marginals.py`。先维持 E42 对全靶点冻结的统一无拟合 ES+SOLV，未来若发现仅凭标签盲态结构特征选择模型的规则，必须新增独立标签冻结组验证，不得在这些已解盲靶点上训练选择器又当外部盲测。

**E44/E55 结合熵长轨迹门槛**：E11 CDK2/17 显式两态 3 重复 × 2 ns 的 free 状态有效样本数最低 `64.7`，split-half JS 最高 `0.412`，明显未达到门槛，绝非真实结合熵。E44 的 10 ns 运行因脚本中途修改导致协议哈希不一致，已由 `protocol_invalidated.json` 作废，4 条已完成轨迹仅作诊断。统一替代试验 E55 已完成同一 CDK2/17、3 重复 bound/free 各 `10,000,000` 步 × `1 fs = 10 ns`、100000 步平衡、每 5000 步取样、同一温度/力场/种子，并以 ETKDG+MMFF 独立初始化 free ligand；产物为 `runs/fmo87/ablation/E55_explicit_two_state_10ns_CDK2_etkdg_free/`。六条轨迹的身份/完整性检查通过，但 bound/free 最大 split-half JS=`0.412/0.289`，仍未通过 `0.15` 门槛，因此不扩至四靶点、不生成熵排名；torsional proxy 还**缺蛋白/水/平移转动/标准态熵**，不能加至 87 例当真实 `−TΔS_binding`。

**E45 无标签排名生产出分已完成**：`scripts/export_fmo_affinity_ranking.py` 在 E42 冻结的全部 87 个通过结构门槛的 case 上，逐个检查 E23 输入哈希、同一 AB 的 FMO2/DFTB3/PCM/选择性 PIEDA 模式、GAMESS 成功及所有 receptor–ligand pair 对数，直接复读原始日志输出 ES+SOLV（无标签、无拟合），所有目标全成功才生成带逐靶点名次的 `runs/fmo87/ablation/E45_label_free_ES_SOLV_ranking/predictions.csv` 和 SHA 审计。**独立的后验评估进程**将 87 个分数与实验标签按身份匹配，复得加权 Pearson `0.6336343767904575`、Kendall `0.5077501046195622`，与 E23 汇总精确一致。未来新队列可复用结构冻结与结果审计接口，但输入 PDB2PQR、固定蛋白 k10 最小化、同靶点 5 Å 并集与 GAMESS AB 的**实际运行仍需先完成**；不能因本次已有 E23 量子结果就宣称对任意新靶点的端到端自动化已验收。E11/E44 姿态仍用旧 E08 先导而非 E23 k10：**即使 E44 收敛，也不能直接与主评分拼成同几何总自由能**。

**外部 Boltz-2 53 例完整完成**：同一 ColabFold/MMseqs2 服务和 Boltz 官方 `compute_msa` 为 Galectin3 与 FXa 各生成一个共享 MSA，首行查询序列严格同原蛋白 A 链一致；E36/E40 标签剔除的 SDF 各转换 26/27 个 Boltz-2 affinity YAML，`runs/fmo87/boltz_external_inputs/{Galectin3,FXa}/manifest.json` 冻结所有身份、序列、MSA、输入 SHA，没有实验标签。`scripts/audit_fmo_external_boltz.py` 对全量预测、模型/checkpoint/协议和逐例 JSON/YAML/CSV 审计后才产出盲预测，独立 `compare` 阶段才合并实验值；结果在 `runs/fmo87/boltz_external_analysis/{Galectin3,FXa}/unblinded_comparison.json`。Galectin3 的 Boltz2 Pearson/Kendall 为 `0.800/0.647`，FXa 为 `0.678/0.479`。这不证明模型训练集不存在公开案例重叠，亦不使当前 FMO 主分自动外推成立。

**2026-10-01 并行执行状态**：E44 已作废；E55 六条轨迹已完成严格身份审计，但 bound/free 最大 split-half JS=`0.412/0.289`，未通过 `0.15` 收敛门槛，因此只记录为不收敛的配体扭转熵代理，不能进入 FMO 评分。外部 Galectin3 26 + FXa 27 已按冻结参数在 GPU0 完整执行；每组均先做盲态全量审计后解盲比较。其指标只作为外部对照，不用来后验调 FMO 因素，也不影响已经冻结的 FMO 评分。

**E44 单条轨迹的早期诊断，不是全两态结论**：已完成的 bound replica 1 单独逐扭转前后半程比较，五个扭转的 split-half Jensen–Shannon 距离中最高 `0.284`，lag-1 近似 ESS 最低 `85.6`（单条按原样 2000 帧统计），均提示该单条仍可能存在慢构象变化；最终门槛必须从六条完整轨迹**联合**按预定 `scripts/audit_fmo87_explicit_entropy_pilot.py` 重新计算，不能现在选择性终止/自称收敛，也不能把这些初步值当为真实结合熵。

**E44/E55 free-start 修订状态**：检查 E44 实现后发现其 free 模拟虽然使用游离体系拓扑，但初始坐标仍从 bound complex 中提取；随后 E44 又在运行中发生脚本修改，无法证明六条轨迹共享同一冻结协议，已由 `runs/fmo87/ablation/E44_explicit_two_state_10ns_CDK2/protocol_invalidated.json` 作废，已完成的 4 条仅作诊断。修订方案在 `scripts/run_fmo87_explicit_entropy_pilot.py`/`run_fmo87_explicit_entropy_matrix.py` 中使用冻结的 `--free-start etkdg`，实际执行 RDKit ETKDG 独立嵌入并用 MMFF 最多 500 步优化，再以同一原子顺序初始化 free ligand，并把 `free_start` 写进 provenance/审计。统一的 E55 六条轨迹已完成，但 bound/free split-half JS=`0.412/0.289` 未通过收敛门槛；该结果只用于说明当前熵代理不稳定，不与 E44 混合，也不称为完整结合熵。

**E46 水–配体原子接触审计完成 87/87**：原始 P38 34 例中 7 例水氧与配体重原子 `<2 Å`，最短 `0.921 Å`；E23 k10 固定受体局部最小化后 87 例全部 ≥`2.450 Å`。七个碰撞水分子 **7/7 在优化后仍保留且氧坐标原样**，对应配体全原子 RMSD 中位 `0.658 Å`；可归因于配体几何改动消除了这些严重接触，不能称“去水治好了它”。完整 87 个标签无关的原 PDB→SDF→k10 原子映射和水身份、几何结果：`runs/fmo87/ablation/E46_water_contact_geometry_audit.json`，脚本 `scripts/audit_fmo87_water_contacts.py`。对应方法的排序是否提高以已配对的 E23/E26 为准，不能从几何质控的成功推断预测指标因果改善。

**E47 FXa 外部靶点纯几何补充臂（FXa 已解盲，事后探索；完成）**：锁定 FXa E40 全部 27 个中性配体和同一个 PDB2PQR 受体、同一 E40 k10 口袋/片段/封端/PCM；仅将 OpenMM/OpenFF 位置约束 k10→k100、固定蛋白/晶水重算全部 27 个 AB。`scripts/analyze_fmo87_fxa_geometry.py` 核查 **27/27 受体 deck 完全相同**，重原子 RMSD 中位 `0.262 Å`，完整 tPIE Pearson `0.473→0.119`、Kendall `0.382→0.126`（配对 k100−k10 的 Pearson 95% CI `[−0.583,−0.136]`）；主 ES+SOLV Pearson `0.529→0.450`、Kendall `0.377→0.399`（Pearson CI `[−0.149,+0.012]`）。逐例 EX 改变量绝对值中位 `16.358 kcal/mol`、ES `2.972`，26/27 例 EX 改变量绝对值大于其他四分项的绝对改变量之和；因此**在这个对照里**完整 tPIE 的几何波动主要经 EX 传播，不能推断 EX 物理上无用。不能将已经解盲的 E47 当第二次独立外部盲验。原 `reference-pocket-poses` 硬编码问题已修复为读取 E40 `frozen_protocol.json` 的 `poses_root`，避免误用新姿态选口袋。完整分数/5000 次配对抽样、逐例 PIEDA/位移及审计证据：`runs/fmo87/ablation/E47_fxa_k100_pocket10/{unblinded_analysis.json,geometry_comparison.json,component_comparison.json}`。

**E48 原队列阈值稳健性完成**：E23 87/87 的三档 ES+SOLV 同协议逐例计分，四靶点分层配对 bootstrap 5000 次；k100/k25/k10 点估计均高于 Pearson 0.55、Kendall 0.38。主 k10 `0.634/0.508` 的条件抽样 CI 为 Pearson `[0.524,0.736]`、Kendall `[0.375,0.622]`，两条 CI 下界都略低于对应门槛；bootstrap 中 93.4% 的重复样本两项同时过线。去任意一个整靶点后的主评分点估计亦同时过线，但这仍只是原 87 化学系列的内部稳健性分析，非外部保证。`runs/fmo87/ablation/E48_threshold_robustness.json`、`scripts/analyze_fmo87_threshold_robustness.py`。

本文件的 `已完成` 仅指已完成的指定对照，**不等于所有可能的自由能因素都计算完毕**。以下项目仍需独立工程和计算资源；不能将 64/128 构象代理、单构象 QM 优化或 PBSA 跨模型求和称为它们的结果。

1. **完整两态结合熵（E11 扩展）**：为每个 case 的同参数化 AB 结合态和游离 A/B 建立显式溶剂/离子两态，逐态固定电荷/质子化、温度、标准态约定；同独立初速度/随机种子做至少三重复。分别报配体和蛋白构象贡献与平动/转动/标准态贡献的适用范围，不能只拿配体单态 `T·S` 充当 `−TΔS_binding`。先做四靶点覆盖的 pilot，记录失败、有效独立样本数、时间分块收敛和重复间区间；若先导轨迹不稳定，不允许扩展 87 例或声称全量熵已完成。当前虽有可用 OpenMM/CUDA 设备，但尚没有满足这些条件的成对长轨迹，不得编造熵值。
2. **统一 QM/PBSA 受体边界（E07 扩展）**：审计已量化 QM A 比 PBSA A 多 96/84/120/96 个 HOP/ACE/NME 相关原子；在能为全部边界新增原子提供可追踪电荷/半径、且两态电荷守恒前，不能擅自把 QM cap 坐标附加进经典 PBSA 或把现有两套能量视作同一个闭合循环。先构造一个 target/ligand 的闭合 pilot 并原子级复核，再进入 87 例全量。
3. **外部迁移（部分已执行，仍不具通用保障）**：E36 Galectin3 26/26、E40 FXa 27/27 已在各自实验标签解盲前锁定无拟合分项；一个靶点 ES+SOLV 为负、另一个为正，说明原 87 例内部达标不保证跨蛋白迁移。E39 FAAH 技术失败、未解盲；不是外部指标。若提出新的按结构/化学特征动态选模型规则，仍须在**新的、标签冻结的靶点**预先锁定分数、口袋、晶水、结构制备、标签映射和评价方式并全部汇报；不能以 E36/E40 已看结果训练后再称同一批的独立验证。Boltz-2 外部 26+27 已在各组先盲审后解盲：分别为 `0.800/0.647`、`0.678/0.479`，这是独立模型排序基线，不是 FMO 自动特征选择的证据；E55 仍待独立审计。
