# 通用 FMO 工作流：可运行的入口、三态/溶剂核查与仍待实现的生产环节

日期：2026-10-13。目标：建立**无需按标签挑 case 的结构/计算入口**，继续分析三态和溶剂为何不能稳定改善排序，将跨片段相互作用解析扩展到蛋白–蛋白、蛋白–RNA 和 RNA–配体；最终对未知靶点的全量配体也可完整运行。**现状是分阶段原型，不是已完成的跨相互作用生产流水线或已证实高 Pearson 的万能修正。**

## 1. 分数能否统一修正：封闭的原87三态与溶剂对照

运行 `python scripts/audit_fmo_three_state_solvent_transfer.py`。输出 `runs/fmo87/ablation/E17_three_state_solvent_closed_audit_20261013/audit.json`、`cases.csv`、`gas_12_case_diagnostic.csv`。使用旧 E17 **共享受体、同姿态 FMO2-DFTB3/PCM** 87 例和其 E04 独立 A/B/AB、E15 native subsystem、E07 PBSA 结果；逐例重新读取 A/B/AB 原始 `energy_hartree`、native 各项和 PBSA 原值，并强制闭合/覆盖。与 2026-10-13 的 E23 wet k10、E87 dry19 或 OpenBind 的输入**不是相同协议**，数值不可横向冒充同一条验证臂。按 CDK2/JNK1/P38/TYK2 内算 Pearson/Kendall 并按例数加权；能量与实验 ΔG 都取“越负越好”，不反转单边符号。

| 同一 E17 体系的无拟合候选 | Pearson | Kendall | 相对 tPIE 的含义 |
|---|---:|---:|---|
| PCM `tPIE` / native `DEint` | 0.308 | 0.246 | AB 跨片段直接作用；已验证 native 与 selected PIEDA 相符 |
| PIEDA `ES` | 0.555 | 0.421 | 仅静电，**是 tPIE 子项，不是独立附加自由能** |
| PIEDA `ES+SOLV` | **0.624** | **0.487** | 有标签的历史开发对照分支；不能后验声称新靶点已验证 |
| Native `DEsolv` 单项 | −0.486 | −0.403 | 仅数学相关方向；具体物理误差来源未单独定位 |
| Native `DEint+DEsolv` | 0.054 | 0.071 | 比直接作用下降 0.254，配对 bootstrap 95% 区间约 `[−0.394,−0.122]` |
| Native `DEint+DEPLd` | 0.232 | 0.227 | 比直接作用下降 0.076，区间约 `[−0.128,−0.027]` |
| Native/独立三态结合能 | −0.018 | 0.043 | 已逐例数值闭合；相对作用项下降约 0.325 |
| PBSA 极性项 | −0.479 | −0.352 | 不可视为 PCM 分解的同一溶剂项 |
| `PCM tPIE + PBSA polar` | 0.030 | 0.078 | **跨模型拼接的诊断代理，不是闭合热力学自由能**，Pearson 比 PCM tPIE 下降 0.278 |

这些是**当前 E17 体系内的回顾性证据**，不是物理上证明去溶剂化“不存在/应删除”。`ES+SOLV` 不能再加一次 tPIE 或把 PCM 的 `SOLV` 与 PBSA 极性项直接等价；输入/电荷/半径/帽片定义不同，溶剂响应还可能重复计数。只有 12/87 成功的 gas PIEDA（11 例 CDK2、1 例 P38），CDK2 11 例 Pearson gas `0.109`、matched PCM `0.456`；**不能称 87 例气相比较，也不能凭 11 例固定气相修正**。

与当前 Pearson 尺度研究 `docs/FMO_PEARSON_OPTIMIZATION_DIAGNOSIS_20261013_CN.md` 共同表明：一个物理上较完整的三态表达式未必预测排序更好；OpenBind 上冻结截尾/MW 校正也未跨集稳定。**因此暂未发现已证明可“统一固化并显著提高所有集合 Pearson”的方法。**现有生产分支仍应分别报告原始 tPIE、PIEDA 子项、独立三态和独立溶剂诊断，不改动历史冻结排序定义。

### 同 deck PIEDA 溶剂消融补充（2026-10-14）

`../fmo87-pipeline-tool/envs/fmo87/bin/python scripts/audit_fmo_pieda_solvent_transfer.py` 生成 `runs/fmo87/ablation/pieda_solvent_transfer_20261014/evaluation.json` 和逐例 `scores.csv`。严格检查两套原 87 例 87/87 身份及 P38 dry19 之外各分数不变，外部 Galectin3 26/26、FXa 27/27 的 QM 完成状态和输入 SHA；每个分项来自**同一 AB deck 的 PCM PIEDA**，`tpie_without_solv = tpie - SOLV`，而不是 `tpie + SOLV` 再把已含的 SOLV 加第二次。统一取能量与实验负号为“越大越有利”；原 87 例以下四靶点内算相关再按例数加权，外部为单一靶点。配对 bootstrap 在每个靶点内重采样 1800 次；是回顾性条件区间，不是前瞻性泛化证明。

| 数据和协议 | n | tPIE Pearson | 去 SOLV Pearson | ES Pearson | SOLV Pearson | ES+SOLV Pearson | ES+CT+SOLV Pearson |
|---|---:|---:|---:|---:|---:|---:|---:|
| 原87 wet k10 | 87 | 0.533 | 0.432 | 0.566 | 0.168 | 0.634 | 0.644 |
| 原87 P38 dry19 k10，其余 wet | 87 | 0.601 | 0.552 | 0.524 | 0.107 | 0.595 | 0.605 |
| Galectin3 自动结构 + k10 | 26 | 0.559 | 0.390 | −0.532 | 0.462 | −0.218 | −0.086 |
| FXa 自动结构 + k10 | 27 | 0.586 | 0.647 | 0.647 | −0.155 | 0.617 | 0.593 |

原87 dry19 的 `tPIE/ES+SOLV/ES+CT+SOLV` 加权 Spearman 分别为 `0.602/0.612/0.612`，Kendall `0.439/0.463/0.457`；Galectin3 分别为 `0.550/−0.151/−0.022`（Spearman），FXa 为 `0.563/0.522/0.517`。wet 87 的 `ES+SOLV − tPIE` Pearson 差 `+0.100`，配对 95% 区间 `[−0.050,+0.284]`；dry19 差 `−0.006`、`[−0.075,+0.062]`；Galectin3 差 `−0.777`、`[−1.175,−0.261]`；FXa 差 `+0.031`、`[−0.069,+0.129]`。删除 SOLV 相对 tPIE 的差在 wet 87 是 `−0.102 [−0.168,−0.022]`，Galectin3 `−0.169 [−0.412,−0.018]`，FXa `+0.060 [−0.028,+0.159]`。

**结论是两种固定规则都不能据此推广为普适 Pearson 改进**：wet87 的 ES+SOLV 增益随 P38 去水处理消失，Galectin3 相同定义甚至显著反向；去掉 tPIE 自带溶剂贡献也明显伤及 wet87、Galectin3。原87每个靶点的相关和所有逐例分数在输出文件中。原87的重跑 k10、外部的结构和受体制备不完全同源，Galectin3/FXa 也早已参与方法开发；这些是**跨既有协议稳健性排查，不是盲测外部验证**。E17 的 native A/B/AB 和 PBSA/PCM 结果仍属另一套几何协议，禁止与本节 k10 分项做逐项加和或当作同 deck 独立变量。现阶段保留 tPIE 为默认分数，分项只作标明协议的诊断；后续要在未看标签的独立系列冻结处理决策后再验统一修正。

## 2. 新目标 + SMILES 的可执行流程（仅输入与结构门）

统一命令 `scripts/fmo_workflow.py`：

```bash
python scripts/fmo_workflow.py intake --target NEW --ligands new_ligands.csv --output runs/new/intake
python scripts/fmo_workflow.py intake --target NEW --ligands new_ligands.csv --protein bound_receptor.pdb --poses aligned_bound_ligands.sdf --output runs/new/intake
python scripts/fmo_workflow.py geometry --intake-dir runs/new/intake --output runs/new/geometry
python scripts/fmo_workflow.py run --geometry-dir runs/new/geometry --output runs/new/full --gamess-dir software/gamess --library-dir OPENFE_LIB --dftb-parameter-dir software/gamess/auxdata/DFTB/3OB-3-1 --dry-run
python scripts/fmo_workflow.py run --geometry-dir runs/new/geometry --output runs/new/full --gamess-dir software/gamess --library-dir OPENFE_LIB --dftb-parameter-dir software/gamess/auxdata/DFTB/3OB-3-1 --workers 1 --ncores 4
```

`new_ligands.csv` 为 `ligand_id,smiles`，不接收含实验效力字段的输入。只给靶点名/分子式或 SMILES **不包含三维结合构象、蛋白结构、核酸结合位点、质子化和水的信息**：入口明确返回 `NEEDS_STRUCTURES`，不能直接计算亲和力；“分子式”还不唯一确定分子图，必须得到明确 SMILES/分子结构。配套输入需要**同一受体参考坐标系**的蛋白 PDB 和每个配体的晶体/外部无标签模型结合构象 SDF；先核对配体身份、形式电荷、显式 H、几何和源哈希。

`geometry` 复用既有 `audit_fmo_external_structure_quality.py`、`decide_fmo_structure_treatment.py` 和 `apply_fmo_structure_treatment.py`：按预定义 severe contact 检测整个靶点的冲突水并**统一移除**；未修复的非水严重碰撞标记 `K10_REPAIR_REQUIRED`，必须在统一 k10 后全部复检通过才能进入 QM；水残留、身份/原子异常直接停在 review。水处理仅修改受体，不移动配体；结构通过仍保持 `ranking_eligible=false`。

`run` 是蛋白–配体执行器 `scripts/run_fmo_new_target.py`：只选 cohort 注册的受体链，PDBFixer 修重原子/质子几何，**清除 PDBFixer 氢（尤其 HXT），让 PDB2PQR/PROPKA 重新放氢**，每一步重新检查严重接触；统一 k10 最小化、完整 GAMESS FMO2-PCM 和 tPIE/ES/SOLV/ES+SOLV 导出。可选 `--repair-peptide-numbering` 必须在同一靶点全例一致，所选参数写入状态及 frozen QM 协议；针对 FXa 的片段间虚拟帽碰撞需审计。任何一步失败写 `REVIEW_INCOMPLETE` 而非部分排序，全部同源 ligand 名称成功才输出 `blind_predictions.json`，状态 `BLIND_FMO_COMPLETE_REQUIRES_EXTERNAL_VALIDATION`。仍缺多构象/带电态采样与不同靶点广泛验证。

实验性 `--pocket-heavy-only` 保持选中**完整受体残基**，仅用重原子距离选 5 Å 队列并集；默认历史含氢口袋不变。针对 FXa 原始 PDB 27/27 可消除 GLN61 所在口袋臂的 SCF 故障，但 FXa Pearson 从历史 E92 的 0.622 降至 0.545；Galectin3 相同结构流水线 26/26 从 0.556 降至 0.421。**不能把它自动设为通用“高精度”回退方案**。`pocket_selection_audit.json` 解释每个差异残基的最近重/全原子距离和直接接触或填补断口来源；全量对比见 `docs/FMO_PEARSON_MECHANISM_AND_STRUCTURAL_REPAIR_20261013_CN.md`。

**实测端到端冒烟（不是新的泛化指标）**：原87 CDK2/17 的晶体姿态从单独的新 target 输入经过上述新 CLI 的 intake→geometry→PDBFixer→PDB2PQR→k10→GAMESS，1/1 完成，k10 后严重碰撞 0，日志/QM 成功。隔离输出 `/tmp/fmo-goal-single-case-20261013/qm_chainA_noH/run_status.json`，tPIE `−71.786 kcal/mol`、ES+SOLV `−25.969 kcal/mol`；**单例分数不能计算 Pearson，也不等同原87的全队列共享口袋协议。**首次冒烟错误暴露了 HXT 与 PDB2PQR 的冲突，已在独立氢重建边界修复并重测成功。

## 3. 非蛋白–小分子相互作用的真实可移植层

```bash
python scripts/fmo_workflow.py partner-intake --complex bound_complex.pdb --interface-type protein_rna --chain-a A --chain-b B --output runs/new/partner_intake.json
python scripts/fmo_workflow.py interface --log gamess.log --partition fragments.json --output runs/new/interface.json
```

`partner-intake` 支持 `protein_protein`、`protein_rna`、`rna_ligand`；RNA–小分子额外要求 `--ligand-sdf bound_ligand.sdf` 并检查重原子元素、数量和与 PDB 的坐标匹配。不从序列、配体 SMILES 或孤立结构**虚构结合构象**；配对结构需预对齐、链不重叠、重原子接触合理。未识别非标准残基、缺失 SDF、强重叠或核酸化学不确定则拒绝。

`interface` 使用 `src/sophosqm_baseline/interface_pieda.py` 按手工或经审计的 **GAMESS fragment ID 两组** 汇总跨组 PIEDA，仅计 A–B，不把 A 内部或 B 内部对能算进来；要求每个跨组 fragment pair 恰好一条、各分项可闭合，日志正常结束，必要时与已有蛋白–小分子 `result.json` 和 GAMESS native selected TIE 对齐。`fragments.json` 例：`{"fragment_ids_a": [1,2], "fragment_ids_b": [3,4], "interface_type": "protein_rna"}`。已有蛋白–小分子的 CDK2、Galectin3、FXa 日志已实测与生产 tPIE 闭合；蛋白–蛋白、蛋白–RNA、RNA–小分子目前仅对解析/分组做了**合成结构测试**，还**没有可复现的生产 GAMESS AB deck、fragment-cap/核酸参数、SCF 和外部亲和力验证**。更不能宣称其 tPIE 是蛋白–蛋白/RNA 结合自由能。保留完整体系统一参数、交互对适应性和类型门控是后续开发要求。

## 4. 达成用户最终目标还缺的工程和证据

| 工作流阶段 | 当前证据 | 完成门及尚需动作 |
|---|---|---|
| 统一 Pearson 修正 | OpenBind/MW/percentile/三态/溶剂全量回顾审计；跨集有负结果 | 在**全新 target/系列**冻结 label-free 规则盲测，对比 raw tPIE、MW-only、Boltz2 与不确定性；未证明则维持 raw 而非篡改指标 |
| 三态/溶剂物理误差来源 | E17 87/87 解析和原始 A/B/AB、native、PBSA 数值闭合 | 冻结一致蛋白/配体几何与电荷、PCM/PBSA 溶剂参考，逐 target 非重复计数对照；补齐气相目前失败的 75 例或明确不选用气相完整比较 |
| 新靶点纯 SMILES → 结合 pose | 无 pose 时工具明确停机 | 外部无标签结构/构象生成、位点指定、手性/互变异构审核、多构象统一聚合和失败门；不能由文字 target 名默认产生晶体级证据 |
| 配体结构进入完整 QM 生产 | CDK2/17 单例成功；Galectin3 新流程 26/26 成功、回顾 Pearson 0.556；FXa 在历史已审计受体加统一肽编号规则下 27/27 成功、Pearson 0.586（旧 E92 0.622） | 原始 FXa 仍需自动识别 GLN61/片段化故障并重建；多构象/质子微态采样和**从未见标签**的靶点验证仍待实证 |
| PPI/protein–RNA/RNA–ligand | 新类型的 bound 结构身份门和通用 PIEDA pair-score；三类仅合成测试 | 体系适用的片段构造/核酸帽/质子化/离子水、GAMESS AB 作业、独立接口参考、完整同类数据集比较和资源预算；未获证据前绝不声称已泛化预测性能 |

本表表示目标仍在推进，当前工作并**未达到“完整完成并在各类 interaction 上稳定改进 Pearson”**。软件验收必须包括失败样本不静默忽略、不同协议不混合、统一 target 受体和水计划、所需物理项闭合、原样本数保留、全新 target 盲态指标；此前报告的相关系数不替代这些验收项。
