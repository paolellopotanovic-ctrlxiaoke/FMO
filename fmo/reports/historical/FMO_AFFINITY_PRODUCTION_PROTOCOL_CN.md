# 标签隔离的 FMO 亲和力排序生产协议

状态：原 87 例的固定评分/审计/排序已验证；Galectin3/FXa 的 Boltz-2 同例外部比较已完成端到端验收；E55 显式溶剂熵诊断已完成但未通过收敛门槛。日期：2026-10-01。

## 输入和决策规则

输入是一组有相同受体、同一测定系列的配体 3D SDF 与一个原始蛋白 PDB；在看该组实验亲和力前固定所有分子、质子化、晶水政策、序列、口袋/配体几何及评价口径。`scripts/freeze_fmo_affinity_protocol.py` 接受一个按 `data/fmo87_cohort.csv` 结构列组织的队列 CSV；可以包含实验标签列，但**代码只投影白名单结构字段且输出中的结构指纹不含标签属性**。入口检查重复配体、SDF/SMILES 身份、非中性、缺失显式 H/三维坐标、蛋白链缺失和来源姿态严重重原子碰撞。失败标记 `REVIEW`，不能丢掉失败分子再按幸存者报告全队列排名。每个靶点独立固定整个队列的协议，**不依据 QM 分数或已知 ΔG 为某一配体改用去水、EX、口袋半径或 SCF 结果更好的方案**。

冻结后的主计算为 PDB2PQR AMBER/pH7.4、保留晶水且蛋白/水位置固定、配体 OpenMM/OpenFF restrained minimization `k=10 kJ mol⁻¹ Å⁻²`、靶点**全配体姿态并集** 5 Å 口袋、统一 HOP/ACE/NME cap/天然二硫键和连续肽键审核、同一参数 GAMESS FMO2-DFTB3/PCM 选择性 PIEDA。排序分数直接取**同一 AB 计算**的配体–蛋白跨片段 `Σ(ES+SOLV)`，单位 kcal/mol，系数都为 1；越负代表预计结合越强。这里 `SOLV` 已包含于完整 tPIE，不能在完整 tPIE 上重复加。不是独立 A/B/AB 三态差、不是以最小化后的单体能量/结合熵代替标准态 ΔG。

### 新队列：先诊断，再全组统一计算

此处的“判断适合什么计算手段”是**仅凭输入结构/计算质量的可行性分流**，不是先看部分配体实验排序，再为目标/配体调出最高相关的参数。确定队列和评价口径后按以下次序执行，所有修订记录版本与哈希：

| 门槛 | 只看标签之外的观测量 | 整组统一处置 |
|---|---|---|
| 入口身份/化学 | SDF 名称/SMILES、显式氢、三维、净电荷、蛋白链和重原子接触；运行 `freeze_fmo_affinity_protocol.py` | 所有 case `PASS` 才继续；目前仅冻结**中性**配体协议，带电、错误结构或不同蛋白/姿态来源需整组另立协议或停下，不能只弃失败例 |
| 蛋白化学 | PDB2PQR 状态、可解析残基/晶水、共价二硫键、编号缺口的肽键 C–N 几何、跨片段异常碰撞、QM 片段电子奇偶性 | 结构原因未解决则停下；修复后对靶点**整个队列**重制蛋白/封端/5 Å 并集，不沿用旧受体下的幸存者 |
| 姿态/水 | 固定受体 k10 最小化后的配体–蛋白/水接触、配体位移、跨配体相同受体坐标和同一口袋 | 若出现严重碰撞/大位移，先查质子化、水位和原始姿态；审计或整组重做，不对单个配体切 k100、去水或放大口袋 |
| QM 与排名 | 全部 AB 的收敛、PIEDA 原始分项、输入 SHA/同靶点受体一致性；运行 `export_fmo_affinity_ranking.py` | 只在**全组成功**后输出一个固定 ES+SOLV 排名；失败即停组/报告覆盖率，不用只成功样本计算 Pearson/Kendall |

ES、EX、完整 tPIE、PBSA、三态差和熵是**同组对照/诊断项**，并非这个生产评分的逐分子“智能”开关；当前 Galectin3 与 FXa 的相反外部趋势也不足以推导无标签、可验证的靶点专用选择器。未来若要引入规则 `g(结构质控)→评分协议`，需先固定无标签输入、阈值、整组一致执行和失败策略，再在**未用于发现规则**的新靶点验证；不能因为 87 例点估计达标就宣称这种自适应规则已验证。生产分数亦不等于有标准态/熵的绝对结合自由能。

E54 又将同 AB 的五项 PIEDA **全部 31 种非空 0/1 组合**在原 87 例三档及两个外部靶点完成逐例评分；在 15 种非空背景中加入 EX，原 87 k10 为 **0/15** 提升 Pearson，而 Galectin3 为 **14/15**。这是重要的因素交互/外推风险而非可供当前队列挑选“最佳组合”的预注册证据。完整探索性表格及 9734 行逐例分数在 `runs/fmo87/ablation/E54_full_pieda_factorial.json`、`runs/fmo87/ablation/E54_full_pieda_factorial_cases.csv`；生产主分数仍固定 ES+SOLV，不从这些已解盲组合中为 Galectin3 后验选 EX 特例。

E54 的更严格汇总是：原 87 k10 的 31 个组合中 8 个同时过 `0.55/0.38`，但 Galectin3 与 FXa 各 31 个组合均为 0 个同时过线。因而当前证据支持“原87上 ES/CT/SOLV 的组合有效、EX 对原87不利且跨靶点交互明显”，不支持任何已知的、无需新标签验证即可自动迁移的因素开关。

对 E36 Galectin3 和 E40 FXa 的两个既有外部队列再用**同一无标签结构准入**回审：`scripts/audit_fmo_external_structure_preflight.py` 检查冻结的原始蛋白、原始和脱敏 SDF 哈希，仅按 SDF 结构生成必要的队列列并调用生产 `freeze()`；`runs/fmo87/ablation/E50_external_structure_preflight/label_free_preflight.json` 为 **53/53 通过**（26+27）。这说明两个队列在该**基础入口检查**下均可计算，**绝不**证明 Galectin3 的 ES+SOLV 倒序可以由这些门槛事先发现或修复：此类门槛判断“能不能按相同协议算”，不是判断“算出来能不能准确排序”。

基础入口后、开始实际碎片生成前，用 `scripts/audit_fmo_affinity_protein.py --policy ... --output ...` 对整个靶点原蛋白检查共价几何：CYS/CYX `SG–SG≤2.3 Å`、跨非连续 PDB 编号及插入码的实测肽键 `C–N 1.15–1.65 Å`；对应报告给出 **全靶点** `merge_disulfides`/`repair_peptide_numbering` 候选开关，而非逐配体调参。E52 的原 87 四个源蛋白均未检出这两类边，Galectin3 亦无；FXa 有 4 对二硫键（含先前 QM 失败的 CYS191–CYS220）、2 处跳号但真实成键的肽键（含 GLY218–CYS220）和 7 处插入码连续肽键，对应 E40 实际统一使用的两个结构修订。若同一个编号与**同名残基**只靠插入码区分、而解析器可能丢弃插入码，则标记人工复查，不能自动认定安全。该原 PDB 结构建议必须在 PDB2PQR 后、实际并集口袋上重新审核，尤其不能切断跨口袋二硫键；质量诊断不等于通过亲和力排序外推。

外部已计算队列还由 `scripts/audit_fmo_external_prepared_protein.py` 对照 E52 原始蛋白和 E36/E40 冻结准备蛋白、全 53 例量子盲态受体哈希：FXa 4 个二硫键、2 处跳号肽键和 7 处插入码肽键在 PDB2PQR 后距离与身份保持不变；Galectin3 无此类边，两个靶点的共价/整组碎片开关一致。见 `runs/fmo87/ablation/E53_prepared_protein_covalent_audit.json`。新靶点仍须复核最终**实际口袋**没有截断任何共价键，不能只凭原 PDB 或准备后报告跳过 FMO deck 逐例审计。

## 运行及审核

```bash
export PY=/path/to/fmo87/bin/python
export PYTHONPATH=src:scripts
$PY scripts/freeze_fmo_affinity_protocol.py \
  --cohort data/fmo87_cohort.csv \
  --output runs/fmo87/ablation/E42_frozen_production_policy.json

# E23 已有同口袋、同姿态的全部 87 个 QM 日志；新队列必须先按冻结政策
# 运行 PDB2PQR→整组固定受体最小化→共享口袋→同参数 GAMESS AB，
# 并将完整成功结果整理为 QM_ROOT/<target>/<ligand>/input.inp,run/...
$PY scripts/export_fmo_affinity_ranking.py \
  --policy runs/fmo87/ablation/E42_frozen_production_policy.json \
  --qm runs/fmo87/ablation/E23_fixed_geometry_k10_pieda \
  --output runs/fmo87/ablation/E45_label_free_ES_SOLV_ranking
```

第二步重新解析**全部病例**的 GAMESS 原始 pair 表、检查冻结的蛋白与 SDF 结构指纹、同靶点逐原子受体坐标相同、输入 SHA/FMO2/DFTB3/PCM/选择性 PIEDA/片段数/原始 tPIE 与结果文件一致；任一例失败则不生成部分 `predictions.csv`。成功输出独立靶点名次、原始分数和逐例 log/input SHA，**不读取实验 ΔG**；随后使用单独评价脚本/流程加入实验数据评价每靶点 Pearson/Kendall，跨靶点汇总以每靶点样本数加权，配对差值按靶点分层 bootstrap。原 87 例 87/87 通过，独立后验标签评价为 `0.633634/0.507750`。上述两脚本冻结/审核出分已经可复用，但 E23 QM 上游执行依赖现有多阶段脚本，未宣称任意新体系都可无人监管自动收敛。遭遇二硫键/肽键缺号等真实结构错误，只能在**没有目标标签**时修复原因、对整组统一重算并保留失败日志。

## 影响因素与何时停下来

| 因素 | 原 87 例同输入/对照结论 | 外部风险与处理 |
|---|---|---|
| 弱约束 k100→k10 | 完整 tPIE `0.308→0.533`，ES+SOLV `0.624→0.634` | E47 FXa 同受体 27 例 k100→k10 完整 tPIE Pearson `0.119→0.473`，ES+SOLV `0.450→0.529`（后者配对区间跨零）；固定 k10，不根据新集标签切换 |
| ES 基线上加 SOLV | k10 ΔPearson `+0.068`，95% CI `[−0.016,+0.137]` | Galectin3 ES+SOLV 为负，FXa 为正；PCM 贡献非通用收益保证 |
| ES 基线上加 EX | k10 ΔPearson `−0.276`，CI `[−0.459,−0.107]` | Galectin3 `+0.481`，FXa `−0.197`；不根据已解盲靶点制定特例切换器 |
| MMFF 同姿态应变 | E35 k10 的 ES+SOLV `0.634→0.313` | 与 QM/PCM 的哈密顿量不一致；不加进冻结主分数 |
| PBSA/PB polar | E25/E28 对照跨界项破坏 87 例排名 | 绝不把 PB 与 AB pair SOLV 当同一可独立增加的溶剂因子 |
| 原始晶水重原子短接触 | E46 的 P38 7/34 例原始 <2 Å，最短 0.921 Å；同一晶水保留且不移动，k10 后 87/87 水距 ≥2.450 Å | 在配体最小化后复查水坐标/距离；不能按原始碰撞为个别分子删水或改评分 |
| 全结合熵 | 两态 E11 2 ns 的 ESS/JS 门槛未通过；E44 作废；E55 六轨迹完成但 bound/free split-half JS=`0.412/0.289` 均未过 `0.15` | 缺蛋白/水/标准态项；不按不成熟代理改变已达标主排名 |

详见 `docs/FMO87_ABLATION_RESULTS_2026-09-30_CN.md` E23/E33–E45 的定义、逐例数据、配对区间与外部反例。这里“无拟合”只指评分的物理分项以原始系数 1 相加，**不宣称** k10/ES+SOLV 在最初 87 个实验值被观测前就已经预注册。

### `tPIE` 与 `ES+SOLV` 的嵌套关系

当前 GAMESS PIEDA pair 表满足：

```text
tPIE = ES + EX + CT + DI + SOLV
ES+SOLV = ES + SOLV
```

所以 `ES+SOLV` 是 `tPIE` 的子集，不是可以直接加到 `tPIE` 上的独立物理贡献。禁止使用 `tPIE + ES+SOLV` 作为“加入更多物理项”的无拟合能量，因为它会重复计算 `ES` 和 `SOLV`。如果需要保留两条信息，只能把它们作为两个候选分支，或明确写成混合分数 `Sλ=(1−λ)tPIE+λ(ES+SOLV)`，并记录 λ 是预先固定的无标签规则还是监督拟合参数。熵、PBSA 或 SA 只有在独立收敛/身份审计通过后，才可能作为这两条 PIEDA 分支之外的增量。

原 87 例的 `0.633634/0.507750` 是**该队列的点估计**，并非任意未来队列的保底指标。E48 靶点分层 5000 次重采样，Pearson 95% 区间约 `[0.524,0.736]`、Kendall `[0.375,0.622]`，两者同时超过 `0.55/0.38` 的重采样比例 `0.9342`；区间下界并未同时达标。不能据此保证迁移成功，亦不能将外部 Galectin3 的失败掩去。

E51 单例删去敏感性检查逐一重算同一 87 例，`87/87` 个删去版本仍双指标过线（最差 Pearson `0.613`、Kendall `0.495`）。这是**同一已观察数据内**的离群点稳健性，不是新靶点上准确率的保证。

## 后续的 Boltz-2 同例比较

原 87 例 Boltz-2 已跑完并通过 `runs/fmo87/final/boltz_production_audit.json` 的 87/87 严格检查。对外部 Galectin3 26 例和 FXa 27 例，**在因素结论和对照协议冻结后**，按现有 `src/sophosqm_baseline/boltz_inputs.py` 方式生成共享蛋白序列/同一 MSA 策略的无标签 YAML；两靶点完整运行 Boltz-2 `seed=0`, recycling=3、affinity steps=200、affinity diffusion samples=5。Galectin3 的 FMO/ Boltz2 指标分别为 `-0.538/-0.362`、`0.800/0.647`；FXa 分别为 `0.529/0.377`、`0.678/0.479`。结果文件分别为 `runs/fmo87/boltz_external_analysis/{Galectin3,FXa}/unblinded_comparison.json`。仅采用 Boltz 的 `affinity_pred_value` 进行靶点内排序/与同一实验值比较，不将其数值与 FMO kcal/mol 直接相加；公开模型潜在训练集重叠仍未被排除。

**输入和推理已完成**：从原公开蛋白 A 链提取 Galectin3/FXa 长度 137/233 的序列，使用与原 87 例相同的 Boltz `compute_msa`/ColabFold MMseqs2 策略生成每靶点一个非单序列 MSA，验证 CSV 首行查询序列与蛋白逐字符一致。MSA 与产物分别位于 `runs/fmo87/boltz_external_msa/` 和 `runs/fmo87/boltz_external_inputs/{Galectin3,FXa}/`；manifest 含 source/MSA/YAML 哈希，没有亲和力标签。修复后的 `scripts/run_boltz_fmo87.py` 会在运行后检查所有合法 Boltz 输出目录，避免把已成功生成的 affinity JSON 误判为 runtime failure；`scripts/audit_fmo_external_boltz.py` 的 `audit`/`compare` 两阶段均已对两个外部靶点通过。**未读模型训练索引，不能保证 53 个公开案例未进入 Boltz 训练集**；后续解释必须注明此类潜在重叠。

以下是已经执行的同例比较复现命令示例（`TARGET=Galectin3`；FXa 需相应更换三个路径）：

原 E44 在 free replica 2 启动前后发生了脚本修改，无法再证明六条轨迹共享同一冻结协议；因此已在 `runs/fmo87/ablation/E44_explicit_two_state_10ns_CDK2/protocol_invalidated.json` 明确作废，已完成的 4 条仅作诊断，**不进入熵分析或排名**。E55（free ligand 用 ETKDG 独立初始构象）六条轨迹已完成并通过身份/完整性审计，但 bound/free 最大 split-half JS 为 `0.412/0.289`，超过 `0.15` 收敛门槛；因此 E55 是“轨迹完整但扭转代理不收敛”的负结果，绝不误称全结合熵，也不加入生产分数。外部 Boltz-2 已在同一冻结 GPU/模型参数下完成，其 26+27 例均在整组无标签审计后才解盲。

```bash
TARGET=Galectin3
HOLDOUT=runs/fmo87/ablation/E36_galectin3_holdout
SOURCE=vendor/public_binding_free_energy_benchmark/fep_benchmark_inputs/structure_inputs/misc/galectin3_extra_ligands.sdf
PY=/path/to/boltz-enabled/python
PYTHONPATH=src:scripts "$PY" scripts/run_boltz_fmo87.py \
  --manifest "runs/fmo87/boltz_external_inputs/$TARGET/manifest.json" \
  --out-dir "runs/fmo87/boltz_external_predictions/$TARGET" \
  --python "$PY" --cache /path/to/boltz/cache --accelerator gpu --devices 1 \
  --seed 0 --recycling-steps 3 --sampling-steps-affinity 200 \
  --diffusion-samples-affinity 5
PYTHONPATH=src:scripts "$PY" scripts/audit_fmo_external_boltz.py --stage audit \
  --manifest "runs/fmo87/boltz_external_inputs/$TARGET/manifest.json" \
  --predictions-root "runs/fmo87/boltz_external_predictions/$TARGET" \
  --holdout "$HOLDOUT" --source "$SOURCE" \
  --output "runs/fmo87/boltz_external_analysis/$TARGET"
PYTHONPATH=src:scripts "$PY" scripts/audit_fmo_external_boltz.py --stage compare \
  --manifest "runs/fmo87/boltz_external_inputs/$TARGET/manifest.json" \
  --predictions-root "runs/fmo87/boltz_external_predictions/$TARGET" \
  --holdout "$HOLDOUT" --source "$SOURCE" \
  --output "runs/fmo87/boltz_external_analysis/$TARGET"
```
