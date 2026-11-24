# FMO 结构处理流水线

## 目标

把“某些结构处理后排序变好”的事后观察，固化为**标签无关、预先冻结、失败即复核**的结构质量流程。该流程只根据输入结构、化学一致性、几何接触和 QM 完整性做决定，不能读取实验亲和力，也不能根据相关系数选择处理方式。

## 固定决策顺序

1. **身份与化学检查**：配体 ID、SMILES、形式电荷、显式氢、3D 坐标、原子顺序和蛋白输入必须一致。
2. **原始接触审计**：检查配体与蛋白全部原子、晶体水及显式氢；不能只检查蛋白重原子。
3. **氢修复**：缺失或明显不一致时，保留配体重原子坐标，只重新生成显式氢坐标；不移动重原子。
4. **晶体水规则**：默认保留所有水。先在 target 的完整 cohort 上汇总严重配体–水碰撞；只要某个水残基被审计为严重碰撞，就把它加入该 target 的统一去水集合，并对该 target 的所有配体使用同一个受体版本。不能逐配体删不同的水，也不能因为排序变好而删除水。
5. **非水严重碰撞**：配体与蛋白非水原子严重重叠时，不自动移动蛋白或配体重原子，直接进入 `REVIEW`。
6. **最小化后复审**：统一力场、统一约束和统一停止条件重新最小化；若仍有严重接触或位移异常，进入 `REVIEW`。
7. **QM 完整性门**：全部预注册配体必须成功，输入哈希、片段身份、PIEDA 分项和结果状态必须闭合；否则不得输出完整排序。
8. **盲态冻结后才解盲**：结构处理计划、输入哈希、QM 参数和失败处理规则冻结后，才允许读取实验标签做回顾性验证。

## 自动决策脚本

先做结构审计：

```bash
$PY scripts/audit_fmo_external_structure_quality.py \
  --protein <prepared_receptor.pdb> \
  --ligands <ligands_no_properties.sdf> \
  --minimized-root <minimized_dir> \
  --output <structure_audit.json>
```

再生成标签无关处理计划：

```bash
$PY scripts/decide_fmo_structure_treatment.py \
  --audit <structure_audit.json> \
  --output <structure_treatment_plan.json>
```

决策脚本位于 `scripts/decide_fmo_structure_treatment.py`，其规则是：

- 只删除被严重配体–水碰撞直接指认的水，并把删除集合统一应用到整个 target；
- 只在缺氢时自动补氢；
- 不自动移动非水重原子；
- 非水严重碰撞、最小化后严重碰撞和结构不完整一律 `REVIEW`；
- 计划中的水残基若不在实际输入 PDB 中，执行器直接失败，不能按残基编号猜测或静默忽略；
- 任何结果都不读取实验标签。

## 当前证据

- P38 的事后去水消融使原 87-case tPIE 达到约 `0.556/0.395`，但这不是独立验证；它说明晶体水/结构质量可能是重要因素。
- Galectin3 修复水/氢后，tPIE 从原始的弱反向变为约 `0.24/0.17`，但 `ES+SOLV` 仍反向，说明结构修复能消除部分污染，不能保证整个能量指标正确。
- 当前 Galectin3 修复后审计：`25/26` 为 `PASS`，`18h` 因非水严重接触为 `REVIEW`；该 case 不应被自动修到通过。

## 生产输出规则

只有在结构计划不是 `REVIEW`、所有 QM case 完整、协议指纹一致、且姿态方向验证独立通过时，才允许进入 `FMO_READY`。否则输出 `FMO_REVIEW` 或 `FMO_FALLBACK`，并保留 tPIE、ES、ES+SOLV 和 ES+CT+SOLV 作为独立诊断分数，不能把它们相加。

## 几何约束与水处理联合 protocol

统一生成器 `scripts/build_fmo_structure_geometry_protocol.py` 固定生成：

- 几何臂：`k=100/25/10 kJ mol⁻¹ Å⁻²`；
- 正常 target：`wet × k100/k25/k10`；
- 源结构水几何门触发的 target：`wet/dry × k100/k25/k10`；
- 每个臂都输出 `tPIE`、`ES`、`ES+SOLV`、`ES+CT+SOLV`；
- 不允许逐 ligand 选择 k 或 wet/dry；
- `ES+SOLV` 是 tPIE 的嵌套子集，不能与 tPIE 相加。

原 87-case 当前自动计划位于 `runs/fmo87/ablation/final_analysis/fmo87_structure_geometry_protocol.json`：

- CDK2：保水，3 个几何臂；
- JNK1：保水，3 个几何臂；
- P38：源水几何门触发，wet/dry 共 6 个臂；
- TYK2：保水，3 个几何臂。

已有同协议探索结果显示，原 87-case 的保水 `k=10` `ES+SOLV` 为 `0.634/0.508`，P38 去水 `k=10` 完整 tPIE 为 `0.606/0.442`。这两项是已完成的历史臂结果，不等于新自动流程已经完成盲态重算；正式部署仍需保留输入哈希、统一受体和所有 arm 完整性记录。

## 协议完整性与结果覆盖

协议生成后必须运行 `scripts/audit_fmo_structure_geometry_protocol.py`。该审计器检查：

- 每个 target 是否按水处理政策注册了完整的 wet/dry 臂；
- 每个 target 是否同时注册 `k=100/25/10`；
- `tPIE`、`ES`、`ES+SOLV`、`ES+CT+SOLV` 是否全部保留；
- dry 臂是否具有 target-level 水残基清单；
- 原87结果矩阵是否完整，并明确标记为 `retrospective_audit_only`。

P38 的 target-level 选择性 dry19 试验现已完成：`34 cases × k100/k25/k10 = 102` 个 full-energy QM 作业和对应的 `102` 个 `MODMOL=1` selected-PIEDA 作业全部成功。full-energy 结果用于能量/输入交叉核验；PIEDA 分项只采用 selected-subsystem 输出。full FMO 日志的宽格式全片段 EX/CT 字段不用于打分，因为其全 pair 表解析会把巨大抵消项误当作配体分项。两种计算的基础 AB 输入哈希一致，所有配体 `tPIE` 差异不超过 `0.02 kcal/mol`。

新结果：`runs/fmo87/ablation/E87_P38_dry19_geometry_analysis/summary.json`。原87上，统一 k10 的 wet→dry19：tPIE 从 `0.533/0.412` 到 `0.601/0.439`；ES+SOLV 从 `0.634/0.508` 到 `0.595/0.463`；ES+CT+SOLV 从 `0.644/0.504` 到 `0.605/0.457`（Pearson/Kendall）。分层配对 bootstrap 显示 tPIE 的增益区间跨 0，而 ES+SOLV 的 Kendall 变化约为负、区间上界接近 0；因此当前支持“修复水碰撞改善 P38 tPIE，但不支持把全 target 去这19个水设为所有分数的默认策略”。wet 仍是默认，dry19 保留为 P38 的审计/候选支路。

因此，“协议审计通过”不等于“所有新臂都已经重新计算”。已有 E23/E22/E26 等结果只作为历史实验覆盖；缺少完整输出、结构指纹或协议指纹的臂必须保持 `FMO_REVIEW`，不能静默纳入生产排序。审计报告写入 `runs/fmo87/ablation/final_analysis/fmo_structure_geometry_protocol_audit.json`。
