# OpenBind 晶体基准：原87例统一流程外测

## 范围与冻结规则

使用 `runs/openbind/openbind_fmo_input_manifest.json` 中全部 639 个独立晶体复合物、490 个化合物。预测之前不读取实验 pKD。每个晶体单独执行 `scripts/fmo_workflow.py` 的 `intake → geometry → adaptive-run`：先严格执行原87例的 PDBFixer、去除其氢、PDB2PQR/PROPKA pH7.4、原始结构/水审计、固定受体 OpenMM/OpenFF k10、5 Å 口袋、GAMESS FMO2-DFTB3/PCM PIEDA，输出 tPIE；**只在已分类的严格失败后**，按预先冻结的顺序触发质量审计和可追溯备用臂。不同晶体的独立受体不会伪装成同一靶点的共同受体；每个化合物若有多个晶体姿势，预先指定对所有晶体 tPIE 求算术平均。

**必要的无标签多链输入审计：** 不能固定选 A 链。检查 639 例原始 PDB 与对应晶体配体的重原子距离发现最近且小于 5 Å 的蛋白链分别为 A 523 例、C 105 例、B 11 例；不存在两个链同时在 5 Å 内或所有链都在 5 Å 外的样本。链选择已抽成通用 `scripts/audit_fmo_receptor_chains.py`；`fmo_workflow.py intake` **默认 `--chain auto`**，对所有同靶点配体逐一检查 `<5 Å` 唯一接触链并要求一致，显式指定错误链也直接拒绝。原87例的自动链为 CDK2 A/JNK1 A/P38 C/TYK2 A，与冻结输入一致。OpenBind 驱动还在 `frozen_inputs.json` 记录逐链距离以便后验核对。曾用固定 A 链启动的首批输出位于 `runs/openbind/new_fmo87_full_pipeline_20261015/`，因 116 例链选错已终止，**不能拿其成功子集计算外测指标**。修正后的目录为 `runs/openbind/new_fmo87_full_pipeline_chain_audit_20261015/`。

## 完整性与统计

仅当 639 个晶体全部经过结构与 QM 质量审计、涵盖全部 490 个化合物后，才运行 `scripts/evaluate_openbind_new_full_pipeline.py --allow-rescue`，计算 `−tPIE` 对实验 pKD 的 Pearson/Spearman/Kendall；该模式逐例审查哪几例改变参数。**不带 `--allow-rescue` 的严格模式仍拒绝出完整指标**。脚本核对输入 SHA、冻存链选择、实际 intake 链、结构处理、严格 k10、备用参数与其额外几何复核，并与以前完整基准做配对 bootstrap。任何仍失败的结构都不容许以成功子集冒充完整外测。

复算方式：

```bash
../fmo87-pipeline-tool/envs/fmo87/bin/python scripts/run_openbind_new_full_pipeline.py \
  --library-dir /inspire/ssd/tenant_predefaa-9a1b-4522-bb10-8850f313be13/global_user/3078-linziqian/ai-vs-physics-ross2023/envs/openfe/lib --workers 12 --ncores 4 --timeout 3600
../fmo87-pipeline-tool/envs/fmo87/bin/python scripts/evaluate_openbind_new_full_pipeline.py --allow-rescue
```

## 实跑结果与边界

严格重放已尝试 **639/639** 例，**635/639** 严格成功；原87例固定流程本身无法提供一个严格的 490/490 全覆盖指标。四例失败：`x2513a` 为 GAMESS 片段距离异常，`x4445a`、`x5836a`、`x6838a` 为单体/二聚体 SCF 不收敛。`scripts/evaluate_openbind_new_full_pipeline.py` 默认模式在首例失败时拒绝出完整指标。

为不漏掉四例，`scripts/fmo_workflow.py adaptive-run` 已自动调用通用 `scripts/rescue_fmo_new_target.py`；旧 `scripts/run_openbind_complete_pipeline.py` 的顺序补救经验仍保留，但**两套流程中最终采用的修复臂不同**：

| 晶体 | 旧流程最终修复 | 新统一流程按失败类别自动选择 |
|---|---|---|
| `x2513a` | 固定受体，配体约束 k0.1 | 重原子 5 Å 口袋；不改变 k10 姿态 |
| `x4445a` | 自由受体 k10 | 自由受体 k10 + DIIS |
| `x5836a` | 原始受体源 + 自由受体 k10 | 新 PDBFixer/PDB2PQR 受体 + 自由受体 k10、SOSCF |
| `x6838a` | 自由受体 k10 | 自由受体 k10 + DIIS |

新流程的规则是严格固定 k10 优先：`FAILED_GAMESS_RUNTIME` 先试重原子口袋；`FAILED_SCF_CONVERGENCE` 依次试 DIIS、增大迭代次数、自由受体 k10 + DIIS、自由受体 k10 + SOSCF，最后才试更弱配体约束；后 k10 几何不合格仅尝试自由受体几何修复。只按失败类别选择，不能看 pKD 或选最佳 tPIE。四例的 PDBFixer、PDB2PQR 输入及来源晶体仍来自这次新流程；备用计算各有独立目录、协议 SHA、失败尝试与几何审计。三个自由受体松弛体系无严重后处理接触，配体 RMSD <2 Å、最大位移 <4 Å；无法修复的新故障仍 `REVIEW_REQUIRED`。**含备用分支的 639/639 并非“全部采用原87例完全一致的固定 k10/QM 参数”**，单独命名为适应性完整结果，绝不能与严格重放混称。

| 490 个化合物，逐化合物对所有晶体取 tPIE 均值 | 晶体覆盖 | Pearson | Spearman | Kendall |
|---|---:|---:|---:|---:|
| 此前 OpenBind 完整适应性基准 | 639/639 | 0.373929 | 0.446519 | 0.292338 |
| 本次新结构流程 + 明确记录的 4 例 QM 补救 | 639/639 | **0.370716** | **0.443552** | **0.290418** |
| 新减旧 | — | −0.003213 | −0.002967 | −0.001920 |

评价方向统一为 `−tPIE` 与实验 pKD；相同 490 化合物配对自助抽样 1200 次得到新减旧的 95% 区间：Pearson `[−0.007521, +0.000216]`、Spearman `[−0.008962, +0.003068]`、Kendall `[−0.006691, +0.003096]`。差别很小，区间跨零，**没有证据表明新结构步骤提高 OpenBind 排序**。两个完整方法都已使用各自无标签补救分支，不能把小差异归因于链、氢、PDBFixer 或去水中的单独一项。OpenBind 是一个目标内部的化合物比较；原87例的四靶点加权相关不能与此处单靶点相关直接等同。

最终指标、逐化合物对照、四例备用尝试与严格失败记录分别位于 `runs/openbind/new_fmo87_full_pipeline_chain_audit_20261015/complete_490_adaptive_tpie.json`、`complete_490_adaptive_tpie.csv`、每例 `rescues/diagnostic_report.json`、`progress.json`。新的主命令一次运行即处理异常；旧单独批量补救脚本仅用于历史失败目录的重试和检查：

```bash
../fmo87-pipeline-tool/envs/fmo87/bin/python scripts/rescue_openbind_new_full_pipeline.py \
  --library-dir /inspire/ssd/tenant_predefaa-9a1b-4522-bb10-8850f313be13/global_user/3078-linziqian/ai-vs-physics-ross2023/envs/openfe/lib --workers 4 --ncores 4 --timeout 3600
../fmo87-pipeline-tool/envs/fmo87/bin/python scripts/evaluate_openbind_new_full_pipeline.py --allow-rescue
```
