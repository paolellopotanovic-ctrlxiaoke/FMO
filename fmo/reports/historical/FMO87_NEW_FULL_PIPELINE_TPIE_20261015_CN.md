# 原87例：新统一端到端流程实跑（仅 tPIE）

日期：2026-10-15。这次是**实际重跑**，不是复用旧 dry19 分数。输入为冻结 `data/fmo87_cohort.csv` 中 CDK2 16、JNK1 21、P38 34、TYK2 16 例；`scripts/prepare_fmo87_full_pipeline_inputs.py` 按其配体 ID/电荷/立体化学逐个核对并生成不含实验标签的各靶点 SDF 与 CSV，JNK1 源 SDF 中非原87的17例不会进入。输入 SHA 和 ID 记录在 `runs/fmo87/new_full_pipeline_20261015/input/input_manifest.json`。

**同一套执行步骤**：四靶点各自 `fmo_workflow.py intake → geometry → run`；来源结构按全靶点原始 SDF 审计严重水接触，同靶点所有配体使用一个去水方案；筛选受体链 → PDBFixer → 清除其氢 → PDB2PQR/PROPKA pH7.4 → 每步复核 → OpenMM/OpenFF 固定受体 k10 → 5 Å 全配体并集口袋 → GAMESS FMO2-DFTB3/PCM 选择性 PIEDA，输出跨受体–配体片段 `tPIE`。k10 消除原始 CDK2 部分非水冲突；各靶点 post-k10 零严重接触；四靶点全队列 QM 状态 `BLIND_FMO_COMPLETE_REQUIRES_EXTERNAL_VALIDATION`。运行配置每靶点 `--workers 4 --ncores 4 --timeout 3600`；四靶点并发。四靶点水规则相同，不代表删水数量相同：CDK2/JNK1/TYK2 原始 PDB 无晶水、删 0；P38 原 PDB 有 343 水、直接严重碰撞并集删 **6**。这**不是旧 P38 删19水的 dry19 协议**，不可套用旧数据的分数。

运行结果（每靶点内部比较 `−tPIE` 与 `−实验 ΔG`；下方加权是按例数，不是把不同靶点的原始能量混池求相关）：

| 靶点 | QM成功 | Pearson | Spearman | Kendall |
|---|---:|---:|---:|---:|
| CDK2 | 16/16 | 0.8745 | 0.8706 | 0.7000 |
| JNK1 | 21/21 | 0.6165 | 0.6093 | 0.4439 |
| P38 | 34/34 | 0.5038 | 0.4722 | 0.3214 |
| TYK2 | 16/16 | 0.7078 | 0.7147 | 0.5333 |
| **全87按靶点例数加权** | **87/87** | **0.6367** | **0.6232** | **0.4596** |

旧 87 `dry19 + k10` 对照为 `Pearson 0.6007 / Spearman 0.6022 / Kendall 0.4392`；本次全流程对应增加 `+0.0360/+0.0209/+0.0204`。逐靶点配对 case-bootstrap 1800 次的 Pearson 差 95% 区间为 **`[−0.0117,+0.1038]`**，Kendall 差 `[−0.0289,+0.0780]`：**不能宣称已经证明显著更好**。旧新之间同时改变了 PDBFixer、氢/蛋白处理和 P38 去水数量，不能将变化单独归因于删水或修氢。原 87 是开发数据；达到阈值 0.55/0.38 是回顾性验证，不能声称跨靶点前瞻性保底。

复算和审计：

```bash
../fmo87-pipeline-tool/envs/fmo87/bin/python scripts/evaluate_fmo87_new_full_pipeline.py
```

`scripts/evaluate_fmo87_new_full_pipeline.py` 校验 87 例身份、所有靶点盲态完整成功、几何状态与 SDF SHA、固定蛋白 k10 元数据和指定 QM 协议；只有全部通过才读取实验标签并与旧案例逐例配对。最终逐例结果及输入 SHA：`runs/fmo87/new_full_pipeline_20261015/complete_87_tpie.csv`、`runs/fmo87/new_full_pipeline_20261015/complete_87_tpie.json`。每个靶点 QM 原始日志和结构审计都留在 `runs/fmo87/new_full_pipeline_20261015/{CDK2,JNK1,P38,TYK2}/`。
