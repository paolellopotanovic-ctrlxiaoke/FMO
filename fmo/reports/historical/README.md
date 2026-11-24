# Project documentation index

## Start here

1. `PROJECT_OVERVIEW_CN.md` — scientific question, exact cohort, final result, and interpretation limits.
2. `EXPERIMENT_PROTOCOL_CN.md` — frozen production protocol and quality gates.
3. `FINAL_REPORT_CN.md` — paper comparison, bootstrap uncertainty, target breakdown, and calibration ablation.
4. `REPRODUCTION_CN.md` — environment, external assets, end-to-end commands, runtime, and artifact regeneration.
5. `CODE_REFERENCE_CN.md` — the current, production-only file and tool inventory.
6. `PROJECT_HANDOFF_CN.md` — current status, paper differences, sensitivity results, and next commands.
7. `FMO_SENSITIVITY_CN.md` — FMO protocol sensitivity definitions and evaluation commands.
8. `ANALYSIS_PLAN_CN.md` — complete hypothesis matrix, stop criteria, and next-stage execution order.
9. `OPENFE_ABFE_ANALYSIS_CN.md` — Boltz-2 中 OpenFE/ABFE 的论文、数据集、开源边界和复现路线。
10. `SOPHOSQM_CLAIMS_AUDIT_CN.md` — 对 SophosQM 70-case 拟合精度与独立预测可用性的逐条审计。
11. `AFFINITY_CALIBRATION_VALIDATION_CN.md` — 文献对照、冻结参数和新化合物/新几何的保留集验证协议。
12. `FMO_TOTAL_ENERGY_RANKING_CN.md` — 完整 FMO 总能量双重作差（无需拟合的靶点内排序）和新流水线使用方法。
13. `FMO_TOTAL_ENERGY_RESULTS_CN.md` — 87/87 无拟合总能量排序、原始日志审计与 tPIE/Boltz-2 全量对照。
14. `FMO_OFFICIAL_BINDING_AUDIT_CN.md` — GAMESS 官方子体系逐数值验证、PCM 分量与 87/87 游离配体优化审计。
15. `TANKYRASE_FMO_SENSITIVITY_CN.md` — Tankyrase 2 MFMO 文献逐项对照与新三项能量差敏感性。
16. `TANKYRASE_PBSA_FOLLOWUP_CN.md` — Tankyrase 2 式独立 PBSA 三态计算、回顾性排序与非等价边界。
17. `TANKYRASE_ORIGINAL_REPRODUCTION_CN.md` — 原论文 23 例原始评分、旋转键与 PCM 消融逐数值复现；从头等价量化计算缺失输入审计。
18. `TANKYRASE_INSPIRED_GAS_FMO_PBSA_CN.md` — 自主气相三态 FMO + PBSA 计算、计分审计、试点失败与非等价边界。
19. `TANKYRASE_GAS_PBSA_87_RESULTS_CN.md` — 全量气相 12/87 收敛阻碍、独立 PBSA 87/87 和非等价 PCM+PBSA 87-case 指标。
20. `ANEWFE_VALIDATION_AUDIT_CN.md` — AnewFEP 的直接物理预测、回顾性协议改进与本项目 FMO 端点分数的边界审计。
21. `FMO_ADAPTIVE_PIPELINE_CN.md` — 结构质控、姿态多臂、PIEDA 分支和熵收敛门组成的自适应生产流程。
22. `FMO_ADAPTIVE_PRODUCTION_DECISION_CN.md` — 可部署的结构/姿态/PIEDA/熵判定树、fallback 和不确定性状态。
23. `ANEWFE_VALIDATION_AUDIT_CN.md` — AnewFEP 的方法边界与本流水线映射。
24. `FMO_ADAPTIVE_METHOD_V2_CN.md` — 结构门、分支选择、尺度聚合和熵收敛的可迁移执行规范。
24. `FMO_FINAL_ADAPTIVE_METHOD_CN.md` — 当前证据、固定处理顺序、判定规则和跨数据集交付总结。
25. `FMO_ADAPTIVE_RUNBOOK_CN.md` — 新 target 的结构门、QM coverage、姿态、分支、尺度和熵执行顺序。

## Implementation and audit

- `ARCHITECTURE_CN.md` — module responsibilities and the FMO/Boltz-2 data paths.
- `METHOD_ANALYSIS_CN.md` — SophosQM method contract and differences from the Boltz-2 in-house FMO baseline.
- `GAMESS_SETUP_CN.md` — GAMESS installation boundary, DFTB3 settings, fragmentation, and failure gates.
- `COMPLETION_AUDIT_CN.md` — 87/87 production checks and summary statistics.
- `PROVENANCE.md` — external datasets, papers, software versions, hashes, and replacement components.
- `ARTIFACT_MANIFEST.md` — which files are intended for git and which are generated locally.
- `REFERENCES_CN.md` — primary literature and tool references.

## Status

- `STATUS_CN.md` is the current completion snapshot.

Label-free 外部化学/QM 一致性审计入口：`scripts/audit_fmo_external_chemistry.py`。
完整 target/结构臂/评分分支矩阵：`scripts/compile_fmo_method_matrix.py`。
自适应分支的独立准确性压力测试：`scripts/audit_fmo_adaptive_selection.py`。
无标签协议判定器：`scripts/derive_fmo_adaptive_protocol.py`；物理因子决策表：
`scripts/compile_fmo_factor_decision_table.py`；尺度聚合消融：
`scripts/analyze_fmo_scale_aggregation.py`；最终 target profile：
`scripts/compile_fmo_target_profiles.py`。
