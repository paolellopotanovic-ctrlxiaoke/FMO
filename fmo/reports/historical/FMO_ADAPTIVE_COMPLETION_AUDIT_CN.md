# 自适应 FMO 方法完成审计

日期：2026-10-02。本文区分“工程规则已实现”“当前数据已验证”和“仍需运行”，不把
后两者混为一谈。

| 要求 | 当前证据 | 状态 |
|---|---|---|
| 统一结构质控、补氢、冲突水处理、最小化 | `audit_fmo_external_structure_quality.py`、Galectin3/FXa `postrepair_audit.json` | 已实现；外部 target 仍需姿态方向验证 |
| 完整 FMO/PIEDA 候选集 | 原 87、E68–E70、E71 的 `unblinded_analysis.json` 和方法矩阵 | 已完成这些队列 |
| tPIE、ES、ES+SOLV、ES+CT+SOLV 独立保存 | `fmo_method_matrix.json`、各 target ensemble audit | 已完成 |
| 禁止嵌套项直接相加 | `adaptive_protocol.py` 的 `prohibited` 规则和单元测试 | 已实现 |
| 臂级结构质量门 | `Arm.structure_gate` 在部署判定中 fail-closed，失败臂不进入聚合 | 已实现并有回归测试 |
| robust scale 只作 rank-scale 诊断 | `scale_aggregation_audit.json`、`audit_fmo_adaptive_ensemble.py` | 已实现；未证明普适增益 |
| 多 restraint 数值复现性 | Galectin3 k1/k10/k100；`galectin_adaptive_deployment_v1.json` | 已完成；仅证明复现性 |
| 独立姿态方向验证 | Galectin3 26/26、FXa 26/27 通过独立 Boltz2 几何门；FXa/8m 单例需复核 | 结构姿态门已通过；不等于能量排序准确性 |
| 熵采样参数按结构复杂度预声明 | `entropy_protocol.py`、`e61_entropy_protocol_v1.json` | 规则已实现；多 target 尚未全部执行 |
| 熵 ESS/JS/block provenance 审计 | `audit_fmo87_entropy_blocks.py` | E61 第一轮 144/144 已完成但未通过；第二轮截至 2026-10-02 21:08 UTC 为 171/288，仍在追加 |
| 外部独立验证 | Boltz2 作为独立 fallback，不能和 FMO 相加 | 已有对照；不能替代 FMO 方向验证 |
| 全 bundle 一致性门 | `validate_fmo_adaptive_bundle.py` | 通过；当前活动目录为 E61 `171/288`，最终熵质量门仍警告 |
| 跨数据集因素/评分统一汇总 | `compile_fmo_universal_factor_summary.py`、`universal_factor_summary.json` | 已生成；24 个评分行、6 个 target 因素决策行，标签不参与选择 |

## 当前允许的生产结论

- 原 87 例：冻结 `ES+SOLV`，87/87 QM 完整；这是当前唯一已达到正式 FMO 生产门的队列。
- Galectin3：当前 E68 结构/QM 臂可部署，但能量准确性仍未建立；不能声称已建立跨靶点通用方向。
- FXa：当前 E92 是标准 union-pocket 结构触发臂；8m 保留 case-level review，能量准确性仍未建立。
- 熵：所有 target 默认 `OMIT`，直到 bound/free 采样和 provenance 审计通过。

## 不能宣称的内容

当前证据不足以宣称：

1. 一个无标签规则已经在所有 target 上达到同一 `0.55/0.38` 水平；
2. robust-scale、PCM、PBSA 或熵的叠加已经普适提升排序；
3. Galectin3 的反向结果已被结构修复完全解释；
4. E71 FXa 参考口袋结果就是正确的新 target 生产姿态。

这些不是代码缺陷，而是当前实验设计仍缺少的独立证据；流水线会保留 review/fallback，
而不是用解盲标签后验选择来掩盖缺口。
