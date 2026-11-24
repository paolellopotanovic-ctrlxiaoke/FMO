# 局部配体应变审计规则：实现记录（2026-10-25）

## 这次实际完成的内容

本轮把 `tPIE + 局部配体应变` 从一个可能改变排序的探索性分数，改造成主流程中的**无标签结构审计层**：

```text
输入姿态 → PDBFixer/PDB2PQR/图与水审计 → k10
        → repair_audit（逐配体、无实验标签）
        → 共享口袋/分片 → FMO/QM → 纯 tPIE
```

主流程新增 `repair_audit.json`，记录：

- 输入 SDF 与 k10 后配体坐标是否逐原子匹配；
- k10 后的几何质量门状态；
- 已登记的共轭 NH2 平面偏离及“输入正常、k10 后新变形”证据；
- 局部应变（如提供）仅作为诊断值，不作为主分数；
- 规则触发原因、输入哈希、候选修复状态和完整队列覆盖率。

## 规则边界

| 情况 | 当前动作 | 是否改主分数 |
|---|---|---:|
| 严重接触、缺失姿态或几何门失败 | fail-closed，要求注册修复或人工复核 | 否 |
| 仅有普通 warning 接触 | 记录 warning，不自动改坐标 | 否 |
| k10 诱发共轭 NH2 变形 | 生成实验候选标记；需要单独重跑 k10/分片/QM | 否 |
| 局部应变为负 | 检查能量参照与 QM 模型 | 否 |
| 局部应变较大但结构检查通过 | 保留数值并做分布审计 | 否 |
| 没有上述客观触发 | 纯 tPIE 继续作为主排序 | 否 |

高/负应变不能单独证明结构错误；否则会把“在某个集合上相关性更好”偷偷变成方法选择。

## 工具入口

审计已接入 `run`/`adaptive-run` 的 k10 后、QM 前阶段，也可单独审计已有完整产物：

```bash
fmotpiefinal repair-audit \
  --ligands-sdf RUN/full_target_fmo/ligands_no_properties.sdf \
  --minimized-root RUN/full_target_fmo/minimized \
  --geometry-report RUN/post_k10_geometry.json \
  --output RUN/repair_audit.json
```

只有显式追加 `--emit-experimental-candidates` 才会在旁路目录写出氢坐标候选。它不会覆盖原始 k10 结构；候选不是 QM 结果，必须重新 k10、重新分片、重新 QM，并完成完整队列回归后才可能进入默认流程。

## 已做的回归

在已有完整产物上运行了不改分数的审计回归：

| 队列 | 覆盖 | 共轭 NH2 候选 | 主分数 |
|---|---:|---:|---|
| BACE | 36/36 | 36 | 纯 tPIE |
| CDK8 | 10/10 | 0 | 纯 tPIE |
| Galectin3 | 26/26 | 0 | 纯 tPIE |
| FXa | 27/27 | 0 | 纯 tPIE |
| CDK2 | 16/16 | 0 | 纯 tPIE |
| JNK1 | 21/21 | 0 | 纯 tPIE |
| P38 | 34/34 | 0 | 纯 tPIE |
| TYK2 | 16/16 | 0 | 纯 tPIE |

BACE 的 36 个候选说明规则能够复现此前观察到的 BACE 几何信号；它不代表 BACE 的新 QM 排名已经自动替换生产结果。此前独立研究对照中，图触发候选的纯 tPIE Spearman 从约 0.177 到约 0.251，但仍属于研究证据，不能据此在其他 target 上自动加应变项。

## 验证状态

- 新增核心规则测试、候选只改被标记氢的测试和 CLI/队列覆盖测试；
- 当前发布包全套测试：`432 passed, 1 skipped`；
- `git diff --check` 和 Python 编译检查通过；
- 909-case 固定应变迁移队列未被修改、未被中断，也未用于选择规则；目前 `909/909` 都已处理，`892` 例严格验收成功，另有 17 例失败（FXa 2、PTP1B 11、OpenBind 4）。因此队列执行结束，但总体报告是 `INCOMPLETE`；FXa、PTP1B、OpenBind 不发布成功子集指标。

## 尚未声称完成的部分

共轭 NH2 候选尚未成为默认自动修复。要晋级，必须冻结触发条件，在旁路对完整 target 重跑 k10/QM，并对原 87、Galectin3、FXa、OpenBind 及已有集合做逐例回归，确认未触发样本坐标逐原子不变、失败全部保留、覆盖率完整，且不按实验亲和力选择规则。
