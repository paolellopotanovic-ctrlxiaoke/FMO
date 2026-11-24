# FMO 自适应方法 v2：可迁移的无标签执行规范

日期：2026-10-02。本文定义的是**执行和报告规则**，不是用当前数据集标签
训练一个新的回归器。它借鉴 AnewFE/AnewFEP 的核心工程思想：输入结构、状态
定义、采样和失败模式必须分层审计；不能把一个后验上看起来最好的分数伪装成
普适物理规律。

## 1. 分数对象的严格定义

对每个配体 `i`，FMO/PIEDA 产生

\[
  tPIE_i = ES_i + EX_i + CT_i + DI_i + SOLV_i.
\]

`ES+SOLV`、`ES` 和 `ES+CT+SOLV` 是同一 PIEDA 分解中的嵌套分支，不是彼此
独立的自由能项。因此生产排序只能选择其中一个分支，不能写成
`tPIE + ES+SOLV`。

PCM、PBSA、三态差值和 native subsystem decomposition 也不能自动相加。只有
在 bound/unbound 状态、参考几何、介电边界、标准态和符号完全一致时，才可作为
独立的自由能代理；否则它们保留为候选分支和诊断量。

## 2. 无标签决策树

### 2.1 输入结构门

1. 检查蛋白全部原子（包括水氢）、配体显式氢、价态、电荷、元素和重原子键连。
2. 任一硬碰撞：统一删除冲突水、重新补氢、保留配体重原子姿态，并对整个 target
   重新最小化；不能只修复排序中“有问题”的配体。
3. 修复后重新生成所有口袋、fragment、封端和 FMO deck，并重新做 QM coverage。
4. 结构修复、系列异质性或位移超阈值时，运行同一输入规则下的 `k=1/10/100`
   restraint arms；不能逐配体选择不同 `k`。

### 2.2 分数选择门

定义候选集合

\[
  \mathcal C = \{tPIE, ES, ES+SOLV, ES+CT+SOLV\}.
\]

- 全部 QM 成功、结构门通过，但没有独立方向验证：生产候选固定为全局冻结的
  `ES+SOLV`；其它分支只作审计。
- 跨 restraint 的 Kendall、rank spread 和逐例波动只测**数值复现性**，不测准确性。
- 只有预先冻结的独立新靶点/盲测证据验证了某分支方向，才能把该分支写入新的
  production profile；已解盲数据上的最佳分支只能标记 `review_candidate_score`。
- 结构或 QM 不完整：不输出完整 FMO 排名，状态为 `FMO_REVIEW` 或
  `FMO_FALLBACK`。

当前实现中，`derive_fmo_adaptive_protocol.py` 强制执行这一规则，避免 Galectin3
这类“跨 arm 很稳定但整体反向”的分支被自动部署。

## 3. 数值尺度处理

若预先声明需要组合多个**非嵌套候选排序**，只能在 target 内进行 robust rank-scale：

\[
 z_j(i)=\frac{S_j(i)-\operatorname{median}(S_j)}
 {1.4826\operatorname{MAD}(S_j)},\qquad
 z_{\rm agg}(i)=m^{-1}\sum_j z_j(i).
\]

这表示排序尺度聚合，不表示 kcal/mol 相加。它只能作为冻结协议中的独立候选，
不能使用实验标签选择权重，也不能用系列内归一化掩盖 Galectin3 的状态/姿态问题。

## 4. 熵项的正确加入条件

对于同一 ligand 的 bound/free 状态，只有在以下条件全部满足后，才允许计算

\[
 \Delta G_{\rm entropy}=-T(S_{bound}-S_{free}).
\]

条件包括：相同力场、质子化和原子映射；至少两个独立 replica；完整轨迹；有效
样本数达标；split-half JS 和 block/replica 重复性通过；并且明确声明这里只是
某些构象自由度的熵，不是包含蛋白、溶剂、平移、转动和标准态的完整结合熵。

若只得到 torsion proxy、某一条轨迹或未收敛 block，熵项必须为 `OMIT`。不能把
“不收敛”通过调大采样权重强行改写为收敛；正确处理是追加预先声明的 block，预算
耗尽仍不通过就保持 `OMIT`。

## 5. 当前各 target 的无标签处理

| target | 结构/QM状态 | 默认候选 | 额外处理 | 发布状态 |
|---|---|---|---|---|
| 原87（CDK2/JNK1/P38/TYK2） | 统一结构、87/87 QM 完整 | `ES+SOLV` | 保留 tPIE/PIEDA/PCM/PBSA/几何消融审计 | `FMO_READY` |
| Galectin3 | E68/E69/E70 统一修复，26/26×3 QM 完整 | `ES+SOLV` | 维持三 arm；系列异质性 review；禁止后验选择 `ES+CT+SOLV` | `FMO_REVIEW` |
| FXa | E71 27/27 QM 完整，但使用旧 E40 参考口袋诊断 | `ES+SOLV` | 标准新 union-pocket 需独立复算；参考口袋结果不可直接冻结 | `FMO_REVIEW` |

原 87 的冻结结果为 Pearson `0.633634`、Kendall `0.507750`。Galectin3 的
反向分支不能靠加权解决；FXa 的 E71 证明修复结构/QM 可闭合，但把失败定位到
新姿态/union-pocket 拓扑，不能把诊断臂冒充生产臂。

## 6. 迁移与验证边界

真正的普适性不能由当前 87 例内部后验选择证明。发布新 target 时必须保留：

- 输入和结构修复 SHA256；
- 全部 QM case 和失败原文；
- 所有 restraint/solvent/PIEDA 分支；
- 结构、姿态、化学异质性和熵 audit；
- 一个未参与分支选择的外部或 prospective 验证集。

如果外部 target 的所有 FMO 分支方向都未被独立验证，系统应诚实输出
`FMO_REVIEW`/`FMO_FALLBACK`，并把 Boltz2 或其它模型作为独立结果，绝不与 FMO
kcal/mol 直接相加。这个降级不是失败，而是避免把结构/状态误差误报为量子化学
精度。

生产接口因此区分两个字段：`default_score` 是 FMO 审计候选，
`recommended_output` 是当前可发布的结果来源。`FMO_REVIEW` target 可以推荐
`INDEPENDENT_FALLBACK_REVIEW`，但 fallback 不得和 FMO 分数做线性相加。

## 7. 可执行的多臂部署门（2026-10-02）

`src/sophosqm_baseline/adaptive_protocol.py` 和
`scripts/build_fmo_adaptive_deployment.py` 将上述规则固化为机器可审计的入口。
它与 AnewFE/AnewFEP 的可迁移思想对应如下：

| AnewFE/AnewFEP 的问题 | FMO 端对应的硬规则 |
|---|---|
| 输入结构/参数化错误 | 先通过 `effective_quality_gate`，否则只能修复或降级 |
| 采样或构象陷阱 | 要求统一 `k=1/10/100` 臂并报告逐分数秩一致性 |
| 统计收敛不等于准确性 | `Kendall >= 0.80` 只表示数值复现；`orientation_verified` 仍必须独立通过 |
| 单个 restraint/结构臂失效 | `Arm.structure_gate` 逐臂 fail-closed；不能用中位数聚合掩盖坏臂 |
| 自由能状态不闭合 | bound/free 熵未通过 ESS、split-half JS 和身份审计时强制 `OMIT` |
| 回顾性挑选模型 | 不读取实验标签；`ES+SOLV` 是冻结默认分支，其余只作审计 |

对每个 target，发布条件现在是：

\[
  \text{FMO\_READY} = Q_{structure}\land Q_{QM}\land
  Q_{orientation}\land Q_{score\ stability}.
\]

其中 `Q_score_stability` 只描述不同计算臂的数值一致性，不能替代
`Q_orientation`。因此 Galectin3/FXa 即使 `ES+SOLV` 的跨臂 Kendall 通过，因姿态
方向尚未独立验证仍输出 `FMO_REVIEW`；这是故意的安全门，而不是把当前数据后验
解释成泛化能力。

物理分数和诊断聚合严格分离：

- `median_arm`：同一物理分数在多个统一臂之间取中位数，单位仍是 kcal/mol；只用于
  数值异常值抑制，不能混合不同 PIEDA 项。
- `robust_z_mean`：每个分数先按 target 内 median/MAD 标准化，再取均值；单位是
  rank scale，只用于诊断，不可标成能量，也不可与 `tPIE` 相加。
- `tPIE`、`ES+SOLV`、PBSA、PCM 和熵仍是不同候选对象；未经状态闭合和独立验证不得
  线性相加。

生成部署审计示例：

```bash
PYTHONPATH=src python scripts/build_fmo_adaptive_deployment.py \
  --target Galectin3 \
  --arm k1 E69/unblinded_analysis.json \
  --arm k10 E68/unblinded_analysis.json \
  --arm k100 E70/unblinded_analysis.json \
  --structure E68/postrepair_audit.json \
  --chemistry galectin_chemistry_audit.json \
  --output galectin_adaptive_deployment_v1.json
```

该入口的输出是部署状态、默认分支、结构/姿态/化学/熵门、各臂秩一致性和下一步
动作，不包含实验标签选择逻辑。

## 8. 熵采样参数的结构自适应规划

`sophosqm_baseline.entropy_protocol.plan_entropy_protocol` 根据体系原子数、配体可旋转键数和配体重原子数计算预声明复杂度指数：

\[
 C = N_{rot}/12 + N_{heavy}/35 + N_{system}/100000.
\]

它只改变采样预算，不改变任何能量权重：

| `C` | replicas | blocks/state | block 长度 | 适用解释 |
|---:|---:|---:|---:|---|
| `<=2` | 3 | 3 | 10 ns | 低复杂度诊断 |
| `2–4` | 3 | 4 | 15 ns | 中等复杂度 |
| `>4` | 4 | 4 | 20 ns | 高复杂度/大体系 |

原复杂度公式只是一条预声明的起始采样启发式，不是已经验证的收敛预测器。`gpu-hours`
不能单独推导能跑多少模拟纳秒；只有输入同一体系、力场、设备和采样设置下测得的
`--calibrated-hours-per-block`，才报告计划 GPU 小时数。没有实测校准时标记
`UNESTIMATED_NO_RUNTIME_CALIBRATION`；实测预算不足时标记
`INSUFFICIENT_FOR_PREDECLARED_PLAN`，不静默减少 replica/block。后续追加采样只由
未通过的预声明 ESS、split-half JS、replica/block 稳定门触发，不按标签调整；不通过
时熵继续 `OMIT`。该规划只针对有限 torsional
bound/free proxy，不能声称包含完整结合自由能的蛋白、溶剂、平移、转动和标准态熵。

当前 E61 的 `e61_entropy_protocol_v1.json` 是按 CDK2/17 的 54,236 原子溶剂化体系、
8 个可旋转键和 32 个配体重原子生成的中复杂度计划；正在运行的 18-block recovery
任务仍按其已冻结的 block 协议完成，不会中途改参数。
