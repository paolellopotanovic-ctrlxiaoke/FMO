# FMO 自适应亲和力排序流水线

日期：2026-10-01。本文把当前原 87 例、Galectin3 和 FXa 的消融结果整理成一套可重复的生产流程。目标不是用实验标签为每个新靶点后验挑分数，而是先用结构和数值质量决定需要哪些计算，再在预先冻结的候选分支中输出结果。

## 1. 先区分两种模式

### 1.1 零拟合物理排序

用于真正未知的目标，不能看该目标实验亲和力：

```text
结构审计 → 必要修复 → 整组统一姿态生成 → AB FMO2/PCM/PIEDA
→ tPIE、ES+SOLV、ES+CT+SOLV 并行诊断
→ 多姿态一致性审计 → 输出一个预先定义的分支或“不确定”
```

`tPIE` 是跨蛋白–配体 pair 的完整 PIEDA 总项：

```text
tPIE = ES + EX + CT + DI + SOLV
ES+SOLV = ES + SOLV
```

因此 `ES+SOLV` 是 `tPIE` 的子集，不能与 `tPIE` 直接相加。熵、PBSA 和 SA 也不能在没有重复计数审计时直接叠加。

### 1.2 有标签校准排序

如果新靶点已经有一批实验数据，可以拟合 `tPIE`、ES、SOLV、clogP 等特征，但必须报告：

- 同靶点训练集结果；
- 同靶点 leave-one-out 或 scaffold split 结果；
- 完全冻结参数在新几何/新配体上的结果。

同靶点拟合不应伪装成无拟合物理预测，也不能把原 87 例上的点估计当作未知靶点保底性能。

## 2. 结构质量门

旧门只检查蛋白重原子，漏掉了晶体水氢与配体的硬重叠。现在所有新目标必须运行：

```bash
python scripts/audit_fmo_external_structure_quality.py \
  --protein <prepared.pdb> \
  --ligands <ligands.sdf> \
  --minimized-root <minimized> \
  --output <audit.json>
```

审计对象包括：

1. SDF 名称、SMILES、净电荷、显式氢、三维构象和重原子坐标；
2. 蛋白重原子、蛋白氢、晶体水氧和水氢；
3. 配体–受体最短距离、vdW 比、异常重复坐标和水冲突；
4. 最小化前后势能、能量下降、配体重原子位移；
5. 质子化、PDB2PQR 残基、二硫键/肽链连接和 QM 片段身份。

审计报告现在同时保留两个门：`quality_gate` 是原始输入坐标门，
`effective_quality_gate` 是实际最小化复合物门。若原始 SDF 仍有可由最小化
消除的近接触，但 `minimized_complex.pdb` 中没有硬碰撞，报告为
`PASS_AFTER_MINIMIZATION`；这不是掩盖原始问题，而是明确记录“问题发生在
输入姿态、还是仍存在于实际 FMO 姿态”。

门槛的物理含义不是把普通氢键误判为碰撞：氢–重原子接触使用硬距离门，重原子–重原子使用 vdW 比；普通约 `1.5 Å` 的氢键不作硬失败，但 `0.7 Å` 的配体氢–水氢/氧接触必须修复。

### 2.1 结构修复顺序

整个 target 队列必须使用同一修复政策，不能只修复表现差的 ligand：

1. **水冲突**：删除与任一队列配体发生硬冲突的完整水残基，写入新受体哈希；不逐 ligand 删除不同水。
2. **显式氢异常**：保留重原子姿态，删除并按重原子几何重新生成所有配体氢；脚本为 `scripts/repair_fmo_external_structures.py`。
3. **重原子碰撞**：不能自动把单个 ligand 平移到看起来更好的位置；应重新确认公共 docking/crystal pose，或整组进入姿态重建协议。
4. **大能量/大位移**：触发 `k=1,10,100` 三档整组姿态对照，不能只为某个 case 改参数。

5. **系列尺度差异**：对外部队列额外记录每个化学系列的净电荷、重原子数、
   分子量、氢键供受体和可旋转键。若系列平均重原子数相差至少 15%，整组进入
   `FMO_REVIEW`，增加姿态/采样审计；这个门只改变审计强度，禁止按系列切换
   分数、做 heavy-atom 归一化或使用实验标签选择分支。

Galectin3 的修复就是这个流程的实例：发现 `HOH403/410/475/493` 冲突水，重新生成配体氢；重算前 25/26 个 case 的最小化初始能量恢复到约 `-6000 kJ/mol`，原先的百万级异常消失。`18h` 仍被标为重原子近接 review case，不能静默删除。

## 3. 姿态与 FMO 分支

### 3.1 无标签的量纲与取消诊断

对每个 target 先在全队列上计算以下只依赖结构/能量的统计量：

\[
 r_{\mathrm{cancel}}=\frac{|EX|+|CT|+|DI|+|SOLV|}{|ES|+\epsilon},\qquad
 r_{\mathrm{solv}}=\frac{\operatorname{MAD}(SOLV)}{\operatorname{MAD}(ES)+\epsilon}.
\]

它们只用于判断“附加项是否主导排序或发生强烈抵消”，不直接决定哪个分数更
接近实验。若不同 arm 的同一分量 rank Kendall 低于预设门槛，说明该分量对
姿态敏感，应报告不确定性并保留所有候选；不能用实验标签后验选择分支。

这一步的数学边界很重要：Pearson/Kendall 对正比例和正平移不变，所以单纯
z-score、平移或单位换算不可能修复排序；只有经过预声明的非线性物理变换、
状态差值或真正的采样平均才可能改变排序。任何拟合权重必须进入独立的
训练/验证协议，不能冒充零拟合 FMO。

只要结构门触发 review，就对整组生成不同 restraint 的 blind arms：

| arm | 用途 |
|---|---|
| `k=100` | 接近原始姿态，诊断“姿态移动是否造成排序变化” |
| `k=10` | 默认折中姿态 |
| `k=1` | 允许更充分的局部松弛，诊断局部极小问题 |

三条 arm 必须保持同一受体、同一水政策、同一 QM 口袋、同一 FMO 参数；只允许改变预先声明的配体 restraint。然后运行：

```bash
python scripts/analyze_fmo_pose_ensemble.py \
  --arm k1=<dir> --arm k10=<dir> --arm k100=<dir> \
  --output pose_selection.json
```

新 target 可以先用 `scripts/plan_fmo_adaptive_protocol.py` 生成无标签计划。
结构修复或审阅时，计划自动要求 `k=1/10/100` 三臂；只有结构门直接通过时
才允许单独 `k=10`。计划还会把熵标为 `ELIGIBLE` 或 `OMIT`，不会因为熵文件
存在就自动把它加入分数。计划中的 `READY_TO_RUN` 只表示“可以按冻结协议
开始计算”，不表示亲和力方向已经验证；无标签方向始终记录为
`UNVERIFIED_LABEL_FREE`，完成新 target 后不能跳过独立验证或 fallback 判定。
熵延长则使用 `scripts/plan_fmo_entropy_sampling.py`；当前 E55 的正式计划为
追加 `20 ns`、bound/free 各 3 个独立 replica，结果写入
`runs/fmo87/ablation/E55_explicit_two_state_10ns_CDK2_etkdg_free/next_sampling_plan.json`。

这个工具只用不同 arm 的预测排序一致性，不读取实验标签。它用于识别数值不稳定和生成 consensus ranking，**不能保证物理分数方向正确**；Galectin3 已经展示了“排序很稳定但方向仍反向”的情况，所以稳定性门不是准确性证明。

## 4. 当前数据的实际结果

| 数据集/协议 | tPIE | ES | ES+SOLV | ES+CT+SOLV |
|---|---:|---:|---:|---:|
| 原 87，k10 | 0.533/0.412 | 0.566/0.477 | **0.634/0.508** | 0.644/0.504 |
| Galectin3 原 E36 | 0.131/-0.127 | -0.581/-0.424 | -0.538/-0.362 | -0.456/-0.313 |
| Galectin3 修复 k10 | **0.241/0.170** | -0.491/-0.449 | -0.432/-0.378 | -0.413/-0.307 |
| Galectin3 修复 k100 | 0.200/0.238 | -0.559/-0.498 | -0.499/-0.412 | -0.464/-0.368 |
| Galectin3 修复 k1 | 0.201/0.183 | -0.692/-0.486 | -0.543/-0.387 | -0.486/-0.294 |
| FXa E40 repaired | 0.473/0.382 | 0.528/0.394 | **0.529/0.377** | 0.509/0.382 |

这里每格为 Pearson/Kendall。Galectin3 的水/氢修复确实改善了 `tPIE`，但三档 restraint 都没有消除 ES 反向；因此不能再把 Galectin3 失败简单归因于 `k=10` 或“多加了 SOLV”。完整结果位于 `runs/fmo87/ablation/E57_galectin3_structure_quality/`、`E58_galectin3_repaired_k100/` 和 `E59_galectin3_repaired_k1/`。

## 5. 熵的处理

熵不能因为“不收敛”就把当前数值强行当作收敛结果。正确规则是分块延长并在每个检查点重新审计：

```text
trajectory → block distributions → ESS/JS/contact/energy checks
```

只有同时满足以下条件，才允许加入 `-TΔS`：

- bound/free 的身份、脚本哈希、随机种子和时间步一致；
- 轨迹完整、蛋白接触和体系能量没有失稳；
- block split-half JS ≤ 0.15；
- 有效独立样本数达到预先门槛；
- bound/free 的 torsion 定义相同；
- 不是只计算配体扭转 proxy 却宣称完整结合熵。

若不满足，追加下一个预先规定的 MD block；达到最大预算仍不收敛，就报告“熵不可用”，不把最后一个 block 当作收敛熵。当前 E55 的 bound/free JS=`0.412/0.289`，所以熵不进入生产分数。

## 6. 为什么不能保证每个未知 target 都达到 0.55/0.38

当前证据已经显示：

- 原 87 内部 `ES+SOLV` 很好；
- FXa 的 ES 分支仍可用；
- Galectin3 修复后 `tPIE` 变为弱正向，但 ES 分支仍反向；
- 仅凭 pose 间稳定性会在 Galectin3 选出稳定但仍反向的子集；
- 这说明“结构质量通过”与“亲和力排序可迁移”是两个不同性质。

真正普适的工程目标应该是：**不掩盖失败，自动识别不可置信体系，并给出可审计的 FMO 分数/不确定性；在有标签时再进行冻结的 target calibration。** 如果一个未知 target 的所有零拟合 FMO 分支都没有方向证据，流水线应输出 `FMO_UNCERTAIN`，而不是伪造一个看似优秀的排序。

Boltz2 可以作为外部对照或 fallback，但不能与 FMO kcal/mol 直接相加。当前 Galectin3/FXa 的 Boltz2 结果已经完成同例审计，见 `runs/fmo87/boltz_external_analysis/`。

## 7. 生产输出约定

当前跨数据集的机器可读汇总为：
`runs/fmo87/ablation/final_analysis/cross_target_summary.json`（同时有 CSV）。
它明确记录原 87、Galectin3 和 FXa 的 FMO 状态及独立 Boltz2 fallback，避免把
不同模型的分数拼成一个无物理意义的总能量。

每个 target 最终必须保存：

```text
structure_audit.json
repair_manifest.json（若修复）
pose_arm_manifests/
blind_predictions.json（每个 arm）
pose_ensemble_selection.json
entropy_audit.json（若尝试熵）
final_ranking.csv
```

其中 `final_ranking.csv` 必须写明：结构协议哈希、受体/配体哈希、FMO 参数、选用分支、是否加入熵、失败/不确定状态和覆盖率。任何只成功一部分 case 的结果都不能称为完整 target 排名。
