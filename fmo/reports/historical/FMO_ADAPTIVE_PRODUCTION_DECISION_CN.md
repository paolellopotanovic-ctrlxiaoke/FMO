# FMO 自适应生产判定协议

日期：2026-10-01。该协议是当前数据证据支持的最强可部署版本。它不承诺
任何未知靶点都能达到某个固定相关性阈值；它承诺在计算前冻结规则、完整记录
失败模式，并且不把不可靠的 FMO 排名伪装成成功。

## 1. 输出不是一个无条件总分

对每个 target 同时生成四个独立 FMO 分支：

\[
S_1=tPIE=ES+EX+CT+DI+SOLV,
\quad S_2=ES,
\quad S_3=ES+SOLV,
\quad S_4=ES+CT+SOLV.
\]

`S3` 和 `S4` 已经是 `S1` 的子集，不能再写成 `S1+S3` 或 `S1+S4`。如果
需要混合，只能在预先冻结的候选公式中写成

\[
S_\alpha=\alpha S_1+(1-\alpha)S_3,
\]

并把它明确标为经验 blend，而不是能量相加；当前 blend 消融没有发现跨
Galectin3、FXa、原 87 都有效的固定 `α`。

## 2. 无标签判定树

### Gate A：结构与化学

1. 检查蛋白链、二硫键、肽键、质子化、净电荷、显式氢、晶水和配体三维坐标。
2. 检查所有 receptor 原子，不只检查蛋白重原子；水氢和配体氢的硬重叠必须修复。
3. 统一删除冲突水、重新生成配体氢并保留重原子姿态；修复后必须检查实际
   `minimized_complex.pdb`，不能只检查原始 SDF。
4. 任意 case 未完成审计，target 输出 `STRUCTURE_INCOMPLETE`，不输出完整排序。
5. 对系列化外部队列运行 `scripts/audit_fmo_external_chemistry.py`，逐例记录
   电荷、重原子数、分子量、供受氢键数、可旋转键和元素组成，并核对 SDF
   ligand 集与 QM 结果、PIEDA pair 闭合。系列差异只触发诊断和整系列复核，
   不允许仅凭分子量/相关性对不同系列使用不同归一化或分数。

### Gate B：姿态与数值稳定性

- 直接通过 Gate A：运行统一的 `k=10` arm。
- 发生修复、重原子 review 或大位移：整组运行 `k=1/10/100`，不按单个
  ligand 后验挑 restraint。
- 记录每个分支的跨 arm Kendall、一致性和失败覆盖率；稳定性只表示数值
  稳定，不表示方向正确。

### Gate C：PIEDA 物理闭合

要求每个 case 满足

\[
|tPIE-(ES+EX+CT+DI+SOLV)|<\delta,
\]

并且 pair 数、SCF、PCM 和输入哈希一致。若

\[
r_{cancel}=\frac{|EX|+|CT|+|DI|+|SOLV|}{|ES|+\epsilon}
\]

很大，只能说明附加项会强烈改变排序，不能据此自动决定保留或删除它们。

### Gate D：熵

只有同时满足以下条件才允许产生 `S−TΔS` 分支：

- bound/free 身份、原子顺序、脚本和参数哈希一致；
- 每个 replica 完整且有效独立样本数达标；
- split-half JS 通过预设阈值；
- 计算的是与声明相符的自由能项，而不是只计算 torsion proxy；
- 达标结果能在独立 block/replica 重现。

否则追加预先声明的 sampling block；预算耗尽仍不通过就标记 `ENTROPY_OMIT`，
不把最后一个未收敛 block 加入排序。

## 3. 最终状态

| 状态 | 条件 | 输出 |
|---|---|---|
| `FMO_READY` | 结构、闭合、覆盖和姿态门通过，至少一个分支有冻结的适用规则 | FMO ranking + 全部分支审计 |
| `FMO_REVIEW` | 结构已修复但姿态/分支稳定性不足 | 多 arm ranking + 不确定性，不选隐藏最优分支 |
| `FMO_UNCERTAIN` | 结构可算但所有分支缺乏方向/迁移证据 | 保留 FMO 分数，禁止宣称可靠亲和力排序 |
| `FMO_FALLBACK` | FMO_UNCERTAIN 且允许外部模型 | 独立 Boltz2/实验校准结果，不能与 kcal/mol FMO 直接相加 |

当前数据的正确归类是：原 87 的 `ES+SOLV` 属于 `FMO_READY`；FXa 的 E71
参考口袋诊断臂已完成 27/27 QM，但口袋策略尚未独立验证，因此为 `FMO_REVIEW`；
Galectin3 的 E68 一致性 arm
结构和 QM 完整，但方向仍未验证，因此属于 `FMO_REVIEW`，Boltz2 只能作为独立 fallback；FXa 的标准新口袋失败历史仍需保留为诊断证据。这个判定比为 Galectin3
后验选择 `EX+DI+SOLV` 更可迁移、更可审计。

Galectin3 的 26 例另有 label-free 化学/QM 审计：18/18 的 17 系列和 8/8 的
18 系列均为中性；QM 均成功且 ligand–receptor PIEDA pair 数闭合。可是 18
系列平均 35.25 个重原子、518.50 Da，17 系列平均 29.61 个重原子、435.86 Da，
且 18 系列可旋转键为 8–10、17 系列均为 8。这证明两组并非同一尺度的微小
取代变化，但**不能**证明大分子/熵或某个姿态误差就是反向的原因。18 系列
原始分数失败、系列内后验归一化有改善而全队列规则仍失败；因此目前将其列为
“系列效应待独立验证”，不据此部署 heavy-atom normalization。
逐例报告：`runs/fmo87/ablation/final_analysis/galectin_chemistry_audit.json`；
系列内标签诊断：`runs/fmo87/ablation/final_analysis/galectin_series_diagnostics.json`，
后者明确仅为 post-hoc，不可作生产选择器。

当前所有已完成评分臂的统一矩阵在
`runs/fmo87/ablation/final_analysis/fmo_method_matrix.json`（同时有 CSV）。矩阵
逐行保留 target、结构臂、四个独立分数、Pearson/Kendall、结构状态和熵状态，
因此不能用某一个已解盲的最好单元格替代完整方法验证。

无标签姿态选择器的回顾性压力测试也必须保留：Galectin3/FXa 的跨 restraint
秩稳定性门均通过，并选择了 `ES+CT+SOLV`，但 Galectin3 该分支的 Pearson/
Kendall 仍为 `-0.413/-0.307`。这说明“稳定”只表示数值可重复，不表示方向正确；
生产规则因此在方向规则未被独立新靶点验证前，将选择结果标为
`UNVERIFIED_LABEL_FREE`，并自动降级为 `FMO_REVIEW`，而不是把稳定分支当作准确
预测器。对应审计：
`runs/fmo87/ablation/final_analysis/galectin_adaptive_selection_audit.json` 和
`fxa_adaptive_selection_audit.json`。

## 4. 为什么不能强行保证所有 target 都好

### 4.1 新增的尺度消融结论

直接把多个 FMO/PBSA/PCM 数值相加隐含了“每个项的数值尺度和零点可直接比较”。
在原 87 例上，对预先固定的分支做 target 内 robust-z 标准化

\[
 z_j=(S_j-\operatorname{median}(S_j))/(1.4826\,\operatorname{MAD}(S_j)),
 \qquad S_{agg}=\frac{1}{m}\sum_j z_j
\]

是一个合法的**排序尺度聚合**，但不是能量相加。`tPIE+ES+SOLV` 的 robust-z
回放 Kendall 约为 `0.391`，四分支等权 robust-z 约为 `0.399`，高于原始单独
tPIE 的约 `0.320`；然而 Galectin3 的对应 Kendall 仍约为 `-0.065` 和 `-0.195`。
因此尺度处理确实减少了一个误差来源，但不能把反向的结构/状态问题修复成正确
排序。该消融保存在 `runs/fmo87/ablation/final_analysis/scale_aggregation_audit.json`，
实现为 `scripts/analyze_fmo_scale_aggregation.py`，当前只标记为
`AUDIT_ONLY_UNTIL_PROSPECTIVE_VALIDATION`，不以后验结果选择权重。

在当前证据中，Galectin3 的 ES、ES+SOLV 和 ES+CT+SOLV 在三种 restraint 下
都反向，而 tPIE 仅弱正向；原 87 和 FXa 的相同分支却为正向。由此可证明：

1. 不是统一加权就能修复，因为 `S_α` 扫描没有跨靶点有效的固定权重；
2. 不是简单的坐标/水氢错误，因为修复后仍反向；
3. 不能用无标签稳定性选择器证明准确性；
4. 要真正消除该失败模式，需要新增可验证信息：更好的 binding pose/rotamer
   集合、完整 bound/free alchemical sampling，或独立训练过的外部模型。

因此“普适”在工程上应定义为**普适的输入处理、质量门、分支审计和失败策略**，
而不是未经证据支持地承诺某一个 FMO 能量表达式对所有靶点同样准确。
