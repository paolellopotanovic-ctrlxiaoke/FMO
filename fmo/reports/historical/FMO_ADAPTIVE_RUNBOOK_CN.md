# FMO 自适应生产 Runbook

日期：2026-10-02。该流程是无标签执行协议；实验亲和力只能在流程冻结后用于独立审计，不能用于选择本次运行的结构、参数或分数。

## 1. 冻结输入

记录 target、蛋白 PDB、配体 SDF、质子化/净电荷、晶水策略、PDB2PQR 版本、口袋半径、GAMESS/FMO 参数和所有输入 SHA256。

## 2. 结构门

```bash
PY=/path/to/fmo87/bin/python
$PY scripts/audit_fmo_external_structure_quality.py ...
$PY scripts/repair_fmo_external_structures.py ...
$PY scripts/audit_fmo_external_structure_quality.py ...
```

必须检查全部 receptor 原子（包括水氢）和配体显式氢。出现硬碰撞时统一删冲突水、重新生成配体氢、保留重原子姿态，然后整组重新最小化和审计。不能只修单个配体后混入同一排序。

## 3. 化学与覆盖门

```bash
$PY scripts/audit_fmo_external_chemistry.py ...
$PY scripts/audit_fmo_qm_coverage.py \
  --qm-root ... --ligands ... --output ...
```

任意 case 缺少成功 `result.json`，target 都不能发布完整 FMO 排名，状态必须是 `FMO_FALLBACK` 或 `FMO_REVIEW`。

### 配体形式电荷与微观状态

- 保留属性剥离后 SDF 中显式氢、键级、形式电荷和三维坐标；不得为了让计算通过而统一中和，也不得按实验亲和力挑选质子化态。
- 结构门验证 SDF 电荷与 cohort 声明一致；FMO deck 必须把配体 fragment 的电荷写入 `ICHARG`，并通过电子数/自旋奇偶检查。外层预检的 `PASS` 只表示结构和输入可计算，不表示混合质子化态的分数已成为 pH 条件下结合自由能。
- 同一数据集若含不同形式电荷，可先完整运行同一冻结 QM 协议做诊断排名，同时报告每个电荷态覆盖；无 assay pH/微观态布居模型时，不能把混合电荷分数升级成普适生产 ΔG。
- 评分前可用去氢/中和后的立体化学图识别“同一共价化学实体”的微观态；按实体聚合预测和评价，不能把一个分子在多个质子化态的记录重复计成多个独立 ligand。标签相等性只在盲预测冻结后作数据完整性核验，不能参与预测态选择。
- 对单碱性位点，若存在可信的微观态 pKa 和 assay pH，可估算游离配体态布居：

\[
f_{BH^+}=\frac{1}{1+10^{\mathrm{pH}-\mathrm{p}K_a}},\qquad f_B=1-f_{BH^+}.
\]

  主态筛选可作 label-free 的结构准备/敏感性分支；它仍忽略结合诱导的 pKa 位移。pKa 缺失、态不全或实验 pH 未核实时，不能声称完成 pH-aware 预测。
- 若未来有可信的溶液态微观态布居，严格的多态结合常数应先在常数空间求和，而非平均 kcal/mol：

\[
K_{\mathrm{app}}=\sum_{r,s}p_r^{\mathrm{R}}p_s^{\mathrm{L}}K_{rs},\qquad
\Delta G_{\mathrm{app}}=-RT\ln(K_{\mathrm{app}}C^\circ).
\]

  这里的 \(p_r^{\mathrm{R}}\)、\(p_s^{\mathrm{L}}\) 是给定 pH 下游离受体/配体微观态布居，\(K_{rs}\) 是对应微观态结合常数。tPIE/PIEDA 不是 \(\Delta G_{rs}\)，因此不能直接代入 Boltzmann 权重；必须另有各微观态自由能模型或实验约束。

## 4. 姿态/几何分支

- 无硬碰撞、位移正常：统一运行 `k=10`。
- 结构修复、系列异质性或几何位移 review：整组运行 `k=1/10/100`。
- 只用跨 arm Kendall、失败率和输入闭合判断数值稳定性；稳定不等于准确。
- 不允许按 ligand 后验切换 restraint 或系列特化归一化。

### 4.1 多构象 ensemble（预注册、全纳入）

多个构象可以作为物理不确定性/构象系综的显式消融，但必须在读取实验标签前冻结
构象生成规则、构象数量、随机种子、受体/配体哈希、权重和聚合公式。每个 ligand 的
所有预注册构象都必须计算并纳入；不能按较低能量、较好排序、较小 RMSD 或成功结果
后验删构象，也不能给不同 ligand 使用不同构象数。

统一入口为 `scripts/aggregate_fmo_conformer_ensemble.py`。默认聚合是固定权重算术
平均：\(\bar E=\sum_i w_iE_i\)，其中 \(\sum_iw_i=1\)。同时报告但不按 ligand
选择的诊断量为 \(E_{ens}=-RT\ln\sum_iw_i\exp[-E_i/(RT)]\)。这里的 FMO
端点分数不能自动解释为严格结合自由能，因此 Boltzmann 量只是预先声明的诊断分支。
缺失构象、混合结构/Hamiltonian fingerprint 或部分 QM 失败时，聚合器直接 fail-closed，
不能用剩余成功构象代替全 ensemble。跨构象稳定性门只表示数值复现性，不表示准确性。

## 5. FMO 分支

所有完整 arm 同时保存：

- `tPIE = ES + EX + CT + DI + SOLV`
- `ES`
- `ES+SOLV`
- `ES+CT+SOLV`

`ES+SOLV`、`ES+CT+SOLV` 是 tPIE 的子集，不能与 tPIE 直接相加。PCM、PBSA、三态差值和 native decomposition 作为独立候选/诊断，除非参考态和溶剂边界完全一致，否则不混成一个物理能量。

## 6. 尺度层

若预先声明要组合多个候选，只能先做 target 内 robust-z 或 rank 聚合，再等权平均：

\[
S=\frac1m\sum_j \frac{S_j-\mathrm{median}(S_j)}{1.4826\,\mathrm{MAD}(S_j)}.
\]

这表示排序尺度聚合，不表示不同模型的 kcal/mol 可相加。当前该分支只作为审计候选，尚未有外部 prospective 验证。

## 7. 熵门

只有 bound/free 身份一致、多个 replica 完整、有效样本数达标、split-half JS 通过且独立 block 重现时，才允许加入熵分支。只有 torsion proxy 或未收敛轨迹时，状态必须为 `ENTROPY_OMIT`，不能用最后一个 block 代替收敛自由能。

当单个长轨迹受到节点 wall-time 或 GPU 稳定性限制时，使用 E61 的 chunked 协议：
每个 replica 由多个短 block 组成，每个 block 使用独立 seed、独立目录和 provenance
hash；block 不是额外 replica，必须先按 replica 拼接，再计算 ESS、split-half JS 和
bound/free 一致性。任何 block 缺失、被外部终止或 hash 不一致，都不能被当作有效样本。
拼接只用于估计合并分布/熵，**不得**把 restart 边界两侧的帧当作连续时间点计算
自相关；lag-1 ESS 要在各 block 内分别估计后再按 replica 汇总。收敛检查比较每个
replica 早/晚 block 组的分布，不能先把 replica 和 block 全混合后再做一次 split-half。

收敛后用 `scripts/score_fmo_entropy_increment.py` 生成增量报告。该入口会再次检查
完整覆盖和收敛标志；未收敛时直接失败。`--allow-torsion-proxy` 只允许导出明确标记
为诊断的 torsion proxy，`ranking_eligible=false`，不能进入生产排序；只有审计报告
明确提供完整 binding entropy 时，才允许生成可用于排序的熵增量。

## 8. 自动判定与交付

```bash
$PY scripts/derive_fmo_adaptive_protocol.py ...
$PY scripts/compile_fmo_factor_decision_table.py ...
$PY scripts/compile_fmo_target_profiles.py ...
$PY scripts/validate_fmo_adaptive_bundle.py \
  --root runs/fmo87 \
  --output runs/fmo87/ablation/final_analysis/adaptive_bundle_validation.json
```

交付至少包含：

- `fmo_target_profiles.json/csv`
- `fmo_method_matrix.json/csv`
- `fmo_target_readiness.json/csv`
- 结构、化学、QM coverage、姿态和熵 audit
- 熵未收敛时的 `entropy_recovery_plan.json`
- 所有失败 case 及其错误原文
- `adaptive_bundle_validation.json` 必须 `passed=true`；warning 必须在交付说明中解释

## 9. 外部标签解盲门

在盲预测、QM 日志和结构/化学审计全部冻结后，才允许读取实验标签。若 benchmark
同时提供带属性 SDF、绝对 `ΔG` 表和相对 `ΔΔG` 边表，先逐 ligand 核验身份、单位和
符号：

```bash
$PY scripts/audit_fmo_external_label_provenance.py \
  --target <Galectin3|FXa> \
  --output runs/fmo87/ablation/final_analysis/<target>_label_provenance_audit.json
```

门槛是 ligand 集一一对应、SDF `r_exp_dg` 与绝对标签表在 `0.002 kcal/mol` 内一致，
且每条边满足 `ddG(1→2)=G2−G1`。该步骤只证明标签身份/口径/方向一致，不证明 FMO
准确性，也不参与分支选择。保留原带标签 source SHA，同时 QM 只能读取已去除全部
properties 的结构副本；盲预测和输入哈希冻结前不得读取 source labels。

## 10. 当前已知结果

- 原 87 例：`ES+SOLV` 是当前冻结的无拟合主分支，Pearson `0.633634`、Kendall `0.507750`。
- Galectin3：结构修复后 26/26 成功，但跨 arm 稳定不等于方向正确，当前 `FMO_REVIEW`。
- FXa：E71 在修复结构上使用显式记录的 E40 参考口袋完成 27/27 QM，当前为
  `FMO_REVIEW`；标准新 union-pocket 仍需独立验证，E71 不能直接冻结为生产协议。
- 熵：E61 仍需完成 bound/free 多 replica 收敛审计；在此之前生产分数中保持 `OMIT`。

## 11. 每个新 target 的统一执行入口

完成结构准备和所有目标内计算臂后，必须生成可执行部署审计：

```bash
PYTHONPATH=src python scripts/build_fmo_adaptive_deployment.py \
  --target <TARGET> \
  --arm k1 <REPORT_K1> --arm k10 <REPORT_K10> --arm k100 <REPORT_K100> \
  --structure <POSTREPAIR_AUDIT> \
  --chemistry <CHEMISTRY_AUDIT> \
  --output <DEPLOYMENT_REPORT>
```

该报告的判定顺序固定为：

1. 结构 `effective_quality_gate` 和全部 QM case 完整性；
2. 同一配体身份及跨臂 `tPIE`/`ES+SOLV` 秩一致性；
3. 独立姿态方向验证；
4. bound/free 熵门（未通过则 `OMIT`）；
5. 输出冻结的 `ES+SOLV` 或 `FMO_REVIEW/FMO_FALLBACK`。

部署报告只有在每个 arm 都显式通过以下完整性检查时才生成完整排序：请求数等于
逐例数、ligand 身份非空且唯一、状态为成功、四个候选分数有限，且各 arm ligand
身份及顺序完全一致。空 arm、缺分数、部分失败或身份错位一律 fail-closed：不输出
聚合排名并返回 `FMO_FALLBACK`。此门仅证明覆盖与数据闭合，不证明预测方向正确。

`k=1/10/100` 的一致性不能单独升级为 `FMO_READY`。最后运行
`scripts/validate_fmo_adaptive_bundle.py`，该检查会拒绝缺少部署报告、未验证姿态却
标成 `FMO_READY`、或不完整 QM 暴露生产分数的 bundle。

## 11. E61 CUDA 运行时恢复

若 OpenMM 8.4 wheel 在当前驱动上报告
`CUDA_ERROR_UNSUPPORTED_PTX_VERSION (222)`，不能把它记录为熵不收敛。当前节点的
驱动为 `550.144.03`，隔离的 `OpenMM 8.2.0 + CUDA 12.4` conda 环境已通过 GPU1
最小 probe。E61 可用以下环境变量续跑；调度器会跳过已完成 block：

```bash
export FMO87_PYTHON=/tmp/fmo-openmm-conda/bin/python
export PYTHONPATH=/tmp/fmo-openmm-conda/lib/python3.12/site-packages:$FMO_ENV/lib/python3.12/site-packages:$PWD/src
export LD_LIBRARY_PATH=/tmp/fmo-openmm-conda/lib:${LD_LIBRARY_PATH:-}
export FMO87_PLATFORM=CUDA
export FMO87_DEVICE=0
export FMO87_CUDA_VISIBLE_DEVICES=1
export FMO87_CUDA_COMPILER=nvrtc
bash scripts/run_e61_chunked_serial.sh
```

这只是运行时修复，不改变力场、轨迹长度、block 划分、随机种子或收敛门。

若完整 3-block 审计只显示 ESS/JS 不足而 provenance 全部通过，可在审计后显式扩展
`block_protocol.json` 的 `blocks` 和 `FMO87_BLOCKS`（例如从 3 扩展到 6），再续跑
block 4–6。扩展必须保留原随机种子公式、状态、replica 和力场；不能把阈值调低来
制造“收敛”。
