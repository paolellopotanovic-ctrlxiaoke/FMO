# CDK8 QM/MM 与 OpenMM 同帧审计

原始单帧审计：2026-10-26；多帧更正：2026-10-27

> **重要更正：** 本文早期写出的 QM 替换修正 `−0.225 kcal/mol` 使用了不匹配的 MM 1–4 缩放参数，现已作废。匹配参数在 `iteration=2000` 的固定帧给出 `−5.823 kcal/mol`。首轮多帧结果还包含 GAMESS 的 QM-center 人工罚能；同坐标重跑后，已将该算法项从正式结果中剥离。请以 [CDK8_QMMM_MULTIFRAME_20261027.zh-CN.md](CDK8_QMMM_MULTIFRAME_20261027.zh-CN.md) 和 [CDK8_QMMM_CENTER_CONTROL_20261027.zh-CN.md](CDK8_QMMM_CENTER_CONTROL_20261027.zh-CN.md) 为准。

## 本轮实际完成

从 CDK8 的 OpenFE 生产轨迹固定取 `iteration=2000` 的物理端态帧，分别构造并运行：

- bound A / bound B；
- solution A / solution B；
- QuanPol DFTB3/AMBER QM/MM 单点；
- 原始 OpenMM XML 端点势能对照。

四个 QuanPol 作业均正常结束，SCF 未报告不收敛，总势能闭合误差为零。

## 数值

| 方法 | A 端态绑定差值 | B 端态绑定差值 | B-A 双差值 |
|---|---:|---:|---:|
| QuanPol 有限壳层 DFTB3/AMBER | −2023.383 | −1875.490 | **147.893** |
| 原始 OpenMM 完整周期体系 | −320786.393 | −320683.999 | **102.394** |

两种双差值相差 `45.499 kcal/mol`。

随后以匹配的 `LJSIGMA=1`、`WT14LJ=0.5`、`WT14CH=0.833333...` 参数运行 MM-only。MM-only 有限壳层双差值为 `153.716 kcal/mol`，而 QM/MM 为 `147.893 kcal/mol`，因此本帧的 QM 替换修正为：

```text
147.893 - 153.716 = -5.823 kcal/mol
```

原来的 `−0.225` 是 1–4 缩放不一致导致的错误比较，不应引用。即使更正后的单帧 `−5.823`，也只代表固定坐标下的有限壳层势能双差值，不能代表自由能修正。

2026-10-27 对 3 条独立重复轨迹、5 个预先固定时点完成了 `60` 组成对 QM/MM–MM 单点，得到 `15` 个四态双差值。`0.030 kcal/mol` 是 **QuanPol MM-only 与导出的 OpenMM FFDATA MM** 的最大绝对势能差，不是 QM/MM 与 MM 的差，也不证明有限壳层 Hamiltonian 等价于 OpenMM PME。

首轮 `SWRB=1001 Å` 结果均值 `−1.598`、SD `8.244 kcal/mol`，含一个 `+25.764` 的极端值。GAMESS 源码审计发现：当 `SWRB² < 1e8 bohr²` 时会启用 QM-center 键型人工罚能；对同帧将 `SWRB` 改为 `6000 Å` 重跑后，60 个 QM 单点的罚能均为零。对照检查确认坐标、环境、MM deck 与除 `SWRB` 外的 QM deck 一致；坐标盒最大直径上界 `122.347 Å`，小于两种 `SWRB`，故没有改变有限壳层的距离截断。

罚能关闭后，双差值均值 `−3.203`、SD `3.595`、范围 `−8.325` 至 `+4.532 kcal/mol`。原极端帧从 `+25.764` 降至 `+2.499`；另外两帧也分别被中心项改变 `0.266` 与 `1.075 kcal/mol`。新旧变化与日志中 QM-center 项的四态双差闭合误差小于 `2.1×10⁻¹² kcal/mol`，MM 能量完全未变。首轮“极端值 +25.8、SD 8.24”的解释因此作废；关闭中心项后仍有固定帧能量波动，但 5 帧/重复不足以估计重加权自由能。

## 解释边界

这不是“QM 比 FEP 差 45.5”的结论。两行不是同一个 Hamiltonian：

1. OpenMM 使用原始周期 PME、完整水盒和 OpenFE 的混合拓扑；
2. QuanPol 目前使用有限局部壳层、非周期截断和 AMBER/QuanPol 力项；
3. QuanPol 导出尚未完全表达 OpenMM 的约束动力学与自定义水相互作用；
4. `FORMS NO BOND` 是原子级拓扑/约束告警，包含零力常数约束键与离子，不能一概解释为“水丢失”；
5. bound/solution 的绝对势能不能直接混合成实验结合自由能。

因此当前 `147.893 kcal/mol` 只能作为“QM/MM 端态能量接口已跑通”的诊断量，不能进入 tPIE 排名，也不能声称接近 FEP。

## 仍未满足的自由能门槛

本轮已通过“同一 FFDATA 导出模型”的 MM-only 能量闭合，但这不等价于 OpenFE 周期 PME Hamiltonian。若要形成可信的 QM/MM 自由能修正，仍须：

```text
先验证周期边界与有限壳层的环境一致性，再对每个端态采样 `ΔU=U_QM/MM−U_MM`，做指数平均/重加权并报告有效样本量
```

建议通过标准：

- 当前重建 FFDATA 模型内 MM-only 对照的最大误差为 `0.030 kcal/mol`；
- OpenFE PME 与 QuanPol 有限壳层的 Hamiltonian 等价性尚未建立；
- 约束键、离子与 QuanPol `FORMS NO BOND` 告警仍需逐类核实；
- 3 个时间点/重复不足以估计指数平均、有效样本量或自由能不确定度。

所以 QM/MM 仍是独立诊断支线；主流程继续使用 `dry + k10 + tPIE`，不把本次双差值并入排序。

## 生成文件

- `tools/build_quanpol_endpoint_pilot.py`
- `tools/export_quanpol_endpoint_forcefield.py`
- `tools/compare_same_frame_openmm_endpoints.py`
- `tools/report_quanpol_endpoint_pilot.py`
- `docs/CDK8_QUANPOL_QMMM_ENDPOINT_PILOT_20261026.json`
- `docs/CDK8_SAME_FRAME_OPENMM_ENDPOINT_AUDIT_20261026.json`
- `docs/CDK8_QMMM_MM_CORRECTION_20261026.json`
- `tools/run_quanpol_multiframe.py`
- `tools/report_quanpol_multiframe.py`
- `docs/CDK8_QMMM_MULTIFRAME_20261027.zh-CN.md`
- `docs/CDK8_QMMM_MULTIFRAME_20261027.json`
- `docs/CDK8_QMMM_MULTIFRAME_WITH_CENTER_20261027.zh-CN.md`（保留首轮、现已更正的结果）
- `docs/CDK8_QMMM_CENTER_CONTROL_20261027.zh-CN.md`
