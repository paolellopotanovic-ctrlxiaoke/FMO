# 项目总览

## 研究目标

本项目实现一套独立、可审计的蛋白–配体 FMO scoring 工具链，并在 Boltz-2 论文 Table 12 使用的四靶点 87 个中性化合物子集上，与论文报告的 Boltz-2 和内部 FMO baseline 比较。核心检验不是复刻私有软件，而是回答：

1. 用公开结构、公开制备工具和 GAMESS 官方 FMO2-DFTB3/PIEDA，能否稳定完成全部真实蛋白–配体 case；
2. 自研 raw tPIE 与论文内部 FMO 的排序指标差异有多大；
3. 差异主要来自 SophosQM 拟合参数，还是 pose/几何/溶剂/私有 scorer 等 FMO 协议差异。

## 冻结队列

`data/fmo87_cohort.csv` 是唯一主队列，SHA-256：

```text
19edfb8e2ccf3ab1029d848e4062d70cb35532d2b4fee63674f160b2a292997b
```

| 靶点 | 化合物数 |
|---|---:|
| CDK2 | 16 |
| JNK1 | 21 |
| P38 | 34 |
| TYK2 | 16 |
| 合计 | 87 |

队列由公开 Ross/Chen structure inputs 重建，以 RDKit canonical isomeric SMILES 去重。原 104 条记录中的 17 条 JNK1 alternate-flip 重复结构合并，代表 pose 选择不使用实验标签。

## 最终结果

| 模型 | 样本加权 Pearson | 样本加权 Kendall | 论文参考 |
|---|---:|---:|---:|
| Local Boltz-2 | 0.651763 | 0.509589 | Boltz-2 0.66 / 0.48 |
| 自研 raw FMO2-DFTB3 tPIE | 0.479233 | 0.320055 | 内部 FMO 0.55 / 0.38 |

Local Boltz-2 与论文值几乎一致，说明 87 化合物队列、单位转换和评估器可靠。自研 FMO 点估计较低，但 10000 次靶点分层 bootstrap 的 95% CI 覆盖论文 FMO 参考值，因此不能声称统计显著劣化。

次要 SophosQM-style 校准显示：

| 模型 | Pearson | Kendall | 解释 |
|---|---:|---:|---|
| raw tPIE | 0.479233 | 0.320055 | 主比较，未拟合 |
| tPIE+clogP in-sample | 0.608842 | 0.499633 | 乐观训练结果 |
| tPIE+clogP LOO | 0.401952 | 0.330867 | 外推稳定性不足 |

官方 SophosQM SI 的 p38 系数迁移到本地 p38 后，Pearson 只从 0.568153 到 0.589925。这说明主要差距不是简单缺少原文拟合系数。

## 方法路线

```text
public FEP benchmark structures
    ↓ frozen 87-compound cohort
    ├── FMO branch
    │   ├─ PDB2PQR/PROPKA protein preparation
    │   ├─ OpenMM AMBER14/OpenFF restrained geometry
    │   ├─ shifted Cα–C′ fragments + ACE/NME + HOP_C
    │   ├─ GAMESS FMO2-DFTB3/3OB-3-1/PIEDA
    │   └─ strict parsing/TIE reconciliation audit
    └── Boltz-2 branch
        ├─ shared target MSA inputs
        ├─ Boltz 2.2.1 affinity predictions
        └─ input/checkpoint/JSON provenance audit
    ↓ exact identity merge
    ↓ per-target, count-weighted Table 12 metrics
    ↓ bootstrap and calibration ablation
```

## 必须遵守的解释边界

- 本项目不是 Boltz 作者私有 FMO scorer 的复现。
- 不使用或复刻 Schrödinger Glide、MacroModel、OPLS3e 或 VolSurf。
- 主 FMO 比较使用 unfitted raw tPIE；SophosQM-style 回归是次要分析。
- raw tPIE 是 enthalpic ranking score，不是绝对结合自由能。
- Ross/Chen reference poses 和 OpenMM/OpenFF 几何是明确声明的替代协议。
- 自研 DFTB3 协议使用 water PCM 稳定极性 dimer；这是与私有实现不同的实现选择。

## 当前完成度

87/87 FMO 与 87/87 Boltz-2 均成功。化学制备、几何、fragment/cap/HOP、闭壳层电子计数、SCF/PIEDA、TIE 对账、输入哈希、checkpoint provenance、样本加权指标和 bootstrap 均通过严格审计。完整通过条件见 `COMPLETION_AUDIT_CN.md`。
