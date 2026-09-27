# 数据流与 87-case 闭环说明

## 输入层

唯一生产队列为 `data/fmo87_cohort.csv`：

| 字段 | 含义 | 后续消费 |
|---|---|---|
| `target` | CDK2/JNK1/P38/TYK2 | 分靶点制备、合并、指标加权 |
| `ligand_name` | SDF record 名称 | SDF 选_record、目录命名、身份合并 |
| `ligand_sdf` | 冻结来源 SDF provenance | 运行时按 `STRUCTURE_ROOT` 重映射 |
| `protein_pdb` | 冻结来源 PDB provenance | Boltz 蛋白序列读取，运行时重映射 |
| `protein_chain` | 目标链 | Boltz sequence |
| `canonical_isomeric_smiles` | RDKit canonical identity | 去重、输入体检、Boltz ligand |
| `formal_charge` | 中性检查 | cohort/audit |
| `experimental_dg_kcal_mol` | 实验 ΔG | 仅评估层使用，不参与 pose 选择 |
| `source_record_count` | 折叠前重复数 | JNK1 alternate-flip 审计 |
| `representative_rule` | 代表 pose 规则 | 队列 provenance |

`ligand_sdf` 和 `protein_pdb` 中的旧 `vendor/...` 字符串是原始构建路径记录，不要求用户复刻。`validate_inputs.py`、`run_fmo87.py` 和 `build_boltz_fmo87_inputs.py` 均按文件名映射到 `STRUCTURE_ROOT`。

## 处理链

```text
8 个源 SDF/PDB + 87 行冻结 cohort
  ↓ check-inputs：哈希、身份、电荷、配体名、标签
  ↓ prepare：PDB2PQR/AMBER/PROPKA pH 7.4、补氢/OXT/H-bond
  ↓ fmo / fmo-parallel：OpenMM AMBER14/OpenFF restrained minimization
  ↓ fragmentation：5 Å pocket、shifted Cα–C′、ACE/NME、HOP_C、水策略
  ↓ GAMESS FMO2-DFTB3/3OB-3-1/PIEDA
  ↓ parser：真实 PIEDA pair table -> ligand tPIE
  ↓ audit：normal termination、SCF/PIEDA、TIE reconciliation、协议哈希
  ↓ Boltz branch：SMILES+protein sequence+shared MSA -> official Boltz-2
  ↓ merge：target+ligand identity 精确合并
  ↓ evaluate：逐靶点 Pearson/Kendall，再按化合物数加权
```

## 质量门

- 输入体检必须为 `PASS`；
- FMO 与 Boltz-2 均必须 87/87；
- FMO 必须通过 selected final TIE 与本地 pair-sum tPIE 的 ≤0.005 kcal/mol 对账；
- 不允许 GAMESS return code 0 但 SCF/PIEDA/geometry/electron failure 的结果进入评估；
- 合并表不能有缺失、重复或标签不匹配；
- 显式水分支必须声明 `retained`、`removed_from_fmo_pocket` 或 `removed_before_preparation`。

## 本地参考结果

`reference_results/final/` 保存已完成 87/87 双分支输出，包括 FMO/Boltz 特征、预测、审计、指标、bootstrap、校准消融和水消融。该目录用于审查和绘图，不依赖原始 GAMESS/Boltz 日志；从零重建仍按 `docs/REPRODUCTION_CN.md` 执行。
