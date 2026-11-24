# FMO87 完成审计（2026-09-26）

## 通过条件

| 审计项 | 结果 |
|---|---:|
| Frozen cohort SHA-256 | 通过 |
| 四靶点计数 16/21/34/16 | 通过 |
| 87 个配体全部中性 | 通过 |
| PDB2PQR/PROPKA pH 7.4 日志 | 通过 |
| 87 个 flexible restrained geometry | 通过 |
| 100 kJ/mol/Å² ligand restraint | 通过 |
| AMBER14 + OpenFF 2.2.0 | 通过 |
| full ACE/NME shifted fragments | 通过 |
| C 末端 NME 取代 OXT | 通过 |
| `$FMOBND HOP_C` | 通过 |
| 全片段闭壳层 | 通过 |
| DFTB3/3OB-3-1 + `DAMPXH` | 通过 |
| PIEDA `MODORB=3` | 通过 |
| 全 dimer `NACUT=0`/`MODMOL=1` | 通过 |
| water PCM `IFMO=-1` | 通过 |
| GAMESS normal termination | 87/87 |
| 无 SCF/distance/electron 错误 | 87/87 |
| tPIE 与 selected final TIE 对账 | 87/87 |
| deck 哈希匹配 | 87/87 |
| Boltz2 输入哈希与 canonical SMILES | 87/87 |
| Boltz2 checkpoint/软件/参数 provenance | 87/87 |
| Boltz2 JSON/CSV 数值一致 | 87/87 |

完整机器可读审计：

- `runs/fmo87/final/production_audit.json`
- `runs/fmo87/final/boltz_production_audit.json`

## 生产统计

- FMO：87/87 成功；无 SCF 失败。
- GAMESS wall time：mean 39.62 s，median 36.59 s，max 61.46 s。
- 片段总数：3089；HOP 边界总数：1755。
- Ligand RMSD after minimization：mean 0.0663 Å，max 0.3530 Å。
- 单原子最大位移：max 1.7511 Å。
- 唯一真实 C 末端 OXT 进入 pocket 的样本是 `CDK2_1oiu`；NME 正确取代其 LEU296 OXT。

## 结论

主生产队列是完整 87 化合物 official-like 自研工具链。所有必须通过的化学、拓扑、执行、解析、身份和 provenance 门均为通过。
