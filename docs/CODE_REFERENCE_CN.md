# 当前代码与工具索引

仓库只保留 FMO87 主生产链。所有 Python 文件都直接服务于队列构建、双分支计算、审计、合并、评估或消融。

## 一条流水线

```text
1  build_fmo87_cohort.py        冻结 87 个化合物身份
2  prepare_protein.py           四靶点蛋白制备
3  run_fmo87.py                 几何优化 + FMO 分片 + GAMESS
4  build_boltz_fmo87_inputs.py  Boltz-2 输入与 MSA manifest
5  run_boltz_fmo87.py           Boltz-2 affinity 推理
6  audit_fmo87_production.py    FMO strict audit
7  audit_boltz87_production.py  Boltz-2 strict audit
8  merge_fmo87_outputs.py       双分支精确合并
9  evaluate_fmo87.py            Table 12 指标 + bootstrap
10 ablate_sophosqm_calibration.py  校准来源消融
11 validate_inputs.py            87-case 输入体检
12 collect_fmo_shards.py         并行 shard 收集
13 summarize_water_ablation.py   显式水消融汇总
```

完整命令按 `docs/REPRODUCTION_CN.md` 执行。

## CLI 入口 `scripts/`

| 顺序 | 文件 | 输入 | 输出 | 说明 |
|---:|---|---|---|---|
| 1 | `build_fmo87_cohort.py` | 公开 benchmark structure root | `data/fmo87_cohort.csv` | canonical SMILES 去重，冻结 16/21/34/16 计数。 |
| 2 | `prepare_protein.py` | 单靶点蛋白 PDB | prepared PDB/PQR/log/provenance | PDB2PQR 3.6.1、AMBER、PROPKA pH 7.4。 |
| 3 | `run_fmo87.py` | cohort、prepared protein、GAMESS | minimized complex、fragment/deck/result | OpenMM 几何、ACE/NME/HOP 分片、GAMESS FMO2-DFTB3/PIEDA。 |
| 4 | `build_boltz_fmo87_inputs.py` | cohort、target MSA | Boltz YAML 与 manifest | 共享 MSA、canonical SMILES、输入哈希。 |
| 5 | `run_boltz_fmo87.py` | input manifest、Boltz cache | affinity JSON/CSV/provenance | Boltz 2.2.1、seed 0、200/5 affinity 设置。 |
| 6 | `audit_fmo87_production.py` | FMO target roots、prepared proteins | production audit | normal termination、SCF、化学协议、TIE、哈希、explicit-water policy 全检查。 |
| 7 | `audit_boltz87_production.py` | Boltz target roots/cache | Boltz audit | 输入、MSA、checkpoint、软件版本、JSON/CSV 一致性。 |
| 8 | `merge_fmo87_outputs.py` | 已审计双分支输出 | merged feature/prediction tables | 按靶点/配体身份精确合并。 |
| 9 | `evaluate_fmo87.py` | merged outputs | metrics/models/REPORT | 逐靶点指标按化合物数加权，10000 次 bootstrap。 |
| 10 | `ablate_sophosqm_calibration.py` | final compounds、官方 SI | calibration ablation | 区分 tPIE、clogP、拟合系数和 FMO 协议误差。 |
| 11 | `validate_inputs.py` | cohort、structure root | input validation JSON | 校验 cohort 哈希、87 身份、SDF/PDB 源哈希、canonical SMILES 和电荷。 |
| 12 | `collect_fmo_shards.py` | parallel shard roots | auditable target roots | 去重、检查 87/87，并收集 minimized/FMO 输出到靶点目录。 |
| 13 | `summarize_water_ablation.py` | retained/scorer/preparation 三组 models 与 compounds | water ablation CSV/JSON/Markdown | 输出全队列加权指标、P38 成对 tPIE 位移和中文/英文报告。 |

## 核心库 `src/sophosqm_baseline/`

| 文件 | 职责 |
|---|---|
| `cohort.py` | 读取四靶点 SDF/PDB，生成 canonical 身份、代表记录，并将冻结来源路径映射到本地 structure root。 |
| `structure.py` | PDB/SDF 坐标、残基、配体和元素/电荷表示。 |
| `preparation.py` | 配体过滤、结晶水保留或制备前去水、PDB2PQR/PROPKA 调用与 provenance。 |
| `minimization.py` | OpenMM AMBER14/OpenFF 灵活复合物 restrained minimization。 |
| `fragmentation.py` | 5 Å pocket、shifted Cα–C′ fragment、ACE/NME、explicit-water 保留/排除、HOP 与电子校验。 |
| `fmo_input.py` | 生成 FMO2-DFTB3/PIEDA GAMESS 输入，包括 charge/multiplicity/FMOBND。 |
| `gamess.py` | GAMESS 安装检查与哈希/manifest 辅助。 |
| `jobs.py` | 组装蛋白–配体 job 与 fragment manifest。 |
| `runner.py` | 隔离环境执行 GAMESS，实施超时、SCF/PIEDA/TIE 质量门。 |
| `parser.py` | 解析真实 GAMESS PIEDA pair table 并聚合 ligand tPIE。 |
| `boltz_inputs.py` | 生成 Boltz YAML、共享 MSA manifest 和输入哈希。 |
| `exact_evaluation.py` | 身份合并、指标、bootstrap、SophosQM-style 校准层。 |
| `__init__.py` | package 声明。 |

## 测试 `tests/`

| 文件 | 覆盖 |
|---|---|
| `test_cohort.py` | 87-case 身份、JNK1 重复折叠、标签无关代表选择。 |
| `test_gamess.py` | 化学分片、ACE/NME/OXT、HOP、闭壳层、deck 与 runner。 |
| `test_parser.py` | GAMESS PIEDA 解析、tPIE 聚合、selected TIE。 |
| `test_exact_evaluation.py` | strict audit、指标、bootstrap、模型层。 |
| `test_preparation.py` | 蛋白输入过滤、制备级显式水移除、来源路径重映射、策略不匹配时禁止复用旧 GAMESS 成功结果、并行 shard 收集。 |

运行：

```bash
python -m pip install -e .[test]
pytest
```

## 代码边界

`src/` 和 `scripts/` 是本项目自研 orchestration 与分析代码。GAMESS-US 与公开 benchmark 输入已内置；OpenMM/OpenFF/PDB2PQR/Boltz 包由 bootstrap 安装；Boltz 官方权重按哈希下载。
