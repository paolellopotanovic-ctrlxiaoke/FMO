# FMO87 / Boltz-2 私有完整流水线

本仓库是内部研究用的 pull-and-run 工具包，面向 Boltz-2 论文 Table 12 的四靶点 FMO 子集：CDK2/JNK1/P38/TYK2，共 87 个中性配体。

## 已内置

- 87 case 的 4 个蛋白 PDB、4 个多配体 SDF 和源哈希；
- 4 个共享 Boltz MSA 及冻结哈希；
- GAMESS-US 2024.2.1 源码、可执行文件、3OB-3-1 参数、MKL/OpenMPI 运行库；
- 严格生产协议、审计、指标评估、显式水消融和最终参考表；
- SophosQM SI 与 Boltz-2 论文 PDF（仅限本私有内部仓库使用）。

Boltz-2 Python 包由 bootstrap 安装。官方 Boltz2 权重单个超过 GitHub 100 MB 文件上限，因此用 `download-boltz` 从官方地址下载到 `engines/boltz-cache/`。

## 一键启动

```bash
./scripts/bootstrap_engines.sh
./run_pipeline.sh verify
./run_pipeline.sh check-inputs
./run_pipeline.sh download-boltz
./run_pipeline.sh all
```

首次完整运行需要 GPU 运行 Boltz-2；FMO 可用 `FMO_JOBS=N ./run_pipeline.sh fmo-parallel` 并行。默认路径已经是 bundled 数据和 bundled GAMESS，不需要再设置 `STRUCTURE_ROOT`、`GAMESS_DIR` 或 `GAMESS_LIB`。

## 入口

- 总控：`run_pipeline.sh`
- 流程地图：`PIPELINE_MAP_CN.md`
- 引擎说明：`docs/ENGINE_SETUP_CN.md`
- 复现指南：`docs/REPRODUCTION_CN.md`
- 最终分析：`docs/FINAL_REPORT_CN.md`
- 代码索引：`docs/CODE_REFERENCE_CN.md`

## 当前结论

本地 Boltz-2 为 Pearson/Kendall `0.6518/0.5096`，接近论文 `0.66/0.48`，用于验证 cohort 和指标口径。自研 raw FMO2-DFTB3 tPIE 为 `0.4792/0.3201`；论文内部 FMO 为 `0.55/0.38`。后者协议未公开，不能称为同一实现的失败复现。
