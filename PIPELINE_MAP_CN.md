# 完整流水线地图

| 阶段 | 入口 | 引擎/数据 | 主要输出 | 成功门 |
|---|---|---|---|---|
| 环境引导 | `scripts/bootstrap_engines.sh` | pip 锁定版本 + bundled Boltz commit | `envs/fmo87/` | OpenMM/RDKit/Boltz 可导入 |
| 输入体检 | `./run_pipeline.sh check-inputs` | bundled 87 case + MSA | `runs/fmo87/input_validation.json` | 87 身份、SMILES、电荷、源哈希、MSA 哈希全部匹配 |
| 蛋白制备 | `./run_pipeline.sh prepare` | PDB2PQR/PROPKA | prepared PDB/PQR/log/provenance | PDB2PQR 3.6.1、pH 7.4、氢键优化 |
| 标准几何 + FMO | `./run_pipeline.sh fmo` 或 `fmo-parallel` | OpenMM + bundled GAMESS | minimized complex、deck、log、result | 正常终止、无 SCF 失败、PIEDA 合理、TIE 对账 |
| Boltz 输入 | `./run_pipeline.sh boltz-inputs` | bundled cohort/MSA | YAML + manifest | 87 个身份闭环 |
| Boltz 权重 | `./run_pipeline.sh download-boltz` | Boltz 官方 URL | `engines/boltz-cache/` | conf/aff/CCD 存在 |
| Boltz 推理 | `./run_pipeline.sh boltz` | Boltz2 2.2.1 | affinity JSON + provenance | checkpoint/软件/输入哈希一致 |
| 审计合并 | `./run_pipeline.sh audit`、`merge` | 自研审计器 | strict audits、features/predictions | 87/87，无替代或跳过 |
| 指标 | `./run_pipeline.sh evaluate` | Table 12 口径 | metrics/bootstrap/报告 | 逐靶点拟合后按化合物数加权 |
| 校准消融 | `./run_pipeline.sh calibrate` | SophosQM SI 系数 | calibration ablation | 原文系数与本地拟合分开 |
| 去水消融 | `./run_pipeline.sh no-water` | PDB2PQR + OpenMM + GAMESS | 完整无 explicit-water 分支 | 重新制备、重新 relax、重新审计 |

## 引擎边界

- GAMESS 二进制和参数已内置；运行器自动挂载本地 MKL/OpenMPI 和 3OB-3-1。
- Python 化学栈由锁定文件安装，不复制 14 GB 的本机 conda 环境。
- Boltz2 官方源码已按 commit 内置；多 GB 权重不能作为普通 Git blob 提交，`download-boltz` 使用官方 fallback URL 并保留哈希审计。
- 原始运行产物、Boltz cache 和虚拟环境不进入 Git；最终可审计汇总保存在 `reference_results/final/`。
