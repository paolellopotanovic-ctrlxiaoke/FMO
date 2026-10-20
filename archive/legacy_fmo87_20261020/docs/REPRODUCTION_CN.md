# FMO87 完整复现指南

## 1. 引导

```bash
git clone git@github.com:paolellopotanovic-ctrlxiaoke/FMO.git
cd FMO
./scripts/bootstrap_engines.sh
./run_pipeline.sh verify
./run_pipeline.sh check-inputs
```

以上不需要外部结构下载或系统 GAMESS。`verify` 运行回归、CLI、冻结协议和最终参考产物检查；`check-inputs` 验证 87 个身份和全部源/MSA 哈希。

## 2. 完整主链

```bash
./run_pipeline.sh download-boltz
./run_pipeline.sh all
```

`all` 依次执行：input validation、cohort 校验、蛋白制备、OpenMM restrained minimization、FMO2-DFTB3/PIEDA、Boltz 输入和推理、strict audit、合并、Table 12 指标、bootstrap 与 SophosQM 校准消融。

推荐把 FMO 改成并行：

```bash
./run_pipeline.sh prepare
FMO_JOBS=8 ./run_pipeline.sh fmo-parallel
./run_pipeline.sh boltz-inputs
./run_pipeline.sh download-boltz
./run_pipeline.sh boltz
./run_pipeline.sh audit
./run_pipeline.sh merge
./run_pipeline.sh evaluate
./run_pipeline.sh calibrate
```

所有阶段可重复执行；成功 case 会被复用，失败 case 保留原始错误供审计。

## 3. 无 explicit-water 消融

```bash
./run_pipeline.sh no-water
```

该分支在 PDB2PQR/PROPKA 前删除水，重新构建复合物、重新 OpenMM minimization、重新跑 FMO、重新审计并使用同一 Boltz 预测评估。它不是简单从已有 FMO pocket 删除水片段。

## 4. 直接复验参考结果

```text
reference_results/final/metrics.json
reference_results/final/production_audit.json
reference_results/final/boltz_production_audit.json
reference_results/final/FINAL_REPORT_CN.md
reference_results/final/WATER_ABLATION_CN.md
```

这些文件是已完成 87/87 FMO 与 87/87 Boltz 的冻结汇总，不包含原始超大日志。原始复算输出在 `runs/fmo87/`，默认不提交。

## 5. 运行成本

参考单核 FMO smoke case `P38_2ee` 为 43.70 秒。已完成全量 FMO 的 wall-time 汇总见 `reference_results/final/` 与最终报告；实际集群时间取决于 CPU、并行 shard 数和 Boltz GPU。论文报告的 Table 12 是 cohort 汇总指标，不是单 case 指标。
