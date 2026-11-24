# FMO 协议敏感性实验

## 目的

本实验区分三类因素：GAMESS DFTB dispersion、PCM、以及 FMO DFTB 的 `MODESD/MODGAM` 分支。所有协议复用相同的 87-case 队列、公开结构输入、已生成的标准最小化几何和 fragment 规则；变化只写入 manifest。

## 协议

- `current_pcm`：生产 baseline，PCM + `DISP=NONE` + 默认 GAMESS DFTB 模式。
- `pcm_uff_default`：PCM + `DISP=UFF`，保留默认 `MODESD/MODGAM`。
- `gas_uff_default`：无 PCM + `DISP=UFF`，保留默认模式。
- `gas_none_default`：无 PCM + 无 dispersion，保留默认模式。
- `candidate_a_gas_uff` / `candidate_b_pcm_uff` / `gas_none`：显式 `MODESD=2, MODGAM=1` 的原始 FMO-DFTB 风格候选；当前 GAMESS build 中不可稳定运行。

## 12-case结果

结果文件：

- `runs/fmo87/sensitivity_12/summary.json`
- `runs/fmo87/sensitivity_12_default_modes/summary.json`

默认模式对照中，`pcm_uff_default` 成功 10/12，`gas_uff_default` 成功 3/12，`gas_none_default` 成功 3/12。gas 变体主要是 SCF 不收敛；显式 `MODESD=2` 变体还触发 GAMESS `lumpesd_` segmentation fault。故不能把 gas 或 MODESD 变体的部分成功相关性当作可靠指标。

## 87-case 全量结果

`runs/fmo87/pcm_uff_default_full/summary.json` 与 `evaluation/summary.json`：

| 协议 | 成功数 | Pearson | Kendall | 结论 |
|---|---:|---:|---:|---|
| current PCM baseline | 87/87 | 0.4792 | 0.3201 | 冻结生产协议 |
| PCM + UFF 默认模式 | 79/87 | 0.2068 | 0.3324 | TYK2 提升，但 CDK2/P38 恶化，不采用 |

候选 per-target 成功子集为：CDK2 `-0.6559/-0.2821`（13/16）、JNK1 `0.6701/0.5107`（21/21）、P38 `-0.1185/0.2299`（29/34）、TYK2 `0.8892/0.7833`（16/16）。由于候选不完整且整体指标更差，不进行论文 Table 12 等价宣称。

## 运行与评估

```bash
PYTHONPATH=src python scripts/run_fmo_sensitivity.py \
  --output runs/fmo87/pcm_uff_default_full \
  --gamess-dir "$PWD/software/gamess" \
  --library-dir /path/to/openfe/lib \
  --dftb-parameter-dir "$PWD/software/gamess/auxdata/DFTB/3OB-3-1" \
  --protocol pcm_uff_default --all-cases
```

完成后：

```bash
PYTHONPATH=src python scripts/evaluate_fmo_protocol.py \
  --results runs/fmo87/pcm_uff_default_full/results.csv \
  --protocol pcm_uff_default \
  --output-dir runs/fmo87/pcm_uff_default_full/evaluation
```

只有 `strict_87_complete: true` 的候选才可与 baseline 或论文 Table 12 做完整比较。

几何协议需要使用已经安装 `openmmforcefields` 的生产环境，例如：

```bash
OPENFE=/path/to/ai-vs-physics-ross2023/envs/openfe/bin/python
PYTHONPATH=src "$OPENFE" scripts/run_fmo_sensitivity.py ...
```

否则会得到 `FAILED_GEOMETRY_ModuleNotFoundError`，这属于环境失败，不是 QM 协议失败。

片段边界诊断还表明，关闭 HOP 会破坏奇偶电子/闭壳层约束；关闭 sequence-gap bridging 可产生 126.777 kcal/mol 的 ligand pair。runner 的单 pair 绝对能量质量门已收紧至 100 kcal/mol，fresh rerun 将该 no-bridge 变体 12/12 拒绝，而 baseline 最大 pair 为 28.019 kcal/mol。
