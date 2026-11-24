# FMO87 完整复现指南

本文描述从公开结构到 Table 12 对比报告的完整路线。`README.md` 中的短命令只复验已有产物；本文是从零重建的规范路径。

## 1. 外部资产

以下资产不随 git 提交，需按各自许可自行获取：

| 资产 | 验证版本 | 用途 |
|---|---|---|
| GAMESS-US | 2024.2.1 | FMO2-DFTB3/PIEDA 数值引擎 |
| GAMESS 3OB-3-1 | bundled with GAMESS | DFTB3 参数 |
| PDB2PQR | 3.6.1 | 蛋白制备 |
| OpenMM / OpenFF | 8.4.0 / 0.18.1 | 标准几何 |
| RDKit | 2025.09.6 | SMILES、SDF、Crippen clogP |
| Boltz | 2.2.1 | Boltz-2 affinity 预测 |
| Boltz checkpoints | local cache | `boltz2_conf.ckpt`、`boltz2_aff.ckpt` 等 |
| public binding free-energy benchmark | vendored snapshot | 蛋白、配体、实验标签 |

设置路径：

```bash
export PROJECT_ROOT="$PWD"
export PY=/path/to/python
export GAMESS_DIR=/path/to/gamess
export GAMESS_LIB=/path/to/gamess/runtime/lib
export BOLTZ_PY=/path/to/boltz/python
export BOLTZ_CACHE=/path/to/boltz/cache
```

安装核心与测试：

```bash
$PY -m pip install -e .[test]
```

生产环境还需安装 PDB2PQR、OpenMM、`openmmforcefields`、OpenFF Toolkit、RDKit 和 Boltz。当前验证版本记录在 `PROVENANCE.md`。

## 2. 获取公开结构

```bash
git clone https://github.com/schrodinger/public_binding_free_energy_benchmark.git \
  vendor/public_binding_free_energy_benchmark
```

如需从源码研究 Boltz，可另克隆 Boltz 仓库；运行预测只需要安装 Boltz 2.2.1 和 checkpoint cache。

## 3. 冻结 87 化合物队列

```bash
PYTHONPATH=src $PY scripts/build_fmo87_cohort.py \
  --structure-root vendor/public_binding_free_energy_benchmark/fep_benchmark_inputs/structure_inputs \
  --output data/fmo87_cohort.csv

sha256sum data/fmo87_cohort.csv
```

期望 SHA-256：

```text
19edfb8e2ccf3ab1029d848e4062d70cb35532d2b4fee63674f160b2a292997b
```

期望计数为 CDK2 16、JNK1 21、P38 34、TYK2 16。

## 4. 制备四靶点蛋白

```bash
mkdir -p runs/fmo87/prepared_proteins

PYTHONPATH=src $PY scripts/prepare_protein.py \
  --input vendor/public_binding_free_energy_benchmark/fep_benchmark_inputs/structure_inputs/jacs_set/cdk2_protein.pdb \
  --output runs/fmo87/prepared_proteins/CDK2_prepared.pdb --ph 7.4

PYTHONPATH=src $PY scripts/prepare_protein.py \
  --input vendor/public_binding_free_energy_benchmark/fep_benchmark_inputs/structure_inputs/jacs_set/jnk1_manual_flips_protein.pdb \
  --output runs/fmo87/prepared_proteins/JNK1_prepared.pdb --ph 7.4

PYTHONPATH=src $PY scripts/prepare_protein.py \
  --input vendor/public_binding_free_energy_benchmark/fep_benchmark_inputs/structure_inputs/jacs_set/p38_protein.pdb \
  --output runs/fmo87/prepared_proteins/P38_prepared.pdb --ph 7.4

PYTHONPATH=src $PY scripts/prepare_protein.py \
  --input vendor/public_binding_free_energy_benchmark/fep_benchmark_inputs/structure_inputs/jacs_set/tyk2_protein.pdb \
  --output runs/fmo87/prepared_proteins/TYK2_prepared.pdb --ph 7.4
```

默认保留结晶水。P38 输入中的 20 个水进入后续 pocket；其他三个输入没有水。不要加 `--drop-water` 重跑主实验。

## 5. 运行 FMO 几何与 GAMESS

以下命令按靶点生成目录结构。实际生产中可以并行四个 target，但每个 GAMESS job 使用单核。

```bash
for TARGET in CDK2 JNK1 P38 TYK2; do
  PYTHONPATH=src $PY scripts/run_fmo87.py \
    --cohort data/fmo87_cohort.csv \
    --prepared-proteins runs/fmo87/prepared_proteins \
    --output runs/fmo87/fmo_rebuild/targets/$TARGET \
    --only-target $TARGET \
    --gamess-dir "$GAMESS_DIR" \
    --library-dir "$GAMESS_LIB" \
    --dftb-parameter-dir "$GAMESS_DIR/auxdata/DFTB/3OB-3-1" \
    --cutoff 5.0 \
    --timeout 1800 \
    --stage all
done
```

成功标准不是 GAMESS return code 0，而是第 8 步的 strict audit 全部通过。

## 6. 生成 Boltz-2 输入

Boltz-2 使用每个靶点共享 MSA。若已有 MSA，将其放到 `runs/fmo87/boltz_msa/` 并命名为 `CDK2_1.csv`、`JNK1_1.csv`、`P38_1.csv`、`TYK2_1.csv` 后执行：

```bash
PYTHONPATH=src $PY scripts/build_boltz_fmo87_inputs.py \
  --cohort data/fmo87_cohort.csv \
  --output runs/fmo87/boltz_inputs \
  --msa-dir runs/fmo87/boltz_msa
```

如安装了兼容的 MSA 生成工具，可尝试追加 `--generate-msa`。主结果必须记录 MSA 来源和哈希。

## 7. 运行 Boltz-2

```bash
for TARGET in CDK2 JNK1 P38 TYK2; do
  PYTHONPATH=src $PY scripts/run_boltz_fmo87.py \
    --manifest runs/fmo87/boltz_inputs/manifest.json \
    --out-dir runs/fmo87/boltz_predictions/$TARGET \
    --cache "$BOLTZ_CACHE" \
    --python "$BOLTZ_PY" \
    --only-target $TARGET \
    --accelerator gpu --devices 1 \
    --seed 0 \
    --recycling-steps 3 \
    --sampling-steps-affinity 200 \
    --diffusion-samples-affinity 5 \
    --timeout 1800
done
```

若共享一张 GPU 同时启动四个 target，总 wall time 受 GPU contention 影响；论文 Table 12 的 `20 GPU sec` 是方法级单预测参考，不是本地下四进程共享 GPU 的总时间。

## 8. 审计、合并与评估

```bash
PYTHONPATH=src $PY scripts/audit_fmo87_production.py \
  --root runs/fmo87/fmo_rebuild/targets \
  --output-dir runs/fmo87/final --require-complete

PYTHONPATH=src $PY scripts/audit_boltz87_production.py \
  --root runs/fmo87/boltz_predictions \
  --cache "$BOLTZ_CACHE" \
  --output-dir runs/fmo87/final --require-complete

PYTHONPATH=src $PY scripts/merge_fmo87_outputs.py \
  --fmo-targets \
    runs/fmo87/fmo_rebuild/targets/CDK2 \
    runs/fmo87/fmo_rebuild/targets/JNK1 \
    runs/fmo87/fmo_rebuild/targets/P38 \
    runs/fmo87/fmo_rebuild/targets/TYK2 \
  --boltz-targets \
    runs/fmo87/boltz_predictions/CDK2 \
    runs/fmo87/boltz_predictions/JNK1 \
    runs/fmo87/boltz_predictions/P38 \
    runs/fmo87/boltz_predictions/TYK2 \
  --output-dir runs/fmo87/final

PYTHONPATH=src $PY scripts/evaluate_fmo87.py \
  --fmo runs/fmo87/final/fmo_features.csv \
  --boltz runs/fmo87/final/boltz_predictions.csv \
  --boltz-input-manifest runs/fmo87/boltz_inputs/manifest.json \
  --output-dir runs/fmo87/final \
  --bootstrap-repetitions 10000

PYTHONPATH=src $PY scripts/ablate_sophosqm_calibration.py \
  --compounds runs/fmo87/final/compounds.csv \
  --sophosqm-si data/references/sophosqm_si.pdf \
  --output-dir runs/fmo87/final
```

## 9. 只复验既有最终产物

如果 `runs/fmo87/final_fmo_v4_standard_geometry` 和 `runs/fmo87/boltz_predictions_final` 已存在，执行 `README.md` 的 audit/merge/evaluate 命令即可，不重跑 QM 或 Boltz。

## 10. 参考运行成本

- FMO scoring mean runner wall：39.62 s/case；median 36.59 s；max 61.46 s。
- 87 case runner wall sum：57.45 min；四靶点并行实际约 30 min 35 s。
- 该统计不包含蛋白制备、OpenMM 几何、Boltz-2、审计和 bootstrap。
- 本地原 Boltz-2 推理四个 target 共享 GPU，约 3 h 37 m 56 s wall；不能直接等同于论文单预测 `20 GPU sec`。

## 11. 产物边界

`runs/`、`vendor/`、`software/` 被 `.gitignore` 排除。仓库提交代码、冻结 cohort、协议和文档；GAMESS/Boltz 原始日志和预测产物留在本地或另行打包。最终报告的核心表格已固化在 `FINAL_REPORT_CN.md`。
