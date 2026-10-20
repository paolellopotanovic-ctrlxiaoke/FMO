# 引擎与数据接入

## 1. 内置资产

```text
data/fmo87_structures/          4 个蛋白 PDB + 4 个 SDF，覆盖 87 case
data/boltz_msa/                 4 个共享 MSA
engines/gamess/                 GAMESS 源码、gamess.00.x、参数和 tests
engines/gamess/runtime/lib/     MKL/GCC 动态库
engines/gamess/runtime/mpi/     OpenMPI 运行时
data/references/                SophosQM SI 与 Boltz-2 论文 PDF
```

`validate_inputs.py` 按文件名把冻结 CSV 中的历史 provenance 路径重映射到 `data/fmo87_structures/`，不需要复刻原始 vendor 目录。

## 2. Python / OpenMM / Boltz 包

```bash
./scripts/bootstrap_engines.sh
```

脚本创建 `envs/fmo87`，安装 `requirements-lock.txt` 中的 RDKit、OpenMM、OpenFF、openmmforcefields、PDB2PQR 和测试依赖，然后 editable 安装内置 `engines/boltz/`（官方 commit `b1ebfc46ecf57f5414e0d1a6f9027bbb122c53bc`）。`run_pipeline.sh` 会优先使用这个解释器。

本流水线不依赖 OpenFE CLI/package；标准几何直接调用 OpenMM/OpenFF，因此不需要复制本机 14 GB 的 OpenFE conda 环境。

如需 conda 方式，`environment.yml` 使用同一份 pip lock；二者不要混用。自定义解释器可用 `PY`/`BOLTZ_PY` 覆盖，但必须满足同一版本合同。

## 3. GAMESS

默认值为：

```text
GAMESS_DIR=$PWD/engines/gamess
GAMESS_LIB=$PWD/engines/gamess/runtime/lib
DFTB_PARAMETER_DIR=$PWD/engines/gamess/auxdata/DFTB/3OB-3-1
```

`rungms` 使用仓库相对路径；Python runner 设置 `OPAL_PREFIX`、本地 OpenMPI、MKL 库和独立 scratch/restart。`gamess.00.x` 无需系统 GAMESS 安装。

便携引擎已在真实 `P38_2ee` FMO2-DFTB3/PIEDA case 上通过：`SUCCEEDED`，tPIE `-83.082 kcal/mol`，wall time `43.70 s`。

## 4. Boltz-2 权重

```bash
./run_pipeline.sh download-boltz
```

官方 `boltz2_conf.ckpt` 与 `boltz2_aff.ckpt` 单文件超过 GitHub 100 MB 限制，因此不作为普通 Git 对象内置。下载器调用 Boltz 官方 `download_boltz2`，再补齐 CCD 字典；推理和 audit 记录 SHA-256。

## 5. 大规模运行

```bash
FMO_JOBS=8 ./run_pipeline.sh fmo-parallel
```

87 case 按 round-robin 分成可断点续跑 shard；每个 GAMESS job 单核。失败时查看 `runs/fmo87/fmo_parallel_logs/`。集群可将 `scripts/run_fmo87.py --shard-index I --shard-count N` 映射为 array job，最后运行 collector。

Boltz 阶段按 target 顺序执行。默认协议为 GPU 1、seed 0、recycling 3、affinity sampling 200、diffusion samples 5；不要在未更新冻结配置和 audit 的情况下改动。
