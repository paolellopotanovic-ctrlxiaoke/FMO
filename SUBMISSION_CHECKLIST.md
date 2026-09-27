# 私有仓库提交检查

## 提交前

- [x] 仓库保持 private/internal，不向机构外暴露 GAMESS 或参考 PDF。
- [x] 87-case PDB/SDF、MSA、cohort 和源哈希全部内置。
- [x] GAMESS 可执行文件、参数、MKL/OpenMPI runtime 内置。
- [x] 便携 GAMESS 真实 FMO smoke 通过。
- [x] `requirements-lock.txt`、`environment.yml` 和 bootstrap 脚本可重建 Python/OpenMM/Boltz。
- [x] Boltz 权重不进入 Git，由官方 URL 下载并做哈希审计。
- [x] `runs/`、`envs/`、Boltz cache 不进入 Git。

## 提交后

- [ ] fresh clone 后执行 `./scripts/bootstrap_engines.sh`。
- [ ] 执行 `./run_pipeline.sh verify` 和 `check-inputs`。
- [ ] 有 GPU 机器执行 `download-boltz` 和一个小 target 的 Boltz smoke。
- [ ] 重算全量前备份 `reference_results/final/`。
- [ ] 不把结果描述为 Boltz 作者私有 FMO 的同代码复现。
