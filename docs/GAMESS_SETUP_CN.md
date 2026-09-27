# GAMESS FMO87 生产配置

## 引擎边界

完整 GAMESS 2024.2.1 源码、`gamess.00.x`、3OB-3-1 参数、MKL/OpenMPI runtime 已在 `engines/gamess/`。不要另行安装或覆盖该目录；如需替换，必须保留相同版本、二进制哈希和参数哈希。

GAMESS-US 由使用方按官方注册条款下载和本地编译。本项目不发布源码、二进制或注册凭据。当前验证使用项目本地 GAMESS `v2024.2.1`。

运行目录必须提供：

- `software/gamess/gamess.00.x`
- `software/gamess/rungms`
- `software/gamess/auxdata/DFTB/3OB-3-1`
- 兼容的共享库路径

## DFTB3 参数

生产方法为：

```text
FMO2-DFTB3 / 3OB-3-1 / PIEDA / water PCM
```

关键 deck 片段：

```text
 $DFTB SCC=.TRUE. DFTB3=.TRUE. DAMPXH=.TRUE. DAMPEX=4.00
 $FMOPRP NPRINT=9 IPIEDA=1 MAXIT=200 CONV=1.0E-7 NGUESS=2 MODORB=3
          CNVDMP=50 MCONV(2)=530 MCONV(4)=530
 $DFT DC=.TRUE. IDCVER=4
 $PCM SOLVNT=WATER IEF=-10 ICOMP=0 ICAV=1 IDISP=1 IFMO=-1
 $FMO NBODY=2 NACUT=0 MODMOL=1 MOLFRG(1)=...
 $FMOBND
 ... HOP_C
```

本地官方 `software/gamess/docs-input.txt` 确认：

- 3OB-3-1 用于 DFTB3 生化和水体系；
- `DAMPXH=.TRUE.` 用于复现发表 DFTB3 行为；
- FMO-DFTB 为闭壳层实现；
- PIEDA 需要 `MODORB=3`；
- 自定义 fragment 必须禁止自动重分区，`NACUT=0`；
- `CNVDMP` 建议范围 10–100，本项目取 50；
- `$FMOBND HOP_C` 是跨片段 C–C 边界轨道修复。

参数文件必须用 job 目录内短相对路径 `../params` 链接到 3OB-3-1，避免 GAMESS 自动参数路径和列宽截断问题。

## 化学分片

生产 fragment 不是整残基、单 H cap 或简化切片：

```text
side chain + Cα + backbone N/H + preceding C=O
cut Cα–C′
internal boundary -> $FMOBND HOP_C
pocket segment terminal -> full neutral ACE/NME
```

真实 C 末端有 OXT 时，NME 取代 OXT。不能同时保留二者；此前 `CDK2_1oiu` 中 LEU296 OXT 与 NME N 距离约 0.15 Å，导致单体能级异常和 dimer SCF 不收敛。修复后该样本正常结束。

## 质量门

`returncode=0` 不算成功。Runner 额外要求：

1. normal termination；
2. 无 SCF unconverged；
3. 无 distance/electron/impossible 错误；
4. PIEDA pair table 可解析；
5. ligand-environment pair 数等于 `NFRAG-1`；
6. pair-sum tPIE 与 GAMESS selected final TIE 差值 ≤0.005 kcal/mol；
7. pair 能量绝对值不超过 1000 kcal/mol；
8. deck SHA-256 与 result 记录一致；
9. 所有片段 multiplicity=1；
10. 生产审计器所有布尔门通过。

当前 87/87 均通过。原始日志和 result 保存在：

```text
runs/fmo87/fmo_rebuild/targets/<TARGET>/fmo/<JOB>/run/
```
