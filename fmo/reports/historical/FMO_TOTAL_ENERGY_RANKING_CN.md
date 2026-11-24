# FMO 总能量双重作差排序（独立探索协议）

## 算的是什么

同一靶点 A、两个配体 B1/B2，在**共享且不移动的**蛋白口袋/固定水/断键封端下分别计算完整 FMO 能量：

\[
\Delta\Delta E_{2-1} = [E_{\rm FMO}(AB_2)-E_{\rm FMO}(B_2)]
                         -[E_{\rm FMO}(AB_1)-E_{\rm FMO}(B_1)].
\]

如另算一次共享受体 A，则各配体有 `ΔE_bind,i = E(AB_i) - E(A) - E(B_i)`；两者的差严格消去同一个 `E(A)`。本实现以 `E(AB)-E(B)` **从低到高**排序，输出相对第一名的能量差（Hartree × 627.5094740631 kcal/mol/Hartree），不读取实验标签参与拟合。受体能量单独计算并输出供审计，不影响排序。

**不是**热力学的 `ΔΔG` 或 OpenFE 的相对自由能：这里只是固定构象 FMO2-DFTB3/3OB-3-1 (+PCM) 的总电子能量作差，不含受体/配体重排、未结合态配体的构象系综、平动/转动/构象熵等。独立配体 B 是从复合物直接剥离的**相同结合构象**，不是优化的自由配体。不要把 kcal/mol 的数值当作实验 ΔΔG 的定量预测；相关性应仅在靶点内按无拟合排名检验。不同靶点的能量基线不可互比。也不应把新实验在旧数据上调参数后称为免拟合。

## 与旧 tPIE 的不可互换之处

旧 `MODMOL=1 MOLFRG(1)=ligand` 只计算配体参与的二聚体，GAMESS 明确不保证完整体系总能量；旧日志的 `Eunco+D(2)` **不能**用于上述作差。本协议对所有 AB 用完整二聚体 FMO2（无 MODMOL），关掉不需要的 PIEDA；孤立单片段 B 用同参数 FMO1（`NFRAG=1 NBODY=1 IPIEDA=0`）；A 用完整 FMO2（如果只有一个片段同样用 FMO1）。三种输入保留相同 SCC-DFTB3、DFTBSK 和 PCM 设置，只接受收敛且正常结束的运行；**启用 PCM 时必须取含溶剂项的 `EunD+so(n)`，不可误用 `Eunco+D(n)`**；仅不启用 PCM 时优先取 `Eunco+D(n)`，无该条则取 `Euncorr(n)+Edisp(n)`（无色散项则仅 `Euncorr(n)`）。拒收选择性 MODMOL 输入及日志复用。

## 结构与可重复性

每靶点先在整个冻结队列的配体姿势上并集选取 5 Å 的残基及可选固定水，再对共享受体做一次单残基断档桥接。保留现有 FACIO 近似、ACE/NME 截断封端与 HOP_C 断键；逐配体检查受体原子、分片、坐标和 HOP 是否完全一致，失败即停。蛋白取 `prepared_proteins/*_prepared.pdb` 的静态坐标，配体取已有 `final_fmo_v3_dftb3_pcm/minimized/*/minimized_complex.pdb` 中的结合姿势，按 SDF 顺序和元素逐一检查；该批原始 pose 通过固定蛋白/水、受约束配体的 OpenMM/OpenFF 最小化得到，仍须检查局部几何和异常大能量。此协议**不是**旧的每配体 5 Å 口袋：旧口袋不同，无法让 `E(A)` 严格相消。`--limit` 只限制实际计算数，口袋仍由该靶点**全队列**决定。

## 操作

在仓库根目录（`$PY` 指向已有环境，`$GAMESS_DIR`/`$GAMESS_LIB` 同旧管线）运行：

```bash
$PY scripts/run_fmo_total_energy.py \
  --only-target CDK2 --limit 2 --prepare-only \
  --gamess-dir "$GAMESS_DIR" --library-dir "$GAMESS_LIB" \
  --dftb-parameter-dir "$GAMESS_DIR/auxdata/DFTB/3OB-3-1" \
  --output runs/fmo87/relative_total_energy

$PY scripts/run_fmo_total_energy.py \
  --only-target CDK2 --limit 2 --timeout 3600 \
  --gamess-dir "$GAMESS_DIR" --library-dir "$GAMESS_LIB" \
  --dftb-parameter-dir "$GAMESS_DIR/auxdata/DFTB/3OB-3-1" \
  --output runs/fmo87/relative_total_energy
```

确认小批量成功、结构合理及预估开销后，去掉 `--only-target`、`--limit` 运行全队列。默认 PCM 和显式水；`--no-pcm` / `--drop-waters` 是**新协议变体**，不得在同一输出目录混用。每靶点 `manifest.json` 记录受体、全部输入 pose/SDF hash、口袋大小、选中对象、参数及几何策略；各 job 的 `input.inp`、`gamess.log`、`result.json` 记录真实 deck、二进制及参数文件哈希，成功缓存仅在这些哈希相同且日志可重解析时复用；`ranking.csv` 是成功样本排序，`failures.json` 记录失败。重跑不同几何、队列顺序、cutoff 或参数时，**使用不同输出目录**以保留审计痕迹。

全量后独立检查每个 GAMESS 正常收敛、deck/可执行文件/参数哈希以及原始日志能量，并与旧 FMO tPIE、Boltz-2、本论文私有 FMO 指标对照；不满 87 例时以非零退出，不生成虚假的完整指标：

```bash
$PY scripts/audit_fmo_total_energy.py --output runs/fmo87/relative_total_energy
```

审计产物 `audit.json` 包含逐靶点和按样本数加权的 Pearson、Kendall、配体对相对能量 MAE；最后一项不等价于自由能方法的热力学相对能量精度，供诊断使用。

全量预计每个靶点一次 A、每个配体一次 AB 和一次 B（总计 4+2×87=178 个 GAMESS 作业）；由于完整 FMO2 包含所有二聚体而非旧选择性计算，并集口袋可能很大，实际运行时间和成功率须先在小批量测量。只有该靶点所有预定样本通过质量闸门时，才报告 Pearson/Kendall 并与原 tPIE/Boltz-2 比较；不得把两例或少数成功样本的排序称为 87/87 新基线。
