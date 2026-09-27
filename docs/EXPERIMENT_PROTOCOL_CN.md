# FMO87 / Boltz-2 冻结实验协议

## 1. 精确队列

Boltz-2 附录 D.2 明确四靶点子集来自 protein-ligand benchmark 的 CDK2、TYK2、JNK1、P38，共 87 个中性化合物；D.2.2 规定回归指标先逐 assay 计算，再按化合物数加权。

本地公开 Ross/Chen 四靶点 SDF 共 104 条记录。按 RDKit canonical isomeric SMILES 去重后，17 条 JNK1 alternate-flip 记录是同一化合物的重复 pose；保留 non-flip 代表后得到：

```text
CDK2 16 + JNK1 21 + P38 34 + TYK2 16 = 87
```

选择规则不使用实验标签。冻结文件：

```text
data/fmo87_cohort.csv
SHA256 19edfb8e2ccf3ab1029d848e4062d70cb35532d2b4fee63674f160b2a292997b
```

## 2. 蛋白制备

每个靶点执行：

1. 只保留蛋白和结晶水，滤出配体/非标准记录；
2. PDB2PQR 3.6.1；
3. AMBER force field；
4. PROPKA 3.5.1；
5. pH 7.4；
6. 补氢、修复缺失 OXT、优化氢键网络；
7. 保存 PDB/PQR/log/provenance。

P38 保留 20 个结晶水；其余三个输入没有水。水质子保留进入 pocket/fragment，不由简单删除处理。

## 3. 标准几何

生产 geometry 是 v4：

- OpenMM 8.4；
- 蛋白 AMBER14 `protein.ff14SB`；
- 水 TIP3P；
- ligand OpenFF 2.2.0；
- ligand 电荷 `openff-gnn-am1bcc-1.0.0.pt`；
- 蛋白、水、ligand 均可移动；
- 所有 ligand 原子加 harmonic positional restraint；
- force constant 100 kJ/mol/Å²；
- minimization tolerance 1.0 kJ/mol/nm。


## 4. Fragment / cap / HOP

每个氨基酸 fragment 采用 shifted FACIO-like 方案：

- 当前残基侧链；
- Cα；
- backbone N/H；
- 前一残基 carbonyl C/O；
- 切断 Cα–C′；
- 内部边界写 `$FMOBND ... HOP_C`；
- pocket segment 端点用完整中性 ACE/NME，而不是单 H；
- 不机械修改 donor/acceptor `ICHARG`；
- HOP 后做价电子奇偶校验；
- 全部生产片段为 RHF singlet。

真实蛋白 C 末端同时存在 OXT 时，NME cap 取代 OXT；两者不能共存，否则会出现约 0.15 Å 的 N–O 原子重叠。该修复由单元测试和 `CDK2_1oiu` 真实失败转成功共同验证。

## 5. GAMESS FMO2-DFTB3/PIEDA

生产 deck 使用：

```text
$DFTB SCC=.TRUE. DFTB3=.TRUE. DAMPXH=.TRUE. DAMPEX=4.00
$FMOPRP NPRINT=9 IPIEDA=1 MAXIT=200 CONV=1.0E-7 NGUESS=2 MODORB=3
         CNVDMP=50 MCONV(2)=530 MCONV(4)=530
$DFT DC=.TRUE. IDCVER=4
$PCM SOLVNT=WATER IEF=-10 ICOMP=0 ICAV=1 IDISP=1 IFMO=-1
$FMO NBODY=2 NACUT=0 MODMOL=1
$FMOBND ... HOP_C
```

依据本地 GAMESS 官方输入手册：3OB-3-1 是 DFTB3 生化/水参数；`DAMPXH=.TRUE.` 用于复现发表 DFTB3 行为；PIEDA 需要 `MODORB=3`；自定义 fragment 必须 `NACUT=0`；`CNVDMP=50` 在建议 10–100 范围内。

## 6. 成功门

GAMESS return code 为 0 不充分。每个结果还必须：

- normal termination；
- 无 monomer/dimer SCF unconverged；
- 无 distance/electron/impossible 错误；
- 输出 selected-fragment PIEDA pair table；
- pair-sum tPIE 与 GAMESS selected final TIE 差值 ≤0.005 kcal/mol；
- 无绝对 pair energy >1000 kcal/mol 的非物理 outlier；
- result 中 input SHA-256 与 deck 一致；
- 全部布尔化学/协议审计通过。

## 7. Boltz-2 协议

- model `boltz2`；
- Boltz 2.2.1；
- seed 0；
- recycling steps 3；
- affinity sampling steps 200；
- affinity diffusion samples 5；
- GPU 1；
- 每靶点共享 MSA；
- 输入 ligand 用 cohort canonical isomeric SMILES。

`affinity_pred_value` 官方定义为 `log10(IC50 μM)`；转换为 kcal/mol 使用：

```text
RT ln(10) * (affinity_pred_value - 6), T=298.15 K
```

## 8. 指标

按论文 D.2.2：

- 逐靶点计算；
- 按靶点化合物数加权；
- Pearson R；
- Kendall tau-b；
- pairwise MAE；
- centered/non-centered MAE；
- centered/non-centered 1 和 2 kcal/mol 命中率。

主比较只允许：

- local raw tPIE vs 论文内部 FMO 0.55/0.38；
- local Boltz-2 vs 论文 Boltz-2 0.66/0.48。

`tPIE+clogP` 回归和 LOO 是次要 SophosQM-style 分析，不与 raw FMO 直接比 delta。

## 9. 显式水消融协议

生产协议保留已有 crystallographic water，不新增水。为检验水处理是否解释与 Boltz-2 论文 FMO 行的差异，执行两个互补消融：

1. **Scorer 级去水**：复用保留水体系的 OpenMM minimized complex，只在构建 5 Å FMO pocket/fragment 时移除水，然后重跑 FMO2-DFTB3/PIEDA。该分支隔离 explicit-water scorer 贡献，不改变蛋白制备、质子化或几何。
2. **制备级去水**：在 PDB2PQR/PROPKA 前移除蛋白输入水，重新执行蛋白制备、OpenMM restrained minimization 和 FMO。生产输入中只有 P38 有 20 个 selected waters；CDK2/JNK1/TYK2 为 0 个，因此该分支的结构效应集中在 P38。

两个消融均保留 GAMESS water PCM（`IFMO=-1`）。因此结论只针对 explicit crystallographic water，不针对 implicit solvation。

消融成功门包括：

- `removed_from_fmo_pocket`：全部 FMO pocket 的 water fragment 数为 0；
- `removed_before_preparation`：全部 prepared protein 的 water residue 数为 0，且全部 FMO pocket 无水；
- 两个分支均满足 87/87 normal termination、SCF/PIEDA、TIE reconciliation、哈希和化学协议审计；
- 指标仍按论文 D.2.2 逐靶点计算后按化合物数加权。

## 10. 已声明偏差

- 不使用也不复刻 Boltz 作者私有 FMO scorer；
- Ross/Chen reference poses 替代不可得的 Glide SP poses；
- OpenMM AMBER14/OpenFF 替代不可得的 MacroModel/OPLS3e；
- 自研 DFTB3 为了极性 dimer 稳定使用 water PCM；
- RDKit Crippen clogP 只用于次要校准模型；
- raw tPIE 不代表绝对结合自由能。
