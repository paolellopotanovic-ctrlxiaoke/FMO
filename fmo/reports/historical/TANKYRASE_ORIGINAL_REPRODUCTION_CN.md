# Okimoto 2018 Tankyrase 2：原论文等价复现审计

日期：2026-09-28。原文：Okimoto 等，*ACS Omega* 2018，DOI `10.1021/acsomega.8b00175`，附录 `ao8b00175_si_001.pdf`。本页区分**原始已发表计算结果的逐数值复现**与**从原始三维结构重新运行整套计算**，绝不将前者称为后者。此前的四靶点 87 例 PCM/DFTB3 + PBSA 混合代理见 `docs/TANKYRASE_PBSA_FOLLOWUP_CN.md`，不是本论文复现，也不能用来验证本论文的 `R=0.856`。

## 原方法与原队列

- **原队列**：同一靶点 Tankyrase 2，共 23 个同系列配体，A 14 个、B 9 个。标签来自文献 IC50，作者使用 `ΔG_exp=RT ln(IC50/M)`、`T=298 K` 转换；IC50 并非直接测出的平衡结合自由能或 Kd。附录 Figure S1/S2 记录各配体 IC50 与旋转键数，Table S3 记录 23 例各项能量和评分。
- **结构**：C01 的 Tankyrase 2 晶体模板 `3KR8`；其他配体依照 Tankyrase 1–C04 模板 `4KRS` 建模；B 系列用 MOE LowModeMD 生成多个姿势。AMBER12 中以 AMBER99SB/GAFF、Gaussian09 RHF/6-31G(d) RESP 电荷和 TIP3P 水进行局部限制 MM 最小化；去水后只放开配体、完整 Gly1032 和 Ser1068 侧链，Gaussian09 的 B3LYP/6-31G(d):UFF ONIOM 精修，其他原子固定。B 系列的多姿势选择影响结果；不能仅用任一 PDB 静态口袋代替。
- **气相 QM**：同一固定几何分别算复合物、游离蛋白和游离配体的 GAMESS FMO 总能量，再做 `ΔE_gas=E_gas(AB)−E_gas(A)−E_gas(B)`。论文使用 FMO4.1 多层 FMO2，低层 HF/STO-3G，高层配体和约 6 Å 内 35 残基 MP2/6-31G(d)；文字末尾另有 `MP2/6-31G(s)` 的不一致简写，**具体轨道基组应以原始输入文件核实**。这不是 PIEDA 的跨片段 tPIE，也不是 FMO2-DFTB3。
- **溶剂与熵**：对固定几何三状态分别用 AMBER12 PB/SA，`ΔG_solv=G_solv(AB)−G_solv(A)−G_solv(B)`；PB 内/外介电 1/80、0.5 Å 网格、Tan–Luo 半径及 AMBER 电荷；非极性 Molsurf SASA、探针 1.4 Å、`SURFTEN=0.00542` 与 `SURFOFF=0.92`。最终 `ΔG_score=ΔE_gas+ΔG_PB,polar+ΔG_SA,nonpolar+α N_rot`，比较原文指定的 `α=0` 和 `1 kcal/mol/rotor`；**没有用 23 例标签拟合回归系数**。单几何与旋转键代理并未求出严格统计力学自由能。
- **指标**：原文 Pearson `R`、Pearlman–Charifson `PI`（按实验能差加权的两两顺序一致度），均为同队列回顾性评估；没有在额外 Tankyrase 外部盲测系列上验证。强相关不代表绝对能量数值正确；原论文亦指出绝对预测明显偏负。

## 原始发表数据的数值复现：完成

附录 Figure S1/S2 与 Table S3 已转录到 `data/tankyrase2_okimoto2018_table_s3.csv`；Table S1 到 `data/tankyrase2_okimoto2018_table_s1.csv`。原附录 PDF SHA-256：`5feb98d4c19bb263d8034922d3c7826b1fb367c559de0682b5727a855b80b8de`（本文不在仓库内分发该 PDF）。脚本**仅使用作者已经发表的能量及 IC50/旋转键**，没有将作者的量子数据说成自己的重算结果。

```bash
python scripts/reproduce_tankyrase_okimoto2018.py
pytest -q tests/test_tankyrase_okimoto2018.py
```

输出 `runs/tankyrase2_original_si/summary.json` 和逐个配体的 `per_compound.csv`。在各自舍入精度下，重复作者公式及评估如下：

| 23 例全部配体 | 本地原表逐数值重算 | 论文 Table S2/S3 |
|---|---:|---:|
| 气相 MFMO 单项 Pearson | 0.70350 | 0.704 |
| PB 极性单项 Pearson | -0.03972 | -0.040 |
| SA 非极性单项 Pearson | 0.47945 | 0.479 |
| 气相 MFMO + PB/SA，α=0，Pearson / PI | 0.85541 / 0.87335 | 0.856 / 0.874 |
| 气相 MFMO + PB/SA，α=1，Pearson / PI | 0.77028 / 0.80000 | 0.770 / 0.800 |
| A 系列 α=0，Pearson | 0.89154 | 0.892 |

表 S3 的三个输入项及公布总分均已截断到 0.001 kcal/mol，单例加和最多差 0.001；实验 `ΔG` 截到 0.01 kcal/mol。因而 PI 的 `0.87335` 与原文 `0.874` 的微差**不能通过改参数来消除**。IC50 转换和公布 ΔG 最大偏差约 0.0061 kcal/mol，吻合舍入量级。

原 Table S1 的 **14 例常规 FMO + PCM** 原分项消融也逐项重现：气相 FMO `R=0.68917`；PCM 极化电子项 `R=0.60740`；仅加 PCM 极性溶剂 `R=0.86977`；再加全部 cav/disp/rep 非极性项 `R=0.50743`，均吻合论文表内三位小数。这是**不同于**主方法 MFMO+PB/SA 的 14 例消融，不能据此断言 PCM 主方法全队列优于 PB/SA。表 S1 其余六种组合的回归也由脚本核查。

## 从头计算的等价复现：尚不能完成

作者附录提供能量与配体二维结构，**没有提供**全部 23 个经过 MOE→AMBER→ONIOM 的最终三维复合物、各 B 配体全部姿势/筛选规则的逐例输入、FMO4.1 三状态输入与残基高层列表、GAFF RESP 电荷/拓扑、Tan–Luo 半径的逐例使用文件或 AMBER12 PB/SA 输入输出。公开晶体 `3KR8`/`4KRS` 仅是**建模模板**，不是 23 个终态几何；由它们重建会得到**另一种实现**，不能用 23 个发表的分项能量校验称作独立重算。

本地有 GAMESS 与 AmberTools `pbsa` 等程序，但本项目现用的 FMO-DFTB3 与多层 MFMO **不兼容**（见 `software/gamess/docs-input.txt` 的 FMO-DFTB 限制）；必须转成非 DFTB HF/MP2 MFMO。当前未找到 Gaussian09 与 MOE 可执行文件；不能用 OpenMM/其它 QM 软件替代并仍宣称逐步骤等价。特别是此前 87 例的气相 SCF 不收敛或 PCM 极化电子能量、PDB2PQR/OpenFF 半径电荷均**不能**冒充本论文 23 例的气相 MFMO+AMBER12 PB/SA。

若要进一步完成**独立从头等价重算**，先向作者取得或合法重新产生 23 个全部姿势及最终 ONIOM 几何、每态的 GAMESS 输入和收敛日志、原 Amber 拓扑/RESP 电荷与三态 PB/SA 设置、原构象筛选记录；确认原论文基组记述与 Molsurf/半径参数。按 C01 起步逐项核对三态原始能量和能量差，核对 PB 极/非极项至作者的舍入精度，最后扩展到 23/23 并独立确认实验标签、α=0/1 的 R/PI。**在这些关键材料可用之前，本项目只能宣称已完成原论文公开结果的数值再分析，而不能宣称原始物理计算等价复现。**
