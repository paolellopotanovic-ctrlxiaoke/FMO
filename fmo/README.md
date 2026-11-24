# FMO：FEP / MD / QM-MM 研究档案

本档案汇总本机已有计算、逐帧与逐重复分析、benchmark 资料和方法学缺口。整理时未新启动 MD、QM 或网络扩展计算。现有结果横跨 2026-09 至 2026-11，不以文件名中的“final”自动判断科学结论已经收敛。

**已完成快速 FMO 排序工具、两条 OpenFE RBFE 变换及多类 MD 诊断；尚未完成收敛的 QM/MM 自由能修正。** 当前最成熟的成果是可追溯的静态 tPIE 流程、CDK8 酰胺制备变形修复、经典 RBFE 首边对照和定位采样/模型问题的证据链。

## 阅读入口和数据分工

| 位置 | 内容 |
|---|---|
| [AUDITED_RESULTS.json](AUDITED_RESULTS.json) | 从最终三重复 JSON 提取的 OpenFE、MD 性质与内部能量汇总 |
| [SOURCE_MANIFEST.json](SOURCE_MANIFEST.json) | 归档文件、本机来源、字节数及 SHA256 |
| [RAW_TRAJECTORY_INVENTORY.json](RAW_TRAJECTORY_INVENTORY.json) | 12 个本机原始 NetCDF 轨迹的清单；是否上传单独记录 |
| [reports/current](reports/current) | 当前维护工程的完整 Markdown、JSON、CSV 报告及证据 |
| [reports/historical](reports/historical) | 历史研究报告；保留原结论与日期，不当成当前默认 |
| [results/openfe](results/openfe) | 首边协议、状态、端态性质、能量、重复稳定性及生产分析结果 |
| [results/qmmm](results/qmmm) | CDK8 多帧、壳层、内部/环境、力验证和真实桥接的结构化结果 |
| [results/main_response](results/main_response) | 13 队列主入口回归与 CDK8 结构响应决定 |
| [benchmarks](benchmarks) | Ross/FEP+、OpenBind、OpenFE Zenodo 元数据及参考结果 |
| [analysis_code](analysis_code) | 本机当前分析源码、工具脚本、测试与环境声明的快照 |
| [reports/INDEX.md](reports/INDEX.md) | 172 篇原始研究 Markdown 的完整导航与同名配套数据 |
| [payloads](payloads) | 最终端态审计 QM 证据的 gzip 分卷，保留逐文件来源哈希 |

数据包含绝对本机路径与原始证据指纹，便于追溯；换机器时须重新配置 QM、参数、制备环境和输入路径。代码快照包含尚未提交到原工程的研究实现，不冒称已发布版本，也不自动覆盖目标仓库的生产代码。

## 1. 项目在哪里

本机工作区根目录记为 `$WORKSPACE`，即 `/inspire/ssd/tenant_predefaa-9a1b-4522-bb10-8850f313be13/global_user/3078-linziqian`。

| 相对本机根目录的位置 | 定位 |
|---|---|
| `projects/FMO-github-release` | 当前维护工具与最新 CDK8 QM/MM 分支；不要被名字误导成“仅发布包” |
| `projects/sophosqm-style-fmo-baseline/FMOtpiefinal_HIV1lab` | 大量实算、OpenFE 原始轨迹、端态与全队列对照 |
| `projects/sophosqm-style-fmo-baseline` | 早期消融、外部 benchmark、QM 引擎与历史分析 |
| `projects/fmo87-pipeline-tool` | 早期工具及冻结参考结果；公式可能与现行主流程不同 |
| `ai-vs-physics-ross2023` | 2026-09 的独立冻结公开 benchmark / OpenFE / FMO 准备工作区 |

最后一个工作区的 `results/` 目前只有占位文件；其 README 的“未跑模拟”只描述该独立工作区，不能推出整个本机未做 FEP。实际 CDK8/BACE MD 在上表第二个目录。

## 2. 已下载的 benchmark 与尚未下载的材料

| 材料 | 本机状态 | 能支持的比较 |
|---|---|---|
| Ross 2023 public binding free energy benchmark | checkout、输入、公开 FEP+ 结果和校验清单已存在；Gate 0 记录 13 类、92 个 system、1343 配体记录，其中 1340 可评价 | 可以重建公开系列并核查 FEP+ 参考；1343 不是本机已计算的配体数 |
| 当前 FMO 对照子集 | 13 个 Ross 队列共 270 配体；原 87 例是 CDK2/JNK1/P38/TYK2 四靶点子集 | 已有完整静态分数与模型比较；不是 270 例都新跑了 FEP |
| OpenBind EV-A71 2A | benchmark、实验标签、多模型公开预测已下载；FMO 覆盖 639 晶体、聚合至 490 化合物 | 有 FMO 与公开 Boltz2 对照；未找到匹配的 FEP+ 全集参考 |
| protein-ligand-benchmark | 仓库和数据目录已存在 | 数据源之一，不等价于已复现 ABFE |
| OpenFE 2024 Zenodo | 38 条记录的 DOI、文件名、大小和 checksum 元数据已下载 | 约 220.74 GB payload 未整体下载；元数据不是可直接分析的原始轨迹 |
| ABINIT-MP 官方材料 | 示例/工具压缩包与解包内容已存在 | 官方样例不代表具有可直接执行的完整 QM 引擎 |
| Boltz2、SophosQM、AnewFEP 等 | 本地有论文/SI/方法审计 | 方法参考；不能把论文指标称为本机重跑结果 |

详见 [Ross 元数据](benchmarks/ross2023_metadata)、[Zenodo 清单](benchmarks/openfe_zenodo_metadata/manifest.tsv)、[FEP+ 参考](benchmarks/fep_plus_21_4) 与 [OpenBind](benchmarks/openbind)。

## 3. 静态 FMO 做了什么，当前结果是什么

当前主流程：受体/配体身份及化学图审计 → PDBFixer / PDB2PQR / PROPKA → 水氢处理 → k10 最小化 → 共享口袋/分片 → GAMESS FMO2-DFTB3/PCM → 收敛、覆盖率与救援审计 → 完整 tPIE 与单调软压缩。

已做的研究包括 PIEDA 五项及其 31 个非空组合、口袋/截断/封端、显式水、水氢、k 值与姿态、受体侧链网络、化学图与质子化、分子量/大小偏差、三态总能差、PBSA、PCM 溶剂替换、配体局部应变及跨靶点迁移。详细正负结果均保留在 reports，而不是仅保留最高指标。

**版本差异必须保留：** 2026-10-20 的模型总表里 CDK8 仍是旧基线 P/S/K = 0.274/0.236/0.156；2026-10-21 有限酰胺响应接入后为 **0.707/0.697/0.511**。修复针对 k10 新诱发的酰胺变形，只按几何和整队列 QM 验收接纳，属于已知开发集上的回顾性改进。

| 集合 | n | 当前静态 FMO P/S/K | 公开 FEP+ 21-4 P/S/K |
|---|---:|---|---|
| CDK2 | 16 | 0.868/0.856/0.683 | 0.649/0.623/0.444 |
| JNK1 | 21 | 0.614/0.628/0.473 | 0.841/0.877/0.711 |
| P38 | 34 | 0.503/0.481/0.329 | 0.702/0.723/0.543 |
| TYK2 | 16 | 0.708/0.729/0.550 | 0.928/0.929/0.812 |
| Galectin3 | 26 | 0.562/0.544/0.387 | 0.874/0.848/0.661 |
| FXa | 27 | 0.613/0.569/0.394 | 0.751/0.684/0.508 |
| CHK1 | 13 | 0.464/0.443/0.348 | 0.785/0.765/0.606 |
| CDK8，酰胺响应后 | 10 | **0.707/0.697/0.511** | 0.663/0.648/0.378 |
| BACE | 36 | 0.212/0.168/0.112 | 0.624/0.557/0.418 |
| HSP90 | 11 | 0.749/0.725/0.598 | 0.829/0.651/0.524 |
| PTP1B | 23 | 0.745/0.758/0.615 | 0.921/0.875/0.669 |
| FAAH | 24 | 0.737/0.702/0.541 | 0.849/0.656/0.505 |
| HIV1PROT | 13 | −0.147/−0.047/0.039 | −0.006/0.055/0.039 |
| OpenBind，490 化合物 | 490 | 0.371/0.441/0.288；软 P=0.416 | 未找到匹配参考 |

静态 FMO 初始覆盖记录为 14 队列、909 结构（270+639），不代表 909 个独立化合物。FEP+ 列是下载的商业方法发布结果，绝不是本机 OpenFE。不同方法的结构制备、采样与预算不同，点估计没有自动构成显著性或前瞻盲测证据。

局部应变不是通用补丁：CDK8 加应变 P 可到 0.800，Galectin3 和 BACE 也有局部改善，但 CDK2、CHK1、FAAH 明显恶化；FXa/PTP1B/OpenBind 的该轮应变试验还未全覆盖。三态总能差在当前 CDK8 上 P 仅约 0.123；更完整的端点表达式不自动产生更好的自由能预测。

来源：[主入口回归](reports/current/MAIN_RESPONSE_AND_BACE_RESULTS.zh-CN.md)、[模型总表](reports/current/MODEL_COMPARISON.zh-CN.md)、[全队列应变](reports/current/LOCAL_STRAIN_ALL_DATASETS_RESULTS_20261025.md)。

## 4. 本机真正做了哪些 FEP

真正执行的是 OpenFE 经典 MM RBFE：CDK8 的 12→20，以及 BACE 的 CAT-13n→CAT-4i，各自只做预先固定的第一条边 edge_000。

每条边有 bound/solution 两腿、每腿 11 λ 窗口、3 个独立重复，每窗口 1 ns 平衡 + 5 ns 生产；名义积分采样量 396 ns/边、两边共 792 ns。该数不是两个体系各跑了 396 ns 连续物理时间，也不是实测 GPU 忙碌时长。温度 298.15 K，源协议为周期 PME、含氢键约束、HMR 3 amu、4 fs；细节以冻结协议为准。

最终状态文件已是 COMPLETE，两边两腿均已完成三重复；早期 interim Markdown 中“还在等待”是历史快照。

| 变换，方向 B−A | 实验 ΔΔG | OpenFE 三重复均值 | 三次 bound 减共同 solution 均值 |
|---|---:|---:|---|
| CDK8 12→20 | −0.450 | **−0.568** | −0.570、−0.638、−0.496 |
| BACE CAT-13n→CAT-4i | +0.280 | **+2.514** | +2.019、+2.724、+2.800 |

单位 kcal/mol。上表三次投影共用 solution 参照，不能把组合视为完全独立实验；bound SD 分别约 0.071/0.431，solution 重复 SD 约 0.033/0.087，也不是完整模型误差。

CDK8 首边接近实验，BACE 高估约 2.23 kcal/mol；不能由一个配体对推导全队列 P/S/K。没有扩展 CDK8 16 边+BACE 58 边的完整 74 边网络；按同协议名义采样约 29.3 μs，原研究因成本未扩展。

证据：[最终结果快照](results/openfe/final_endpoint_audit_20261023/free_energy_snapshot.json)、[最终状态](results/openfe/final_endpoint_audit_20261023/status.json)、[FEP+ 协议比较](reports/current/FEP_PROTOCOL_COMPARISON.zh-CN.md)。

## 5. MD 相关性质分析完成了什么

从交换后的真实热力学端态取帧，而不是把固定 replica 误当成固定配体；排除虚原子、核查非恒等 analysis_particle_indices、按化学键图展开周期边界。最终各腿/端态/重复取生产 1–5 ns、每 0.1 ns 一帧，共 984 个端态观测；这些不是 984 次独立重复。

已分析：自身对齐内部 RMSD、局部受体对齐姿态 RMSD、重原子非质量加权回转半径、蛋白–配体最小距离、残基接触占比、侧链 χ1 与圆统计、配体杂原子 3.5 Å 水接触数、距离判据水桥、逐帧接触与 PIEDA 的联系、早晚窗与独立重复稳定性。

| 配体 | 内部 RMSD：bound/solution，Å | 姿态 RMSD，Å | 水接触：bound/solution |
|---|---:|---:|---:|
| BACE CAT-13n | 0.712/2.103 | 0.982 | 2.821/9.033 |
| BACE CAT-4i | 0.635/1.978 | 1.327 | 3.528/8.602 |
| CDK8 12 | 0.340/1.483 | 0.742 | 1.455/8.211 |
| CDK8 20 | 0.512/1.455 | 0.850 | 0.732/6.691 |

数值由最终三重复 complete_20261026 JSON 的逐重复均值等权汇总，故与早期只含 1–2 条 bound 重复的 Markdown 不完全相同。

主要发现：游离态内部构象变化更大；BACE 的 ASP93/ASP289 关键近接触总体保持，不能把偏差归结为配体离开催化位点。CDK8 LYS52 的距离和 PIE 随重复明显变化，即使 χ1 都被归为 trans；单个粗粒度转子分类不足以概括协同几何。

这些是结构/平衡分布性质，不是动力学速率常数。当前档案没有形成 k_on/k_off、驻留时间、MSM、扩散系数、长时间蛋白 RMSF/PCA/自由能面或完整结合熵的已验收结果。距离水桥不是含角度判据的氢键，也不是水驻留寿命或去溶剂化自由能。

证据：[完整性质 JSON](results/openfe/endpoint_structure_diagnostics_complete_20261026.json)、[端态分析方法](reports/current/OPENFE_ENDPOINT_DIAGNOSTICS.zh-CN.md)。

## 6. 同帧能量及溶剂研究做了什么

### 孤立配体内部能量

完成 984/984 个最终固定帧 DFTB3 单点；MM 则分解 bond、angle、torsion、Coulomb、LJ，核对总能闭合与原子对静电重构。这里比较的是相同 MM 采样构象在 MM/DFTB3 下的内部势能，不是重新采得的 QM 平衡分布。

| 三重复最终内部能双差 | CDK8 | BACE |
|---|---:|---:|
| MM | +0.817 | +1.401 |
| DFTB3 | −1.209 | −0.144 |
| DFTB3−MM | −2.026 | −1.545 |
| 模型差重复标准误诊断 | 0.320 | 1.093 |

单位 kcal/mol。BACE 重复与早晚时间窗敏感；不能将 −1.545 直接加到 +2.514 充当“QM 修正后的 FEP”。早期 738 帧结果是中间快照，完整结果见 [repeat_stability.json](results/openfe/final_endpoint_audit_20261023/repeat_stability.json)。

### 复合物同帧 FMO 与水

完成 492/492 项结合态同帧 FMO，36 项复用、456 项新增；每体系两端态×三重复×41 帧。tPIE 的 B−A 均值：CDK8 −3.870（重复 SD 1.966）、BACE +4.134（SD 2.343）。BACE 三帧时曾为 +0.451，加密后偏离，说明“偶然贴近实验”经不起取帧变化。BACE 早晚窗差约 6.525 kcal/mol，尚不能宣布收敛。

另做保持蛋白/配体不变的局部水替换、覆盖与水身份映射审计。只换水也会改变蛋白片段 PIE；显式水与 PCM 响应存在耦合，不能在旧 tPIE 上简单加一个独立水修正。

证据：[加密同帧 FMO](reports/current/DENSE_SAME_FRAME_FMO_RESULTS.zh-CN.md)、[水对照](reports/current/LOCAL_WATER_FMO_CONTROL_RESULTS.zh-CN.md)、[内部能定义与旧快照](reports/current/ENDPOINT_INTERNAL_ENERGY.zh-CN.md)。

### 熵代理

历史 E55 的显式溶剂扭转熵代理有六条轨迹，但分半 JS 约 0.412/0.289，高于 0.15 门槛，未通过稳定性。不把 RMSD、水数或未收敛二面角边缘分布当作可补加的结合熵。历史默认曾为 ES+SOLV，现行主流程是完整 tPIE；两时期的“最终结论”不能混用。

## 7. CDK8 QM/MM 自由能修正推进到哪里

研究对象始终只有 edge_000 的 bound/solution×A/B 四端态。局部 QM 区为配体，蛋白/环境为固定电荷 MM；不是整蛋白 FMO 或全周期可极化 QM/MM。

1. 三重复×五固定时点，共 60 对 QM/MM 与匹配 MM 单点、120 个 GAMESS 作业。发现 GAMESS QM-center 惩罚项污染，关闭后能量差严格闭合。安全设置 SWRB=6000 Å，而不是旧试跑 1001 Å。
2. 四态双差均值由 −1.598、SD 8.244 变为 −3.203、SD 3.595 kcal/mol。8/12/16 Å 壳层双差均值约 −3.458/−3.203/−3.368；壳层改动不是均值异常的唯一解释，也不恢复周期长程静电。
3. 总修正重加权 ESS 只有约 1.71–2.33/15。内部/环境双差 −1.303/−1.901，二者有协方差；不能删内部项来美化 ESS。
4. 全部 60 构型仅氢松弛，另复算 60 个原基线，共 120 QM；内部双差 SD 从 3.810 增到 4.134，未支持“MM 松一下氢即可解决”。
5. 完整嵌入力四态 56 项双步长差分核验通过，最大误差 0.004044 kcal mol⁻¹ Å⁻¹；组合周期 MM 时还解释并核验原参考势硬截断跳变。
6. 已真实采样 λ=0/0.05/0.10：四态各窗 4 次预热+12 次生产提案，总计 192 提案、196 个新 QM/MM、144 条含拒绝驻留的生产记录、60 个完整周期快照。不是仅在旧帧上重算 λ。

真实桥接的关键负结果：bound/0 在 λ=0.10 的生产接受为 **0/12**，bound/10 为 1/12；溶液态分别 3/12、10/12。bound/0 记录 12 个相同驻留状态时，权重 ESS 仍可能接近 12，不能据此认为有 12 个独立样本。

匹配提案诊断显示结合态主要先卡在 MM 提案本身，而 solution/0 有明显 QM 修正接受代价；步长减小到 0.1 fs 的有限对照并未消除结合态问题。

目标明确为 `U_MM,PME + λ(U_QMMM,shell − U_MM,shell)`。实际传播是受约束 MM 提案 + QM/MM Metropolis 能量校正，并非每步 QM 力驱动的 MD。尚未到 λ=1，也未解决本轮固定盒 NVT 与源 NPT FEP 的循环一致性。

证据：[最新桥接结论](reports/current/CDK8_QMMM_BRIDGE_CONCLUSIONS_20261101.zh-CN.md)、[力验证](reports/current/CDK8_QMMM_EMBEDDED_FORCE_VALIDATION_20261101.zh-CN.md)、[逐态采样](reports/current/CDK8_QMMM_BRIDGE_SAMPLING_20261101.zh-CN.json)。

## 8. 当前缺口与可验收的后续工作

| 优先级 | 缺口 | 完成的判据 |
|---|---|---|
| 1 | 结合态 MM 提案效率不足 | 在同一目标势下验证可逆性/详细平衡；比较独立有效样本每次 QM 调用，不能仅看原始权重 ESS |
| 2 | QM 修正态连接未完成 | 到 λ=1；足够独立正反向重复、时间相关与重叠/滞后检查 |
| 3 | Hamiltonian/系综一致性 | 明确有限壳层局部替换与全周期 QM/MM 的区别，处理边界连续性及 NVT/NPT 循环 |
| 4 | BACE 化学状态与机制未解决 | 实验 pH、配体/催化微观状态及其权重、力场/扭转与溶剂响应的可分离对照 |
| 5 | 单边诊断不能证明排序提高 | 冻结协议后增加多个边/完整网络，报告覆盖、重复误差、成本及跨靶点验证 |
| 6 | 应变并非普适修复 | 补齐未完成队列，并对 CDK2/CHK1/FAAH 的恶化给出一致机制；不逐靶点选公式 |
| 7 | 全结合熵与水交换自由能未得出 | 使用明确定义、充分采样的统计力学估计，防止与已有采样/溶剂项重复计数 |
| 8 | 文档中的历史状态易混淆 | 以本页+最终状态 JSON 导航，保留旧快照和协议版本，不悄悄覆盖历史 |

当前不应输出：收敛 QM/MM 四态自由能修正、QM 修正后的新全 CDK8/BACE P/S/K、已经超过 FEP 的结论、已复现完整 ABFE/RXRX/Boltz-ABFE2 的声明。

## 9. 归档完整性与大文件边界

本档案完整复制当前 docs 和历史研究 docs，并收录结构化结果、可用参考 CSV、协议、分析源码与来源哈希。最终端态审计的 QM 子目录另完整收录为分卷压缩证据包，每卷小于 50 MiB；执行 `python verify_archive.py --reconstruct` 可重建并逐项核验，不自动解包或运行模拟。仅校验运行 `python verify_archive.py`。

它是研究分析档案，不是原运行目录的逐字节全盘镜像。未收录的其他原始 QM 日志/输入、生产 checkpoint、全部中间结构与软件二进制依旧可按本机来源追溯；复制范围由 SOURCE_MANIFEST 精确列出。OpenFE 相关完整本机运行树约 38.21 GB，其中日志约 29.51 GB、checkpoint 约 5.57 GB。不要将本档案称为这 38.21 GB 的全部原始文件。

已确认原始生产 `simulation.nc` 共 12 文件、约 2.70 GB。是否上传原始轨迹及大文件目的地需单独确认；未上传时明确标识，不用元数据冒充轨迹已上传。OpenFE Zenodo 的 220.74 GB 远程原始 payload 同样不在本归档内。

用户已确认将本档案上传到公开仓库 `paolellopotanovic-ctrlxiaoke/FMO` 的 `fmo/` 目录。仓库可见性保持公开；归档包含研究报告、结构化结果、分析源码及指定证据分卷，原始大轨迹仍按清单留在本机。
