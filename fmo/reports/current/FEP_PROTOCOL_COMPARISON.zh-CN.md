# FEP+对照比当前tPIE多做了什么

核查日期2026-10-21。来源为本地官方benchmark checkout的README、运行与后处理脚本，
不是本项目新运行的FEP/OpenFE，也不是仅凭方法名称推断全部设置。

## 总体差别

当前tPIE：结构制备→固定受体条件下配体k10最小化→共享口袋→FMO2-DFTB3/PCM→跨蛋白配体PIEDA总和。
这是指定制备姿态下的相互作用分数，不是溶液结合自由能。软压缩只改变尺度。

FEP+：在配体变换网络上计算结合态和溶液态的相对自由能差，
涉及环境/构象的统计采样、力场、网络估计和适用的状态后处理。
不是在tPIE上多加一个熵或溶剂数字；FEP也不保证采样收敛或所有熵贡献都准确。

| 项目 | 当前tPIE | 对照FEP+ |
|---|---|---|
| 主估计量 | 单个制备复合物的跨片段相互作用总和 | 配体间相对结合自由能 |
| 构象 | 默认一套k10姿态，受体固定 | 动力学/自由能采样，而非单个最小化点 |
| 溶剂 | PCM连续介质及保留结构水 | 溶液/结合态模拟；发布命令带muVT水采样设置 |
| 参数 | 当前AMBER/OpenFF制备与DFTB3 QM | OPLS系参数、Force Field Builder或已有定制目录 |
| 化学状态 | 保留输入明确H/电荷，未建实验pH下的状态种群模型 | 部分map显式建多个结合模式/状态并作后处理；不是所有集合都需要全部修正 |
| 网络 | 各配体独立打分 | 原始Bennett边估计、cycle-closure corrected估计与节点结果 |
| 实验锚点 | 主分数不拟合实验值 | 发布绝对节点值整体平移到实验均值；不改变P/S/K |

## 对这13个Ross队列实际查到的设置

对应子目录：jacs_set的CDK2/JNK1/P38/TYK2/BACE/PTP1B；waterset的CHK1/HSP90；
misc的Galectin3/CDK8/FAAH/HIV1PROT；opls_stress的FXa。
不要与同名靶点的其它目录混淆，比如jacs_set/TYK2不是charge_annhil/TYK2，
misc/CDK8不是merck/CDK8，jacs_set/BACE不是janssen_bace。

1. 四类目录的公开启动脚本均可见`-time 20000.0 -vacuum -ensemble muVT -seed 1`，
   通常还指定`-ffbuilder`。`-vacuum`不能解读成只做真空结合能，完整工作流仍需按FEP+定义解释。
   这里记录原始参数，不用该值推断每边实际总GPU时间、lambda数、平衡时间或所有REST细节。
2. **Galectin3：** misc脚本为其单独指定`-OPLSDIR galectin3_ffb_oplsdir`；目录含
   `custom_2021_4.opls`与`custom_2022_3.opls`。脚本注释提到22-3 Force Field Builder问题。
   这是公开复现设置存在专门参数目录的证据，**不是已经证明21-4每个原始作业使用了哪项定制参数**。
3. **FXa Yoshikawa：** 专门命令比通用命令多`-modify_dihe`，属于额外二面角相关设置。
   其具体改动原子、历史21-4实际启用情况和贡献大小，须进一步查fmp/作业元数据；不能从开关断言提分量。
4. **JNK1：** jacs_set运行map名含`manual_flips`，文档要求binding-mode correction；
   我们匹配的最终结果名为`jnk1_manual_flips_symbmcorr_out.csv`，属于已处理结果。
   手工设置额外结合模式不等于已证明按实验标签挑最佳姿态；也不应称完全自动单姿态基线。
5. **其它当前队列：** 在此文档中未找到要求它们全部进行同一种额外pKa/互变异构修正的证据。
   不能把merck/CDK8的symmetry correction套给当前misc/CDK8。

### 重要版本限制

当前对照取`21_4_results/ligand_predictions/`，官方README说明这些是构象、旋转/对称性、
pKa/互变异构等适用后处理后的节点预测；具体适用项依赖map。
输入目录的提交脚本包含22-3相关更新，因此脚本不能代替逐个归档21-4作业日志。
本轮核查没有用商业FEP+ API解析所有二进制fmp，不能给出每个case所有底层设置的完整清单。

## 使用实验标签了吗，比较是否公平

README明确说预测值有未知整体常数，发布时将预测均值移到实验均值。
所以这些数字不能直接当完全独立的绝对结合自由能预测；但常数平移不改变Pearson、Spearman或Kendall，
**不能用均值对齐解释FEP+相关性较高。** 网络cycle closure也不是用实验亲和力拟合相对排序。
目前没有证据说明这些集合的相对分数是用实验标签拟合出来的；是否有更早的回顾性方法选择偏差不能仅凭CSV排除。

这是相同化合物/标签下的端到端方法对照，不是相同计算预算、采样或结构制备的单变量比较。
FEP+较高说明更完整的自由能方法在这些系列可有优势，不证明差距全部来自熵或某个特设开关。
OpenBind没有匹配的FEP+参考，本表结论不能扩展成其已有FEP+对照。

## 可复核来源

研究工作区`vendor/public_binding_free_energy_benchmark/`下：
- `21_4_results/README.md`：边、cycle closure、后处理节点与实验均值对齐。
- `fep_benchmark_inputs/fep_plus_inputs/jacs_set/setup_and_run.sh`。
- `fep_benchmark_inputs/fep_plus_inputs/waterset/setup_and_run.sh`。
- `fep_benchmark_inputs/fep_plus_inputs/misc/setup_and_run.sh`：Galectin3单独OPLSDIR。
- `fep_benchmark_inputs/fep_plus_inputs/opls_stress/setup_and_run.sh`：FXa额外modify_dihe。
- `fep_benchmark_inputs/fep_plus_inputs/post_simulation_corrections.md`：按map区分修正。
- `fep_benchmark_inputs/fep_plus_inputs/jacs_set/run_all_corrections.sh`：JNK1结合模式后处理。
