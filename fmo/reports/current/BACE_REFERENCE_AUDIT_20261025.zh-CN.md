# BACE文献与FEP输入核对：2026-10-25

## 本36例的身份

本地发布数据`vendor/public_binding_free_energy_benchmark/fep_benchmark_inputs/benchmark_metadata.csv`将当前36例列为
`jacs_set/bace`、参考4DJW、原方法文献DOI 10.1021/ja512751q，而不是`janssen_bace`。
4DJW的原结构研究为Cumming等，Bioorg. Med. Chem. Lett. 2012, 22, 2444–2449，DOI 10.1016/j.bmcl.2012.02.013；
RCSB登记1.90Å晶体结构，iminohydantoin配体，并列有酒石酸TLA。不能把TLA当作偶然出现就静默删去。

当前制备记录pH7.4不是已确认的实验测定pH。查到的`janssen_bace/bace_retrospective_pka.txt`的pH4.5属于另一套输入，
不得借来当作当前36例的实验条件。`post_simulation_corrections.md`中的Janssen BACE质子校正也不能自动套用给本队列。
本次未取得Cumming论文完整实验部分及补充材料，当前36例的assay pH仍记为未核实，不据此擅改化学态。

## 原FEP方法真正提供的启示

Wang等，JACS 2015, 137, 2695–2703，DOI 10.1021/ja512751q，作者机构存档全文已核对：

- Figure 2所示配体变换分别在受体和溶液两条路径计算，不是只对复合物做更精确的单点。
- 采用REST增强采样；闭环误差主要反映采样问题，并不消除力场误差。
- 文中指出刚性受体评分可能漏掉侧链或配体R基团重排；这支持做匹配结构对照，而非证明某个新能量项必然提分。
- 文中Figure 4C的氯取代排水例子属于thrombin，不是当前BACE，不能误当本队列的已确认病因。

原论文Table 2的BACE相关系数约0.78，和当前本地21-4完整36例重算对照的Pearson约0.624来自不同发布版本，
不可在同一行把前者Pearson和后者Spearman/Kendall拼接。原论文结果也不是本轮重跑。

## 对当前实验的约束

本轮H归一化只回答MM/QM几何模型差异是否污染局部应变，不宣称替代REST、自由能闭环或溶液人口。
后续若研究真正的溶液态参照，必须保持化学态及结合/游离态模型一致，并明确采样与加权，而不是直接给tPIE减一个游离配体总能量。
主流程现有结构处理保持不变，候选只有经过完整覆盖与跨集合回归后才能考虑接入。
