# 缺失六组Boltz2补跑协议

2026-10-20冻结：CHK1 13、CDK8 10、BACE 36、HSP90 11、PTP1B 23、FAAH 24，共117配体。

- 使用最新FMO全量运行登记的原始受体与无标签配体SMILES清单，整组保留，不挑病例。
- 原始受体中所有连续蛋白链转换为序列；CDK8保留A/B两链，不拼成虚构单链。
- 使用先前对照的Boltz2结构/亲和力checkpoint、CCD和相同运行依赖，冻结SHA256。
- seed=0、recycling=3、结构200步/1样本、affinity 200步/5样本；不换参数挑结果。
- 每个受体生成共享、查询序列匹配的真实多序列MSA；多链使用greedy配对，不静默退回单序列。
- 不提供晶体坐标、模板或FMO分数；Boltz2自行生成复合物，不是晶体姿态重打分。
- 独占本批任务的单GPU串行队列；不干扰另一块已使用GPU。每组MSA和推理最多各两次同协议尝试。
- 不成功的记录保留日志；任一组不完整则不生成冒充全六组完成的补充结果。
- 自动验证全队列身份、输入/MSA/模型哈希、GPU成功日志、亲和力JSON、复合物链长度与配体重原子数。
- 所有六组验收通过后才读取实验标签，自动更新`MODEL_COMPARISON.json`、`MODEL_COMPARISON.csv`与中文比较表。
- 记录MSA与推理尝试墙钟，包括启动/加载和重试；缓存命中时间为零不冒充真实GPU成本。

启动前先冻结：

```bash
python tools/run_missing_boltz2.py --root RUN_DIR --workspace RESEARCH_WORKSPACE \
  --cache BOLTZ_CACHE --gpu 1 --prepare-only
```

后台执行`python tools/run_missing_boltz2.py --root RUN_DIR`即可；它同步等待子任务，
完成后自行验收和汇总，不需要交互代理反复轮询。
状态文件为`RUN_DIR/status.json`，总日志为后台启动时指定的`supervisor.stdout`。
每组保留`msa_attempt*.log`、`inference_attempt*.log`及逐例GPU日志。

后台完成不等于自动提交Git；报告与研究输出保存在本机，待核查后再发布。
