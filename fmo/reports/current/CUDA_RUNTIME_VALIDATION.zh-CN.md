# 本机OpenFE GPU运行环境验证

日期：2026-10-22。两张NVIDIA H800，驱动550.144.03，系统CUDA12.8。

## 已定位与验证的环境问题

项目Python环境包含NVRTC12.9；首次OpenFE CUDA试跑在创建GPU上下文时出现 `CUDA_ERROR_UNSUPPORTED_PTX_VERSION (222)`。这次失败不能解释成采样不收敛或结构化学错误。保留原始ProtocolDAGResult和错误记录，不覆盖为成功。

通过进程局部加载系统CUDA12.8 NVRTC修复，无需修改驱动或共享Python环境：

```bash
export LD_PRELOAD=/usr/local/cuda-12.8/lib64/libnvrtc.so.12
export LD_LIBRARY_PATH=/usr/local/cuda-12.8/lib64:${LD_LIBRARY_PATH:-}
```

在该环境下，OpenMM CUDA平台已成功创建、执行自定义势能GPU内核并完成10步积分。`/proc/self/maps`确认实际加载系统12.8.61的NVRTC及builtins。此环境修复只证明GPU内核可以执行，不证明OpenFE两腿和自由能收敛。

## 当前试运行

研究目录 `FMOtpiefinal_HIV1lab/runs/matched_openfe_20261022_cuda128`：GPU0运行CDK8，GPU1运行BACE；每张卡按溶液腿→结合腿顺序执行首条边的工程短跑。生产5ns/3重复及全网络须在两腿试验、输出审计和成本测量后安排，不把短跑估计值作为亲和力预测。

CPU已完成两组溶液腿短跑；原CPU结合腿中止并记录 `INTERRUPTED_FOR_GPU_MIGRATION`，不是完成。首次CUDA失败在 `matched_openfe_20261022_cuda`，修复后独立新目录保留输入与执行证据。

该配置是当前主机特定的兼容性处理，不写入通用包的默认环境；其他主机先运行CUDA内核检查，不能无条件预加载本机路径。高层级RHF/MP2对照仍使用CPU，不宣称已GPU加速。
