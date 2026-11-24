# FMOtpiefinal：独立发布候选版

新增 [FEP / MD / QM-MM 完整研究档案](fmo/README.md)，汇总已完成计算、
benchmark、动力学轨迹性质分析、结果与未完成项；全部原始研究报告见
[172 篇报告索引](fmo/reports/INDEX.md)。档案内分析代码为研究快照。

这是从研究工作区中分离出的完整主流程代码包，名称和命令均为
`FMOtpiefinal` / `fmotpiefinal`。不要把外层研究目录一起上传。

本 GitHub 仓库根目录就是当前工具；之前的完整旧工程保存在 `archive/`，
不要再用历史启动脚本作为当前主流程。安装后先运行 `fmotpiefinal-doctor`，
也可加上三个 QM 路径参数进行依赖及路径检查，再启动 `blind-run`。

## 主流程

原始受体与已对齐配体姿态 → 输入/组分/链审计 → PDBFixer →
PDB2PQR/PROPKA → 图审计及已验证的自动修复 → 水与氢处理 →
OpenMM k10 → 全队列共享口袋与分片 → GAMESS FMO2-DFTB3/PCM →
收敛/覆盖率检查及注册救援 → 完整 tPIE → 原始分数与稳健软压缩。

这里不默认加入三态差值、PBSA、熵项、多构象采样或按亲和力标签挑选姿态。
未知化学模型不会靠删除原子强行给分。补全结构也不代表保证任意结构都能成功。

## 文件分工

- `src/fmotpiefinal/pipeline/`：统一 CLI、输入、制备、QM 编排和失败恢复。
- `src/fmotpiefinal/core/`：结构解析、化学图、局部修复、最小化、分片、GAMESS 与收敛。
- `src/fmotpiefinal/scoring.py`：完整队列校验、tPIE 排序与软压缩。
- `tests/`：独立回归测试，不依赖历史实验目录。
- `docs/`：原理、安装、架构、验证范围、指标表和发布检查表。
- `examples/`：输入格式示例和参数化启动脚本；不包含第三方晶体结构。
- `tools/`：完整性检查和 wheel/sdist 构建。
- `SOURCE_PROVENANCE.json`：来源版本、模块重命名及文件哈希。

历史实验、旧 HF/三态/PBSA 分支、数据集、QM 日志、环境和二进制均不打包。
原研究目录保留，不因整理发布包而删除。

首次下载建议使用英文 README 的 sparse checkout 命令，仅拉取当前工具文件，
避免把历史大体积引擎和数据也展开到本地。需要追溯时再查看 `archive/`。

## 使用与验证

完整安装与启动命令见 `README.md`，科学流程见 `docs/WORKFLOW.md`。
完整计算需要独立安装的 GAMESS、DFTB 参数和 OpenMM/OpenFF 制备环境；
不是只用 pip 安装基础包就能完成 QM。

底层计算来源版本已验证 909/909 个结构；该结果与发布候选版自身的测试分开记录。
各靶点并非都预测良好，详见 `docs/VALIDATION.md`。不把开发集称作全新盲测。

**公开上传前仍需你确定源码许可证、作者/仓库信息和第三方使用权限。**
当前没有替你选择许可证。此次集成会提交到现有仓库，但不改变仓库可见性。
