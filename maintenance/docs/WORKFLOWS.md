# 工作流入口

## 晶界构造

唯一现行入口是 [grain_boundary/current](../../grain_boundary/current/README.md)。标准链条为：

**淬火 → 10 K 下 aniso 松弛 → 检查应力、盒长和结构稳定性 → 保存标准起点。**

完整构建、初始弛豫、退火及淬火参数见[标准协议](../../grain_boundary/current/PROTOCOL.md)。默认模型为 14400 原子。已有淬火 restart 可直接进入 aniso 阶段，但要保留输入来源；结果先称为 candidate，通过明确审查后才保存为 standard_start。

旧的 3600 原子基准、早期 2×2 脚手架、阶段一 iso 实验和控压对照统一收在 [history](../../grain_boundary/history/README.md)，不作为新运行的默认流程。它们是不同实验记录，不能拼成一次实际运行。

## SOAP 方法

SOAP 方法线仍从 [soap/](../../soap/README.md) 进入，独立于 Li 的阶段编号。

```bash
python maintenance/tools/prepare_soap_workspace.py --dest runs/ni-soap-001
```

准备工具把 Stage01–06 的脚本和 CFG 复制到独立目录，生成 `RUN_ORDER.md`，按其中顺序执行。目标目录必须不存在；不复制历史输出、不启动计算。依赖列表位于 `soap/06_环境与实现/requirements_stages.txt`，尚未冻结包版本。

跨阶段的数据、报告和逐文件索引见 SOAP 项目入口。缺失轨迹及矩阵见[数据恢复清单](../records/DATA_MANIFEST.md)。
