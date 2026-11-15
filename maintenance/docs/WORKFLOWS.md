# 工作流入口

## 晶界构造

老师要求的高温制备、淬火与低温保持原目标已完成，见[完成记录](../../grain_boundary/current/TEACHER_PLAN_2026-11-05.md)。下文运行入口是2026-10-31版复现与aniso优化代码。历史实跑高温保持25 ps、低温保持5 ps，满足原指导的时间要求。10 K、零压稳定超过1 ps是继续研究判据；相变问题暂不作为本阶段前置研究。

唯一现行入口是 [grain_boundary/current](../../grain_boundary/current/README.md)。标准链条为：

**淬火 → 10 K 下 aniso 松弛 → 检查应力、盒长和结构稳定性 → 保存标准起点。**

完整构建、初始弛豫、退火及淬火参数见[标准协议](../../grain_boundary/current/PROTOCOL.md)。默认模型为 14400 原子。已有淬火 restart 可直接进入 aniso 阶段，但要保留输入来源；结果先称为 candidate，通过明确审查后才保存为 standard_start。

旧的 3600 原子基准、早期 2×2 脚手架、阶段一 iso 实验和控压对照统一收在 [history](../../grain_boundary/history/README.md)，不作为新运行的默认流程。它们是不同实验记录，不能拼成一次实际运行。

## SOAP方法

当前参数研究从[soap/current](../../soap/current/README.md)进入：

```bash
python soap/current/prepare.py --dest runs/ni-soap-radial-001
```

只准备经输入哈希核对的脚本与结构，按新目录的 `RUN_ORDER.md`逐步执行，不复制旧结果。准备的是20000对分层审计的复现基线，不是PPT的3600对统计。

已完成扫描与诊断见[实验索引](../../soap/experiments/README.md)。早期六阶段整体迁入[历史归档](../../soap/history/README.md)，确需复现时才使用 `maintenance/tools/prepare_soap_workspace.py`。

原始轨迹和大矩阵继续按[数据恢复清单](../records/DATA_MANIFEST.md)及实验发布清单追溯；当前结论统一看[研究状态](RESEARCH_STATUS.md)。
