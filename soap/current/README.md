# 当前SOAP工作入口

当前任务是解决径向表示的数值可靠性与收敛问题。现有参数扫描和补充诊断已经完成，不代表已经获得最终参数。进度统一看[研究状态](../../maintenance/docs/RESEARCH_STATUS.md)。

## 准备独立研究副本

从仓库根目录执行：

```bash
python soap/current/prepare.py --dest runs/ni-soap-radial-001
```

命令只复制经输入哈希核对的Ni结构、原始脚本和运行说明，不启动科学计算，也不复制旧结果。目标目录必须不存在。

查看新目录中的 `RUN_ORDER.md`，在独立环境中逐步运行。原始科学脚本保持不变，通过保留它们的嵌套工作区结构解决路径与导入依赖；输出仅写入新工作区。`requirements-recorded.txt`记录原实验环境版本，不宣称已完成新的环境适配。

这是后续分层审计的复现基线，采用20000个环境对，不能直接与中期PPT的3600对通过率混合。下一次修改方法时，从新的独立副本开始，登记变动后再计算。

## 查结果

- [中期PPT统计口径](../experiments/2026-11-01_midterm/REPORT.md)
- [分层审计及扩阶补查](../experiments/2026-11-01_audit/REPORT.md)
- [径向基函数与受控结构对照](../experiments/2026-11-01_radial_followup/REPORT.md)

早期六阶段分析位于[历史归档](../history/README.md)，不是当前方法已经完成的六个必做步骤。
