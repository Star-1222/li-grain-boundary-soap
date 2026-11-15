# 本地工作区与GitHub的关系

GitHub只发布代码、报告、必要证据和数据索引；完整轨迹、矩阵、临时构建与整理工作树保留本地。

```text
工作区/
  README.md                 本地总入口
  repository/               指向实际GitHub检出的目录链接
  runs/                     完整运行记录与续算数据
  output/
    midterm/current/        当前可用汇报材料
    midterm/history/        旧汇报版本
    ni_soap_*/              完整SOAP结果与矩阵
  archive/
    source-snapshots/       原工作脚本快照
    worktrees/              整理用Git工作树
    organization/           本地移动记录与哈希
  tmp/                      保留原路径的构建依赖、OVITO环境
  .cache/                   图表中间包与新缓存
```

日常代码从 `repository/`进入；它是目录链接，不是复制出来的另一份代码。工作区根目录不是GitHub提交目标，具体位置以本地README和AGENTS为准。

原始runs与完整SOAP输出保留原路径，避免破坏续算和源文件引用。tmp中的运行环境也保留原路径，并从日常视图隐藏；本次没有删除原始轨迹、检查点或矩阵。

本次本地27组移动共1400个文件逐一核对SHA256，移动清单位于本地 `archive/organization/2026-11-15_moves.json`。Git工作树使用Git命令移动，保留分支与注册信息。旧绝对路径仅作为来源记录，当前研究使用仓库中的新入口。

当前中期PPT在 `output/midterm/current/`；Word仍为待核对稿，保留其原始内容。GitHub中的[研究状态](RESEARCH_STATUS.md)是动态进度，汇报文件按原生成日期查阅。

新试验从仓库根目录准备到 `runs/<新名称>`，或从本地工作区调用：

```bash
python repository/soap/current/prepare.py --dest runs/ni-soap-radial-001
```

该操作只准备脚本与输入，不启动计算。已有结果不覆写；新汇报构建使用tmp，图表中间包使用.cache。
