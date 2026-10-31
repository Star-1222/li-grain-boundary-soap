# 现行晶界构造流程

唯一现行入口。模型参数见 [model.json](model.json)，完整协议见 [PROTOCOL.md](PROTOCOL.md)。

```text
构建并校验 14400 原子双晶
→ 0 K 初始弛豫
→ 363 K 退火 25 ps
→ 按 10 K/ps 淬火至 10 K
→ 10 K aniso 松弛
→ 检查应力、盒长和结构稳定性
→ 保存 standard_start.data / standard_start.restart
```

## 运行

要求 Python 3.11+、NumPy；运行 MD 还需要带 OpenKIM 模型的 LAMMPS 和 MPI。以下命令从仓库根目录执行。

```bash
# 只显示运行计划，不启动计算、不创建结果目录
python grain_boundary/current/run.py --out runs/li-standard-001 --dry-run

# 完整流程；目标目录必须不存在
python grain_boundary/current/run.py --out runs/li-standard-001 --lmp lmp --np 6

# 若已有符合协议的淬火检查点，可直接进入 aniso 松弛
python grain_boundary/current/run.py --from-quench /path/to/quenched.restart   --out runs/li-aniso-001 --lmp lmp --np 6
```

`--np 1` 时直接启动 LAMMPS，无需 MPI 启动器。`--steps` 指定本轮 aniso 步数，默认 5000 步（5 ps）；运行结束不等于稳定性验收通过。失败的目录保留诊断证据，重新运行使用新目录，不按旧 `.done` 自动跳过。

## 检查并保存

```bash
python grain_boundary/current/checkpoint.py inspect --run runs/li-standard-001
```

此命令生成 `review/thermo_summary.json`、候选文件哈希和默认不通过的 `review/review.json`。按[协议](PROTOCOL.md)检查应力、盒长、结构，在审查记录中填写标准、结果、证据文件和审查者。结构证据可放入本次运行目录下的 `review/`，不能只依据整体 BCC 比例作结论。

完成审查后保存：

```bash
python grain_boundary/current/checkpoint.py save --run runs/li-standard-001   --review runs/li-standard-001/review/review.json --dest runs/li-standard-start-001
```

未批准、缺检查项、证据缺失或候选文件在检查后改变，均拒绝保存。保存包包含结构、restart、运行源文件、日志、轨迹及审查证据，并生成 `STANDARD_START.json`。

仓库当前未发布标准起点，状态见 [STATUS.json](STATUS.json)。历史结果从[历史入口](../history/README.md)查阅。
