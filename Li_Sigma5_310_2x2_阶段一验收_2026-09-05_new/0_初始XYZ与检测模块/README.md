# 初始 XYZ 与检测模块

本文件夹同时存放初始构型的 XYZ 文件和结构检测脚本：

- `Sigma5_310_K6_J6_L20.xyz`：扩展 XYZ，含晶格与 grain_id；
- `Sigma5_310_K6_J6_L20_grain_id.xyz`：带 grain_id 的初始 XYZ；
- `cna_phase1.py`：阶段一 CNA-lite 结构检测脚本。

说明：`cna_phase1.py` 当前读取 LAMMPS `.data` 快照文件；初始 XYZ 文件
作为同一交付包中的初始构型基线保留。
