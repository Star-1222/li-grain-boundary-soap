# RUN_META — 阶段一零压复现元数据

运行日期：2026-09-02 ~ 09-03。阶段一零压 NPT 全流程重跑（定版）。
上游：2026-09-01 交付的旧管线（`Li_Sigma5_310_2x2_最终结果_2026-09-01/`），
本目录为其余应力问题的修正重跑，复用同一初始构型（Stage 0 构建校验 PASS）。

## 软件

| 组件 | 版本 |
|---|---|
| LAMMPS | 10 Dec 2025（`/usr/local/bin/lmp-kim`，含 KIM 包） |
| KIM API | libkim-api 2.4.1-3 |
| 势函数 | SNAP_ZuoChenLi_2019_Li__MO_732106099012_000（ZuoChenLi 2019 SNAP） |
| MPI | OpenMPI，8 ranks |
| 分析脚本运行时 | WSL Ubuntu python3 + numpy + scipy |

## 关键标定与输入

| 项 | 值 |
|---|---|
| a0 标定（4×4×4 体块 fix box/relax aniso 0 + cg） | 3.50618647533 Å（与原构建一致） |
| 初始构型 | 旧管线 Stage 0 输出（14400 原子，K6/J6/L20，2×2） |
| timestep | 2 fs（metal 单位） |

## 运行参数与随机性

| Stage | 脚本 | 关键变量 / 随机性 |
|---|---|---|
| A 0K 零应力弛豫 | `in.p1_min0` | fix box/relax aniso 0 + cg 分块（确定性），fmax = 8.5e-6 eV/Å |
| B 363K 零压退火 | `in.p1_anneal` | NPT iso 0 bar（**Pdamp = 1.0 ps**），Tdamp 0.1 ps，25 ps，velocity seed = 12345 |
| C 末帧 0K 检查 | `in.p1_casecheck` | 定框 cg 分块极小化，fmax = 2.29e-5 eV/Å（诊断用，不进生产轨迹） |
| D 淬火 | `in.p1_quench` | NPT iso 0 bar，363→10 K @ 10 K/ps（0.02 K/步），35.3 ps |
| E 10K 零压保持 | `in.p1_hold` | NPT iso 0 bar，10 K，5 ps；fix ave/time 输出 `p1_hold_pressure.stats` |
| F 终态 0K 复核 | `in.p1_finalcheck` | 定框 cg 分块极小化，fmax = 3.20e-5 eV/Å（诊断用） |

## 脚本版本

| 脚本 | 版本 | 说明 |
|---|---|---|
| `run_phase1.sh` | 2026-09-03 | 幂等 runner：done 标记防重跑，checkpoint 每 5000 步，串行驱动 A→E |
| `cna_phase1.py` | 基于 cna_v2.py（2026-09-01） | 用已验证的 cna_v2.py 替换文件列表生成；对称等宽 GB 采样（|d_GB|<8 Å，两条 GB 各 3600 原子）；输入 initial / p1_min0 / casecheck / p1_final（0K 力收敛终态） |

## 中断与修正记录（影响复现）

| 时间 | 事件 | 恢复 |
|---|---|---|
| 09-02 深夜 | Stage E 首版脚本 `v_pxx` 未定义崩溃 | 写 `in.p1_hold` 修复版重跑 |
| 09-03 上午 | **Pdamp 单位 bug**：`iso 0.0 0.0 1000.0`（= 1000 ps）致 barostat 失效（10K 时均 −323 bar、体积冻结） | 三个脚本修正为 `iso 0.0 0.0 1.0`；无效日志归档 `.barofail.lammps`；清理 B–E 输出后全量重跑（Stage A 保留） |
| 09-03 14:39 | Stage E 完成 | 无需续跑 |
| 09-03 17:20 | Stage F 终态 0K 复核（追加，1 分 46 秒） | — |

重跑（Pdamp 修正后）B→E 一次通过，无中断、无拼接。

## 复现命令

```bash
cd ~/Li_meta_2x2
bash run_phase1.sh            # 幂等：已完成阶段自动跳过
# 终态 0K 复核：
mpirun -np 8 /usr/local/bin/lmp-kim -in phase1_inputs/in.p1_finalcheck \
  -var kim_model SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
# CNA（在 ~/cna_p1，需 initial.data / p1_min0_sel.data /
# p1_casecheck_0K.data / p1_held_0K.data）：
python3 cna_phase1.py
```

注意：LAMMPS restart 不保存 NPT 随机数流状态；若从 checkpoint 续跑，
轨迹与连续运行有微小随机性差异（对本流程的结构结论无影响）。
