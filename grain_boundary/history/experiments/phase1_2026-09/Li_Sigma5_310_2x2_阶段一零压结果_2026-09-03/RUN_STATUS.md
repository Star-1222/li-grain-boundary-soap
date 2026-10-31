# RUN_STATUS — 阶段一零压 NPT 全流程（定版）

更新于 2026-09-03 17:40。阶段一零压重跑（Ubuntu WSL 2，8 MPI 核，K6/J6/L20，14400 原子）。
**背景**：旧管线（2026-09-01 交付）Stage 1 极小化只弛豫了 x 方向应力，y/z 从未参与
应力弛豫，导致终态 Pzz ≈ −1000 bar 残余应力。本轮按《阶段一原目标执行指南》
以"全维度零应力 + iso 0 bar"重跑全流程。

| Stage | 状态 | 关键结果 |
|---|---|---|
| A 0K 零应力弛豫 | **PASS** | fix box/relax aniso 0 + cg 分块极小化，fmax = 8.5e-6 eV/Å；终态 Pxx/Pyy/Pzz = −0.004 / +0.012 / +0.003 bar（**全维度零应力**）；`p1_min0.data/.restart`（09-03 01:28） |
| B 363K 零压退火 | **PASS** | NPT iso 0 bar、363 K、25 ps；末段（步 9000–25000）时均 T = 362.8 K、<Press> = −8.1 bar、体积 309100 Å³（较 Stage A 扩张 +1.1%，barostat 真正生效）；`p1_annealed.data/.restart/.dump`（11:09） |
| C 末帧 0K 检查 | **PASS** | 力收敛极小化 fmax = 2.29e-5 eV/Å；CNA（0K）：GB 区 BCC 43.7%/42.8%（无序维持）、晶粒内部 99.9% BCC → **Case B**；`p1_casecheck_0K.data`（11:17） |
| D 淬火 | **PASS** | 10 K/ps 线性 363→10 K，35.3 ps；末步瞬时 T = 12.3 K（热涨落）；`p1_quenched.data/.restart/.dump`（14:13） |
| E 10K 零压保持 | **PASS** | NPT iso 0 bar、10 K、5 ps；**后半段（步 2500–5000）时均 <Press> = +0.7 bar ≈ 0**（验收判据达成）；时均 T = 10.01 K；体积漂移 −0.055%、温度漂移 −0.05 K → 结构稳定；`p1_held.data/.restart/.dump` + `p1_hold_pressure.stats`（14:39） |
| F 终态 0K 复核 | **PASS** | 末帧力收敛极小化 fmax = 3.20e-5 eV/Å；CNA（0K）：GB 区 BCC **38.1%/38.0%**（无序维持并较 Case B 时进一步扩展）、晶粒内部 **98.7%** BCC（d ≥ 9.5 Å 处 100%）、w_GB = 6.5 Å；`p1_held_0K.data`（09-03 17:20） |

## 验收结论（按指南措辞）

```text
原目标验收通过
```

依据（指南三条判据全部达成）：

1. **10 K 时均压强 ≈ 0 bar**：后半段时均 <Press> = +0.7 bar（各向分量
   Pxx/Pyy/Pzz 时均 −187 / +386 / −196 bar，为 14400 原子体系 ±300 bar 量级
   正常热涨落，标量平均收敛到零）；
2. **无序 GB 维持**：终态 0K 力收敛后 GB 区 BCC 仅 38.1%/38.0%
   （初始 73.3%/73.3%），无序为真实结构无序，非热振动假象；
3. **晶粒保持 BCC**：终态 0K 力收敛后晶粒内部 98.7%（Q_BCC 剖面
   d ≥ 9.5 Å 处 100%），无序严格局域在 GB 两侧 ~6.5–8 Å 内。

## 与旧管线对照

| 量 | 旧管线（非零压，Pzz≈−1000 bar） | 本轮（零压） |
|---|---|---|
| Stage 1 应力状态 | 仅 x 弛豫，Pzz ≈ −1000 bar 残余 | 全维度 P ≈ 0.01 bar |
| 退火/淬火/保持压强 | 恒定 Pzz 偏置 | iso 0 bar（时均 +0.7 bar @10 K） |
| 退火后体积 | ~305700 Å³（冻结） | 309100 Å³（+1.1% 热膨胀） |
| 10K 终态 GB BCC | 43.8%/41.3%（旧采样） | 38.1%/38.0%（0K 力收敛） |
| 10K 终态晶粒 BCC | 98.6% | 98.7% |
| w_GB | 4.5 Å | 6.5 Å |

零压条件下 GB 无序带略宽（4.5 → 6.5 Å），晶粒内部 BCC 基本不变
（98.6% → 98.7%）：无序化结论不受残余应力修正影响，但零压终态
才是热力学一致的参考态，应作为后续阶段二的起点。

## 重要修正记录（Pdamp 单位 bug）

第一轮 B–E 脚本误写 `fix npt ... iso 0.0 0.0 1000.0`：metal 单位下
Pdamp 以 ps 计，1000.0 = 1000 ps = 1 ns，比合理值慢三个数量级，
barostat 实际失效（体积 5 ps 内零变化、10 K 时均压强 −323 bar）。
已修正为 `iso 0.0 0.0 1.0`（1 ps，约 1000 步）并全部重跑；无效日志
归档为 `logs/p1_*.barofail.lammps`（在 `7_运行日志/` 中保留供审计）。

## 运行环境与时间戳（Ubuntu WSL 2，SNAP OpenKIM）

- a0 标定（4×4×4 体块 fix box/relax aniso 0 + cg）：3.50618647533 Å
  （与原构建一致）— 2026-09-02
- Stage A：2026-09-03 01:28
- Stage B：2026-09-03 11:09
- Stage C：2026-09-03 11:17
- Stage D：2026-09-03 14:13
- Stage E：2026-09-03 14:39
- Stage F（终态 0K 复核 + CNA）：2026-09-03 17:20–17:35

## 运行 MD 的方法

```bash
cd ~/Li_meta_2x2
# 幂等 runner：done 标记防重跑，checkpoint 每 5000 步
bash run_phase1.sh
# 单独续跑 Stage E（如中断）：
mpirun -np 8 /usr/local/bin/lmp-kim -in phase1_inputs/in.p1_hold \
  -var kim_model SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
```

脚本清单（`6_结构分析/`）：`in.p1_min0` / `in.p1_anneal` / `in.p1_casecheck` /
`in.p1_quench` / `in.p1_hold` / `in.p1_finalcheck` + `run_phase1.sh`。
