# Li Σ5 [001](310) 晶界 —— 残余应力成因 × 控压方案对照实验

归档日期：2026-10-15　｜　体系：BCC Li Σ5 [001](310) 对称倾转晶界，14400 原子
起点构型：`stage3b_quench/stage3b_quenched.restart`（NPT 363 K 退火 → 淬火后的同一构型）

## 一句话结论

在同一个淬火构型上、相同时长（5 ps）内只改变控压方式，各向异性控压（aniso）
把残余偏应力从 **627.6 ± 13.9 bar 压到 15.9 ± 6.3 bar**（约 40 倍），
代价是把 [001] 倾转轴方向压缩 ε_zz = −0.99%，CSL 周期在该轴上被破坏。

同时发现原协议存在一个 bug：`units metal` 下 `fix npt ... 1000.0` 的阻尼单位是 ps，
即 Pdamp = 1000 ps，而保持阶段只有 5 ps —— 活塞重到几乎不动，
旧流程测到的残余应力把"耦合方式"和"活塞过重"两个因素混在了一起。

## 目录结构

```
00_双周报告/      双周报 docx / pdf（已回填真实数值）+ 第一部分独立版
01_数据产物/      结果报告.md（完整数据表）、summary.json / .md、
                 三份结构分析 JSON、固定阈值晶界宽度、两分支热力学 CSV
02_图/            fig_stress.png  应力分量与偏应力随时间演化
                  fig_strain.png  三个方向的应变演化
                  fig_profile.png BCC 分数 / q6 沿 GB 法向剖面
03_原始日志/      两分支完整 LAMMPS 日志（baro_xonly_fast / baro_aniso_fast）
04_脚本/          LAMMPS 输入 in.baro_hold、启动 run_barostat.sh、
                 后处理与成图全部 Python 脚本 + Steinhardt 校验脚本
```

## 关键数值（末段 2.4 ps 平均）

| 指标 | x-only（仅 GB 法向控压） | aniso（三向独立活塞） |
|---|---|---|
| P_xx (bar) | +1.3 | −9.6 |
| P_yy (bar) | −708.8 | −11.8 |
| P_zz (bar) | −1249.8 | −35.2 |
| **偏应力** | **627.6 ± 13.9** | **15.9 ± 6.3** |
| ε_xx | −0.19 % | +0.20 % |
| ε_yy | 0.000 %（锁定） | −0.25 % |
| ε_zz | 0.000 %（锁定） | −0.99 % |
| CN==8 原子占比 | 80.6 % | 78.1 % |
| 角度判据 BCC 占比 | 90.3 % | 89.9 % |
| w_GB（Q<0.5 固定阈值） | 11.0 Å | 13.2 Å |

起点（淬火态）的 CN==8 / 角度判据分别为 80.6 % / 90.3 %，与 x-only 分支一致。

## 复现步骤

环境：WSL-Ubuntu + LAMMPS 12（`lmp-kim`，/usr/local/bin），KIM 模型
`SNAP_ZuoChenLi_2019_Li__MO_732106099012_000`，24 逻辑核（约 12 物理核）。

```bash
cd /home/lenovo/Li_meta_2x2
cp 04_脚本/in.baro_hold   exp_barostat/lammps/
cp 04_脚本/run_barostat.sh .

# 1) MD：两分支各跑 5000 步（5 ps），Pdamp = 0.1 ps
NSTEPS=5000 NCPU=6 TAGS=xonly_fast,aniso_fast bash run_barostat.sh

# 2) 热力学对照表
python3 exp_barostat/compare_barostat.py --tags xonly_fast,aniso_fast \
        --outdir exp_barostat/analysis

# 3) 结构分析（三份结构统一口径，rcut 取 g(r) 第一谷 3.25 Å）
python3 exp_barostat/analyze_structure.py stage3b_quench/stage3b_quenched.data \
        --name start_quenched    --rcut 3.25 --outdir exp_barostat/analysis
python3 exp_barostat/analyze_structure.py exp_barostat/xonly_fast/final.data \
        --name final_xonly_fast  --rcut 3.25 --outdir exp_barostat/analysis
python3 exp_barostat/analyze_structure.py exp_barostat/aniso_fast/final.data \
        --name final_aniso_fast  --rcut 3.25 --outdir exp_barostat/analysis

# 4) 固定阈值晶界宽度 + 成图
python3 exp_barostat/gb_width_fixed.py
python3 exp_barostat/make_figures.py
```

原始轨迹与中间结果保留在 WSL：`/home/lenovo/Li_meta_2x2/exp_barostat/`
（含 `slow_run_evidence/` 下 Pdamp=1000 ps 那组作为 bug 证据）。

## 方法口径说明（重要）

- **BCC 判据**：项目 generator/verify.py 里的"8 近邻 + 6 次近邻 ±2%"准则只适用于
  0 K 完美晶格，在 MD 热快照上给出全 0。本轮改用
  ① 配位数（CN==8，rcut = g(r) 第一谷 3.25 Å）② 键角/Steinhardt q4–q6 两个判据。
- **Steinhardt 实现已校验**：scipy 新版 `sph_harm_y` 的角参数顺序是 `(polar, azimuthal)`，
  与旧接口相反。修正后完美 FCC 给出 q4=0.1909 / q6=0.5745，
  BCC(14 近邻) 给出 q4=0.0364 / q6=0.5107，与文献一致。见 `04_脚本/validate_q.py`。
- **体系大小这一半变量本轮刻意固定住了**（2×2 超胞 14400 原子），
  要单独回答它还需要跑 3×3（K6 / J9 / L30）。

## 未闭合的问题

1. aniso 引入的 ε_zz = −0.99% 宏观应变是否影响后续 M–τc 映射，尚未评估。
2. 单条轨迹，未做多种子统计；协议中写明"单条成功轨迹不是系综结论"。
3. 建议在 Stage 5 之前只补一段 2 ps 的短 aniso NPT，
   把"应力为零"与"CSL 严格保持"拆成两版结构分别归档比较。
