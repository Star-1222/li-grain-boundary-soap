# TEAM_REVIEW - Li_Sigma5_310_2x2_run_2026-08-30.zip (Bilingual / 双语版)

本文件为双语版：中文版在前，英文版在后。This file contains the Chinese version first, followed by the English version.

# 中文版

# TEAM_REVIEW — Li_Sigma5_310_2x2_运行结果_2026-08-30.zip

审查日期：2026-09-01

审查对象：`Li_Sigma5_310_2x2_运行结果_2026-08-30.zip`

## 总评 / Verdict

```text
CONDITIONAL — MD 主链真实完成且协议大体合规，但结构分析存在区域定义缺陷，
当前不能作为已核实的 SUCCESS 结论提交
```

这份交付物确实包含真实运行的 Stage 0–4 结果：日志完整、步数/温度/盒子
对得上、淬火两次断点续跑斜率正确、10 K 保持 5 ps 完成。方向上支持
“候选亚稳无序晶界”，但 CNA 分析把两条等价 GB 的区域取成了不同宽度，
并且“晶粒内部 BCC 完好”的证据偏弱，需修正后再定论。

## 问题清单 / Findings

### P1 — CNA 的 GB 区域定义不对称，两条等价 GB 的比较无效

同一份初始结构里，GB1 区域取了约 3600 个原子，GB2 区域只取约 1920 个
原子（对应区域宽度约 16.6 Å vs 8.9 Å）。周期双晶的两条 GB 是严格等价
的，因此报告里初始 GB 的 73.3% vs 75.0%、mincheck 的 46.1% vs 59.2%、
终态 43.8% vs 48.0% 的差异主要来自采样宽度不同，不是物理差异。

影响：GB 区 BCC 分数的绝对值不可靠，至少其中一条 GB 的区域要么混入了
晶粒内部、要么没覆盖完整晶界无序带。

修正：对两条 GB 使用完全相同的、周期一致的对称区间（例如各自
`|x| < w` 和 `|x - Lx/2| < w`，或统一使用沿 x 的等宽 slab），重新出
CNA 报告。

### P1 — “晶粒内部保持 BCC”证据不足，且“完好”表述过强

mincheck（0 K 极小化）后晶粒内部 BCC 为 86.6%，10 K 终态为 86.5%，
也就是仍有约 13.5% 的内部原子不是 BCC。这算“大体保持 BCC”可以，
但叫“完好”是过强表述。

另外包内没有 `Q_BCC(x)` 剖面、没有 GB 宽度估计、没有局部能量剖面，
无法证明无序是局域在晶界、没有泄漏到晶粒内部，也无法给出 `w_GB`。
建议补 Stage 5 的沿 x 结构/能量剖面。

### P2 — mincheck 极小化未达到力收敛

mincheck 日志显示停止原因是能量容差，最终最大力分量仍有
`0.0342 eV/A`。用来去掉 363 K 热振动是有效的，但这不是一个严格收敛的
0 K 局域极小，不能当作“决定性验证”。建议以 `fmax < 1e-4 eV/A` 重跑，
或至少在报告中写明剩余力。

### P2 — CNA 分析不可复现

包内只有 `cna_*.txt` 结果，没有分析脚本，也没有 `initial.data`、
`bulk363.data`、`mincheck.data` 等输入/中间文件。第三方无法复核这些
BCC 百分比。建议随包提供脚本与所用快照。

### P2 — Stage 2 的 363 K 有效停留时间约 1.1 ps，不是 2 ps

stage2 日志显示：step 0 设置 363 K 后，0.1 ps 降到约 183 K，到约
0.9 ps 才回到 355 K，随后才在 363 K 附近波动。即 2 ps 里约一半时间在
升温，真正处于 363 K 的停留约 1.1 ps。若严格执行“363 K 保持 2 ps”，
应先做 363 K 预平衡，或把 run 加长到约 3 ps。

### P2 — 淬火 dump 拼接规则不唯一

`stage3_quench_resume.dump` 含 30000–33000 帧，`resume2.dump` 又从
checkpoint.30000 重跑 30000–35300。直接“三段拼接”会让 30000–33000
重复且可能互相冲突。最终构型来自 resume2，正确，但应明确拼接以
resume2 为准，或裁剪废弃段。

### P2 — 缺复现输入与溯源

包内没有 LAMMPS 输入脚本、PARAMETERS/PROVENANCE。日志可确认 LAMMPS
为 10 Dec 2025，但 OpenKIM 版本、随机种子、脚本版本未记录。建议补上
一份 `RUN_META.md`。

## 已验证正确的部分 / What checks out

- Stage 0：14400 atoms，K6/J6/L20，`STAGE0_CHECK_OK`；
- Stage 1：cg 极小化 `fmax=2.49e-5 < 1e-4`，x 弛豫 66.5252→66.0292 Å，
  y/z 保持 CSL 周期；
- Stage 2：`fix npt temp 363 363 x 0 0`，2000 步（2 ps），dump 覆盖
  step 0–2000，终态 T≈361 K；
- Stage 3：总 35300 步（35.3 ps）；15000 步处续跑温度取 213 K、30000
  步处取 63 K，剩余步数 20300/5300，均与 10 K/ps 线性斜坡一致；最终
  step 35300 T≈12.3 K，`STAGE3_QUENCH_DONE`；
- Stage 4：10 K NPT 5000 步；T=9.91–10.24 K；总能量 -27243.5→
  -27247.6 eV（约 0.00028 eV/atom 漂移）；体积 308008→307987 A3
  （-0.007%）；`STAGE4_10K_HOLD_DONE`；
- mincheck 思路正确：对 Stage 2 终态做 0 K 弛豫以区分热振动与真实结构
  无序，并使用同一 SNAP 势。

## 建议 / Recommended action

1. 修正 CNA 区域定义（等价 GB 等宽采样）并重算；
2. 补 `Q_BCC(x)`、`E_local(x)` 与 `w_GB`，量化无序局域性；
3. mincheck 以力收敛重跑并报告剩余力；
4. 提供 CNA 脚本和所有中间文件；
5. 淬火拼接规则写明以 resume2 为准；
6. 补复现元数据（LAMMPS/OpenKIM 版本、种子、脚本版本）；
7. 修正后再按决策规则给出最终结论。

## 一句话总结

MD 是真跑完的，方向也支持候选亚稳无序晶界；但 CNA 区域 bug 和缺失的
结构剖面让“SUCCESS”目前只能算初步候选，不能作为已核实结论。


# English Version

# TEAM_REVIEW - Li_Sigma5_310_2x2_run_2026-08-30.zip

Review date: 2026-09-01

Reviewed item: `Li_Sigma5_310_2x2_运行结果_2026-08-30.zip`

## Verdict

```text
CONDITIONAL - the MD main chain was genuinely completed and the protocol is
mostly compliant, but the structural analysis has a region-definition defect,
so it cannot yet be submitted as a verified SUCCESS conclusion
```

This deliverable does contain real Stage 0-4 results: complete logs, matching
steps/temperature/box, correct slopes for the two checkpoint resumes during
quenching, and a completed 5 ps 10 K hold. The direction supports a
"candidate metastable disordered GB", but the CNA analysis assigned different
widths to the two equivalent GB regions and the evidence for "intact BCC
inside grains" is weak. The conclusion should be finalized after correction.

## Findings

### P1 - CNA GB region definition is asymmetric, so comparison of the two
equivalent GBs is invalid

In the same initial structure, GB1 used about 3600 atoms while GB2 used about
1920 atoms (region widths around 16.6 A versus 8.9 A). The two GBs in a
periodic bicrystal are strictly equivalent, so differences such as the
reported initial 73.3% vs 75.0%, mincheck 46.1% vs 59.2%, and final 43.8% vs
48.0% come mainly from different sampling widths, not physical differences.

Impact: the absolute BCC fractions in the GB region are unreliable; at least
one of the two GB regions either mixes in grain interior or does not cover
the full disordered GB band.

Fix: use exactly the same, periodically consistent symmetric interval for
both GBs (for example `|x| < w` and `|x - Lx/2| < w`, or uniform equal-width
slabs along x) and regenerate the CNA report.

### P1 - Evidence for "grains remain BCC" is insufficient, and "intact" is
too strong

After mincheck (0 K minimization), the interior BCC fraction is 86.6%; at the
10 K final state it is 86.5%, meaning about 13.5% of interior atoms are not
BCC. Calling this "broadly BCC" is acceptable, but "intact" overstates it.

The package also has no `Q_BCC(x)` profile, no GB-width estimate, and no local
energy profile, so it cannot prove disorder is localized at the GB or has not
leaked into the grains, nor can it provide `w_GB`. Stage 5 structural/energy
profiles along x should be added.

### P2 - mincheck minimization did not reach force convergence

The mincheck log shows the run stopped on energy tolerance, with a final
maximum force component of `0.0342 eV/A`. This is useful for removing 363 K
thermal vibration but is not a strictly converged 0 K local minimum and cannot
be treated as decisive validation. Rerun with `fmax < 1e-4 eV/A`, or at least
report the residual force.

### P2 - CNA analysis is not reproducible

The package only contains `cna_*.txt` results, with no analysis script and no
`initial.data`, `bulk363.data`, `mincheck.data`, or other input/intermediate
files. Third parties cannot recheck these BCC percentages. Scripts and the
used snapshots should be provided with the package.

### P2 - Effective residence time at 363 K in Stage 2 is about 1.1 ps, not
2 ps

The stage2 log shows that after setting 363 K at step 0, the temperature drops
to about 183 K at 0.1 ps and only returns to 355 K around 0.9 ps, after which
it fluctuates near 363 K. About half of the 2 ps is spent heating up, so the
actual 363 K residence is about 1.1 ps. If "hold at 363 K for 2 ps" is to be
enforced strictly, pre-equilibrate at 363 K or extend the run to about 3 ps.

### P2 - Quench dump concatenation rule is not unique

`stage3_quench_resume.dump` contains frames 30000-33000, and `resume2.dump`
reruns 30000-35300 from checkpoint.30000. Simply concatenating three segments
duplicates 30000-33000 and can create conflicts. The final configuration
comes from resume2 and is correct, but the concatenation rule should explicitly
state resume2 is authoritative or discard the obsolete segment.

### P2 - Reproduction inputs and provenance are missing

The package has no LAMMPS input scripts, PARAMETERS/PROVENANCE. Logs confirm
LAMMPS is 10 Dec 2025, but OpenKIM version, random seed, and script version
are not recorded. Add a `RUN_META.md`.

## What checks out

- Stage 0: 14400 atoms, K6/J6/L20, `STAGE0_CHECK_OK`;
- Stage 1: cg minimization `fmax=2.49e-5 < 1e-4`; x relaxes 66.5252 to
  66.0292 A; y/z keep the CSL period;
- Stage 2: `fix npt temp 363 363 x 0 0`, 2000 steps (2 ps), dump covers steps
  0-2000, final T around 361 K;
- Stage 3: total 35300 steps (35.3 ps); resume at step 15000 starts at 213 K,
  at step 30000 starts at 63 K, with 20300/5300 remaining steps, both matching
  the 10 K/ps linear ramp; final step 35300 T around 12.3 K,
  `STAGE3_QUENCH_DONE`;
- Stage 4: 10 K NPT 5000 steps; T=9.91-10.24 K; total energy
  -27243.5 to -27247.6 eV (about 0.00028 eV/atom drift); volume
  308008 to 307987 A3 (-0.007%); `STAGE4_10K_HOLD_DONE`;
- mincheck concept is correct: 0 K relaxation of the Stage 2 final state to
  separate thermal vibration from true structural disorder, using the same
  SNAP potential.

## Recommended action

1. Fix the CNA region definition (equal-width sampling for equivalent GBs)
   and recompute;
2. add `Q_BCC(x)`, `E_local(x)`, and `w_GB` to quantify disorder
   localization;
3. rerun mincheck with force convergence and report residual force;
4. provide CNA scripts and all intermediate files;
5. state that resume2 is authoritative for quench concatenation;
6. add reproduction metadata (LAMMPS/OpenKIM versions, seed, script version);
7. finalize the conclusion under the decision rule after correction.

## One-sentence summary

The MD chain was genuinely completed and the direction supports a candidate
metastable disordered GB, but the CNA region bug and missing structural
profiles mean "SUCCESS" can only be a preliminary candidate for now, not a
verified conclusion.

