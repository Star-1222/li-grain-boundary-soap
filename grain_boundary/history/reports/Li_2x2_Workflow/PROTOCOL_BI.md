# Protocol: Li Sigma5 [001](310) 2x2 Metastable GB (Bilingual / 双语版)

# 中文版

# 协议：Li Σ5 [001](310) 2×2 亚稳晶界（中文版）

## 0. 目的

从已验证的 Li Σ5 [001](310) 双晶出发：

```text
旧模型 3600 atoms (K=6, J=3, L=10)
   -> 2×2 面内扩展 (K=6, J=6, L=20, ~14400 atoms)
   -> Stage 1 初始弛豫
   -> Stage 2 363 K NPT / 2 ps 高温重构
   -> Stage 3 363 K -> 10 K，10 K/ps 淬火（35.3 ps）
   -> Stage 4 10 K / 0 bar 保持 > 1 ps（计划 5 ps）
   -> 结构与能量审计
```

宏观晶界身份必须是 Σ5 [001](310)，不允许变成其他 Σ 值、其他转轴或其他
宏观晶界平面。微观原子排列允许强烈重构。

## 1. 坐标系与“2×2”的物理含义

旧版权威生产项目的坐标约定（已由 v2 生成器与验证报告确认）：

```text
x = GB 法向，Lx = K * sqrt(10) * a
y = GB 面内方向 [130]，Ly = J * sqrt(10) * a
z = GB 面内方向 [001]，Lz = L * a
```

任务书要求“两个平行于 GB 的方向各放大 2 倍，GB 法向不变”，因此：

```text
K = 6 不变（x，GB normal）
J = 3 -> 6（y，in-plane ×2）
L = 10 -> 20（z，in-plane ×2）
```

注意：导师口头表述“模型可在 x 与 y 方向上加长”中的“x/y”若按本项目坐标
系理解会与旧版约定冲突。本实验以旧版权威约定为准：真正加倍的是 GB 面内
两个方向（y、z）。此差异已记录在 `PROVENANCE.md`。

## 2. Stage 0 - 生成与校验

命令：

```bash
python -m generator.generate \
  --a 3.506186475312 --k 6 --j 6 --l 20 \
  --outdir stage0_build_validate
```

已完成。结果：

```text
atoms       = 14400 (g1=7200, g2=7200)
box         = 66.525211 x 66.525211 x 70.123730 A
GBs         = x=0 与 x=33.262605 A，两条等价
min dist    = 2.217507 A
bulk BCC    = 10080/10080 sampled perfect
misorient.  = 36.86989765 deg about [001]
```

## 3. Stage 1 - 初始弛豫

沿用旧版验证过的极小化逻辑：`box/relax x 0.0`（GB 法向 x 零压弛豫，
y/z 固定在零压 CSL 周期），CG 极小化，分块 checkpoint。

输入：`stage0_build_validate/Sigma5_310_K6_J6_L20.data`
输出：`stage1_relax/stage1_minimized.{data,restart}`，log、trajectory、
thermodynamic history。

## 4. Stage 2 - 高温 NPT 重构（关键新实验）

```text
ensemble   = NPT
T          = 363 K（= 0.8 * 454 K）
P          = 0 bar（x 方向，GB normal；y/z 保持 CSL 周期）
duration   = 2 ps
timestep   = 0.001 ps
steps      = 2000
```

Li 熔点约 454 K，80% 约为 363 K。目的不是熔化：晶界应当强烈无序化/重构，
而晶粒内部仍可辨认 BCC。

### Stage 2 判定

必须在淬火前先检查完整轨迹（t=0, 0.5, 1, 1.5, 2 ps）：

- Case A：晶界基本不变 -> `Stage 2 = insufficient`，不进入 Stage 3；
- Case B：晶界明显无序/重构且晶粒保持 BCC -> `Stage 2 = PASS`，进入
  Stage 3；
- Case C：晶粒内部大面积失去 BCC -> `Stage 2 = FAIL`，不是想要的亚稳
  晶界。

## 5. Stage 3 - 快速淬火

```text
T_start     = 363 K
T_end       = 10 K
cooling     = 10 K/ps
duration    = (363-10)/10 = 35.3 ps
steps       = 35300
ensemble    = NPT（x 零压；y/z 保持 CSL 周期）
```

保存完整轨迹与 restart。淬火过程中晶界允许部分弛豫，关键是“重构/无序
特征是否保留到低温”。

注意：原子在 10 K 附近几乎不动，不代表亚稳。必须对比 Stage 1 初始 GB、
Stage 2 高温 GB、Stage 3 淬火 GB、Stage 4 最终 GB。

## 6. Stage 4 - 10 K / 0 bar 稳定性

```text
T          = 10 K
P          = 0 bar
minimum    = > 1 ps
planned    = 5 ps（10 ps 可选）
steps      = 5000（5 ps 计划）
ensemble   = NPT（x 零压；y/z 保持 CSL 周期）
```

### 最小成功判据

```text
1. GB 无序/重构保持；
2. 晶粒内部保持 BCC；
3. 无明显结构塌缩；
4. 10 K、0 bar 下存活超过 1 ps。
```

## 7. Stage 5 - 结构、能量与 SOAP 审计

只有 Stage 4 通过后才做：

- CNA / adaptive CNA / PTM 的晶粒内部 BCC 分数；
- 沿 x 的局部有序度剖面 `Q_BCC(x)`、局部能量剖面 `E_local(x)`；
- 有效晶界宽度 `w_GB`；
- 能量随时间的演化（仅作辅助证据）；
- 如需要，再做 SOAP 环境对比（初始 / 高温 / 淬火 / 10 K 终态 / bulk）。

逻辑顺序必须是 `MD stability -> structure -> energy -> SOAP`，不能用
SOAP 差异代替动力学稳定性。

## 8. Stage 6 - 尺寸效应

2×2 成功后再跑 3×3（K=6, J=9, L=30），重复同一协议并对比。第一次 2×2
成功只代表“存在一个候选路径”，不代表尺寸无关的物理结论。

## 9. 必须遵守的纪律

1. 不修改旧权威项目；
2. 不复用旧输出冒充新实验结果；
3. 不混用旧晶格常数（生产值 `a=3.506186475312 A`）；
4. 不把低温冻结当作亚稳证明；
5. 不用单张终态快照作证明，保留轨迹；
6. 不用总能作为唯一判据；
7. 不用 SOAP 替代动力学；
8. 暂不做相变分类；
9. 周期双晶有两条 GB，解释 GB 能时必须除以 2 个界面；
10. 单条成功轨迹不是系综结论。


# English Version

# Protocol: Li Sigma5 [001](310) 2x2 Metastable GB (English)

## 0. Objective

Starting from the validated Li Sigma5 [001](310) bicrystal:

```text
old model 3600 atoms (K=6, J=3, L=10)
   -> 2x2 in-plane extension (K=6, J=6, L=20, ~14400 atoms)
   -> Stage 1 initial relaxation
   -> Stage 2 363 K NPT / 2 ps high-temperature reconstruction
   -> Stage 3 363 K to 10 K at 10 K/ps quench (35.3 ps)
   -> Stage 4 10 K / 0 bar hold > 1 ps (planned 5 ps)
   -> structure and energy audit
```

The macroscopic GB identity must remain Sigma5 [001](310); it must not become
another Sigma value, another rotation axis, or another macroscopic GB plane.
Strong local atomic rearrangement is allowed.

## 1. Coordinate system and the physical meaning of "2x2"

Coordinate convention of the authoritative legacy production project
(confirmed by the v2 generator and verification reports):

```text
x = GB normal, Lx = K * sqrt(10) * a
y = in-plane [130], Ly = J * sqrt(10) * a
z = in-plane [001], Lz = L * a
```

The task requires "doubling both directions parallel to the GB while keeping
the GB normal unchanged", so:

```text
K = 6 unchanged (x, GB normal)
J = 3 -> 6 (y, in-plane x2)
L = 10 -> 20 (z, in-plane x2)
```

Note: the advisor's oral "elongate in x and y" would conflict with the legacy
convention if interpreted with this project's coordinate labels. This
experiment follows the authoritative legacy convention: the two in-plane GB
directions (y, z) are doubled. The difference is recorded in
`PROVENANCE.md`.

## 2. Stage 0 - generation and validation

Command:

```bash
python -m generator.generate \
  --a 3.506186475312 --k 6 --j 6 --l 20 \
  --outdir stage0_build_validate
```

Completed. Results:

```text
atoms       = 14400 (g1=7200, g2=7200)
box         = 66.525211 x 66.525211 x 70.123730 A
GBs         = x=0 and x=33.262605 A, two equivalent GBs
min dist    = 2.217507 A
bulk BCC    = 10080/10080 sampled perfect
misorient.  = 36.86989765 deg about [001]
```

## 3. Stage 1 - initial relaxation

Use the legacy-validated minimization logic: `box/relax x 0.0` (zero-pressure
relaxation of the GB-normal x, while y/z stay at the zero-pressure CSL
periods), CG minimization, chunked checkpoints.

Input: `stage0_build_validate/Sigma5_310_K6_J6_L20.data`
Output: `stage1_relax/stage1_minimized.{data,restart}`, log, trajectory, and
thermodynamic history.

## 4. Stage 2 - high-temperature NPT reconstruction (key new experiment)

```text
ensemble   = NPT
T          = 363 K (= 0.8 * 454 K)
P          = 0 bar (x direction, GB normal; y/z keep CSL periods)
duration   = 2 ps
timestep   = 0.001 ps
steps      = 2000
```

The Li melting point is about 454 K; 80% is about 363 K. The goal is not
melting: the GB should strongly disorder/reconstruct while grain interiors
remain recognizable BCC.

### Stage 2 decision

Inspect the full trajectory (t=0, 0.5, 1, 1.5, 2 ps) before quenching:

- Case A: GB essentially unchanged -> `Stage 2 = insufficient`, do not enter
  Stage 3;
- Case B: GB clearly disordered/reconstructed and grains remain BCC ->
  `Stage 2 = PASS`, enter Stage 3;
- Case C: grain interiors lose BCC over a large region -> `Stage 2 = FAIL`,
  not the desired metastable GB.

## 5. Stage 3 - fast quench

```text
T_start     = 363 K
T_end       = 10 K
cooling     = 10 K/ps
duration    = (363-10)/10 = 35.3 ps
steps       = 35300
ensemble    = NPT (x zero pressure; y/z keep CSL periods)
```

Save the full trajectory and restart. Some GB relaxation during quenching is
allowed; the key question is whether the reconstruction/disorder features
survive to low temperature.

Note: atoms barely moving near 10 K is not evidence of metastability. Compare
Stage 1 initial GB, Stage 2 high-temperature GB, Stage 3 quenched GB, and
Stage 4 final GB.

## 6. Stage 4 - 10 K / 0 bar stability

```text
T          = 10 K
P          = 0 bar
minimum    = > 1 ps
planned    = 5 ps (10 ps optional)
steps      = 5000 (5 ps plan)
ensemble   = NPT (x zero pressure; y/z keep CSL periods)
```

### Minimum success criterion

```text
1. GB disorder/reconstruction is retained;
2. grain interiors remain BCC;
3. no obvious structural collapse;
4. survives more than 1 ps at 10 K and 0 bar.
```

## 7. Stage 5 - structure, energy, and SOAP audit

Only after Stage 4 passes:

- CNA / adaptive CNA / PTM BCC fractions in grain interiors;
- local order profile `Q_BCC(x)` and local energy profile `E_local(x)` along
  x;
- effective GB width `w_GB`;
- energy evolution over time (supporting evidence only);
- if needed, SOAP environment comparison (initial / high-T / quenched /
  10 K final / bulk).

The logical order must be `MD stability -> structure -> energy -> SOAP`;
SOAP differences cannot replace dynamical stability.

## 8. Stage 6 - size effect

After the 2x2 succeeds, run 3x3 (K=6, J=9, L=30) with the same protocol and
compare. A first successful 2x2 run only shows "one candidate path exists",
not a size-independent physical conclusion.

## 9. Mandatory discipline

1. Do not modify the authoritative legacy project;
2. do not reuse old outputs as new experimental results;
3. do not mix the old lattice constant (production value
   `a=3.506186475312 A`);
4. do not treat low-temperature freezing as metastability proof;
5. do not use a single final snapshot as proof; keep trajectories;
6. do not use total energy as the only criterion;
7. do not replace dynamics with SOAP;
8. do not perform phase-transition classification for now;
9. the periodic bicrystal has two GBs; divide GB energy by two interfaces;
10. one successful trajectory is not an ensemble conclusion.

