# Li_Sigma5_310_metastable_2x2 (Bilingual / 双语版)

# 中文版

# Li_Sigma5_310_metastable_2x2（中文版）

基于已验证的旧版 Li Σ5 [001](310) 双晶模型，建立一个新的、隔离的 2×2
面内扩展实验工作流，用于检验高温重构晶界在 10 K / 0 bar 下能否稳定存在
超过 1 ps。

## 当前状态

| Stage | 内容 | 状态 |
|---|---|---|
| 0 | 2×2 模型生成与几何校验 | **PASS**，14400 atoms |
| 1 | 初始弛豫 | NOT RUN，脚本已备好 |
| 2 | 363 K NPT / 2 ps 高温重构 | NOT RUN，脚本已备好 |
| 3 | 363 K -> 10 K，10 K/ps 淬火 | NOT RUN，脚本已备好 |
| 4 | 10 K / 0 bar 稳定性保持 | NOT RUN，脚本已备好 |
| 5 | 结构与能量审计 | NOT RUN |
| 6 | 尺寸效应（3×3） | NOT RUN，待 2×2 成功后再做 |

本机当前没有可用的 LAMMPS（`lmp` 不存在，WSL 不可用），因此 Stage 1-4
的 MD 计算需要在一台装有 LAMMPS + OpenKIM 的机器上执行。所有输入脚本、
参数和校验文件已放在对应目录与 `lammps/` 中。

## 结论判定规则

只有同时满足以下条件才可判定为：

```text
SUCCESS - candidate metastable disordered GB
```

1. 宏观晶界仍为 Σ5 [001](310)；
2. 高温下晶界发生明显无序化/重构；
3. 晶粒内部保持 BCC；
4. 以 10 K/ps 淬火；
5. 10 K、0 bar 下保持；
6. 晶界无序结构持续超过 1 ps。

否则一律写为 `FAIL / INSUFFICIENT`，不得因为能量平坦、原子几乎不动或
SOAP 距离大而宣称成功。

完整协议见 `PROTOCOL.md`，实时执行记录见 `RUN_STATUS.md`。


# English Version

# Li_Sigma5_310_metastable_2x2 (English)

This is a new isolated production workflow built from the validated legacy
Li Sigma5 [001](310) bicrystal. It doubles the two in-plane GB dimensions
and tests whether a thermally reconstructed GB survives at 10 K / 0 bar for
more than 1 ps.

## Current status

| Stage | Content | Status |
|---|---|---|
| 0 | 2x2 model generation and geometry validation | **PASS**, 14400 atoms |
| 1 | Initial minimization | NOT RUN, scripts ready |
| 2 | 363 K NPT / 2 ps high-temperature reconstruction | NOT RUN, scripts ready |
| 3 | 363 K to 10 K at 10 K/ps quench | NOT RUN, scripts ready |
| 4 | 10 K / 0 bar stability hold | NOT RUN, scripts ready |
| 5 | Structure and energy audit | NOT RUN |
| 6 | Size effect (3x3) | NOT RUN, after 2x2 succeeds |

LAMMPS is not currently available on this machine (`lmp` is missing and WSL
is unavailable), so Stage 1-4 MD must be run on a machine with LAMMPS +
OpenKIM. All input scripts, parameters, and validation files are ready.

## Decision rule

The result may be written as:

```text
SUCCESS - candidate metastable disordered GB
```

only when all of the following hold:

1. the macroscopic GB remains Sigma5 [001](310);
2. the GB clearly disorders/reconstructs at high temperature;
3. the grain interiors remain BCC;
4. quenching is performed at 10 K/ps;
5. the system is held at 10 K and 0 bar;
6. the disordered GB structure survives for more than 1 ps.

Otherwise write `FAIL / INSUFFICIENT`. Do not claim success just because
energy is flat, atoms barely move, or SOAP distances are large.

See `PROTOCOL.md` for the full protocol and `RUN_STATUS.md` for the live
execution log.

