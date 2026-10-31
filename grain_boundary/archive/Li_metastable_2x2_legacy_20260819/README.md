# Li_Sigma5_310_metastable_2x2

基于已验证的旧版 Li Σ5 [001](310) 双晶模型，建立一个新的、隔离的 2×2 面内扩展实验工作流，用于检验高温重构晶界在 10 K / 0 bar 下能否稳定存在超过 1 ps。

This is a new isolated production workflow built from the validated legacy
Li Σ5 [001](310) bicrystal. It doubles the two in-plane GB dimensions and
tests whether a thermally reconstructed GB survives at 10 K / 0 bar for
more than 1 ps.

## 当前状态 / Current status

| Stage | 内容 Content | 状态 Status |
|---|---|---|
| 0 | 2×2 模型生成与几何校验 | **PASS**，14400 atoms |
| 1 | 初始弛豫 Minimization | NOT RUN，脚本已备好 |
| 2 | 363 K NPT / 2 ps 高温重构 | NOT RUN，脚本已备好 |
| 3 | 363 K → 10 K，10 K/ps 淬火 | NOT RUN，脚本已备好 |
| 4 | 10 K / 0 bar 稳定性保持 | NOT RUN，脚本已备好 |
| 5 | 结构与能量审计 | NOT RUN |
| 6 | 尺寸效应（3×3） | NOT RUN，待 2×2 成功后再做 |

本机当前没有可用的 LAMMPS（`lmp` 不存在，WSL 不可用），因此 Stage 1–4
的 MD 计算需要在一台装有 LAMMPS + OpenKIM 的机器上执行。所有输入脚本、
参数和校验文件已经放在对应 `stage*_*/` 与 `lammps/` 目录中。

## 目录结构 / Layout

```text
Li_Sigma5_310_metastable_2x2/
├── README.md
├── PROMPT_FOR_CODEX.md          ← 原始任务书
├── PROTOCOL.md                  ← 完整 Stage 0–6 协议
├── PARAMETERS.json
├── RUN_STATUS.md
├── PROVENANCE.md
├── legacy/                      ← 旧版权威项目只读参考
│   ├── README_v2_final.md
│   ├── MIGRATION_MAP.md
│   ├── README_LEGACY_CORE.md
│   ├── generator_old/
│   ├── lammps_v2/
│   ├── tests_v2/
│   └── reports_v2/
├── generator/                   ← v2.0.0 生成器副本（未修改源码）
├── lammps/                      ← 新协议 LAMMPS 输入
├── scripts/
├── stage0_build_validate/
├── stage1_relax/
├── stage2_highT_NPT/
├── stage3_quench/
├── stage4_10K_hold/
├── stage5_analysis/
└── stage6_size_effect/
```

## 结论判定规则 / Decision rule

只有同时满足以下条件才可判定为：

```text
SUCCESS — candidate metastable disordered GB
```

1. 宏观晶界仍为 Σ5 [001](310)；
2. 高温下晶界发生明显无序化/重构；
3. 晶粒内部保持 BCC；
4. 以 10 K/ps 淬火；
5. 10 K、0 bar 下保持；
6. 晶界无序结构持续超过 1 ps。

否则一律写为 `FAIL / INSUFFICIENT`，不得因为能量平坦、原子几乎不动或
SOAP 距离大而宣称成功。

See `PROTOCOL.md` for the full protocol and `RUN_STATUS.md` for the live
execution log.
