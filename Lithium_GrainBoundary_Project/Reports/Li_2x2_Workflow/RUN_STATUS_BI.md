# RUN_STATUS (Bilingual / 双语版)

# 中文版

# RUN_STATUS（中文版）

更新于 2026-08-30。本文件是实验的实时执行记录。

| Stage | 状态 | 关键结果 / 备注 |
|---|---|---|
| 0 | **PASS** | 14400 atoms；box 66.525211 × 66.525211 × 70.123730 Å；两条等价 GB（x=0, 33.2626 Å）；bulk BCC 10080/10080；min dist 2.217507 Å；Sigma5/misorientation 校验通过 |
| 1 | NOT RUN | `lammps/in.1_minimize` 已备好；本机无 LAMMPS |
| 2 | NOT RUN | `lammps/in.2_highT_npt` 已备好；363 K NPT / 2 ps |
| 3 | NOT RUN | `lammps/in.3_quench` 已备好；363→10 K @ 10 K/ps，35.3 ps |
| 4 | NOT RUN | `lammps/in.4_10K_hold` 已备好；10 K / 0 bar，5 ps |
| 5 | NOT RUN | 待 Stage 4 通过 |
| 6 | NOT RUN | 待 2×2 成功 |

## 在有 LAMMPS + OpenKIM 的机器上运行 MD

Linux/WSL：

```bash
cd Li_Sigma5_310_metastable_2x2
OUTDIR="$PWD" \
LOGDIR="$PWD/logs" \
MPI_NP=12 \
OMP_NUM_THREADS=1 \
bash lammps/run_all.sh
```

Windows 本地 LAMMPS：

```powershell
lmp -in lammps/in.0_check -var outdir stage0_build_validate -var kim_model SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
lmp -in lammps/in.1_minimize -var outdir stage1_relax -var input_data stage0_build_validate/Sigma5_310_K6_J6_L20.data -var kim_model SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
lmp -in lammps/in.2_highT_npt -var outdir stage2_highT_NPT -var seed 12345 -var kim_model SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
lmp -in lammps/in.3_quench -var outdir stage3_quench -var kim_model SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
lmp -in lammps/in.4_10K_hold -var outdir stage4_10K_hold -var kim_model SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
```

注意 `in.2_highT_npt` 读取 `stage1_relax/stage1_minimized.restart`，
`in.3_quench` 读取 `stage2_highT_NPT/stage2_highT.restart`，
`in.4_10K_hold` 读取 `stage3_quench/stage3_quenched.restart`，因此必须
按顺序执行。恢复中断时可使用 `lammps/resume.sh stageX`。

## Stage 2 判定清单（淬火前必须人工检查）

- [ ] t=0, 0.5, 1, 1.5, 2 ps 帧是否已查看；
- [ ] 晶界是否明显无序/重构（Case B）；
- [ ] 晶粒内部是否仍为 BCC；
- [ ] 若 Case A → 不得进入 Stage 3；若 Case C → FAIL。

## 最终结论模板

```text
SUCCESS — candidate metastable disordered GB
```

或

```text
FAIL / INSUFFICIENT — do not overinterpret
```


# English Version

# RUN_STATUS (English)

Updated 2026-08-30. This file is the live execution log of the experiment.

| Stage | Status | Key results / notes |
|---|---|---|
| 0 | **PASS** | 14400 atoms; box 66.525211 x 66.525211 x 70.123730 A; two equivalent GBs (x=0, 33.2626 A); bulk BCC 10080/10080; min dist 2.217507 A; Sigma5/misorientation checks passed |
| 1 | NOT RUN | `lammps/in.1_minimize` ready; no LAMMPS on this machine |
| 2 | NOT RUN | `lammps/in.2_highT_npt` ready; 363 K NPT / 2 ps |
| 3 | NOT RUN | `lammps/in.3_quench` ready; 363 to 10 K @ 10 K/ps, 35.3 ps |
| 4 | NOT RUN | `lammps/in.4_10K_hold` ready; 10 K / 0 bar, 5 ps |
| 5 | NOT RUN | after Stage 4 passes |
| 6 | NOT RUN | after 2x2 succeeds |

## Running MD on a machine with LAMMPS + OpenKIM

Linux/WSL:

```bash
cd Li_Sigma5_310_metastable_2x2
OUTDIR="$PWD" \
LOGDIR="$PWD/logs" \
MPI_NP=12 \
OMP_NUM_THREADS=1 \
bash lammps/run_all.sh
```

Windows local LAMMPS:

```powershell
lmp -in lammps/in.0_check -var outdir stage0_build_validate -var kim_model SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
lmp -in lammps/in.1_minimize -var outdir stage1_relax -var input_data stage0_build_validate/Sigma5_310_K6_J6_L20.data -var kim_model SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
lmp -in lammps/in.2_highT_npt -var outdir stage2_highT_NPT -var seed 12345 -var kim_model SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
lmp -in lammps/in.3_quench -var outdir stage3_quench -var kim_model SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
lmp -in lammps/in.4_10K_hold -var outdir stage4_10K_hold -var kim_model SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
```

Note: `in.2_highT_npt` reads `stage1_relax/stage1_minimized.restart`,
`in.3_quench` reads `stage2_highT_NPT/stage2_highT.restart`, and
`in.4_10K_hold` reads `stage3_quench/stage3_quenched.restart`, so they must
run in order. Use `lammps/resume.sh stageX` for interrupted runs.

## Stage 2 manual checklist before quenching

- [ ] have t=0, 0.5, 1, 1.5, 2 ps frames been inspected?
- [ ] is the GB clearly disordered/reconstructed (Case B)?
- [ ] do grain interiors remain BCC?
- [ ] if Case A, do not enter Stage 3; if Case C, FAIL.

## Final conclusion template

```text
SUCCESS - candidate metastable disordered GB
```

or

```text
FAIL / INSUFFICIENT - do not overinterpret
```

