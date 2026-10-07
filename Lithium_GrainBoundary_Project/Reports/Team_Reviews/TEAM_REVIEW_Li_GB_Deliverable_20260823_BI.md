# Peer Delivery Review: Li_GB_Deliverable_20260823.zip (Bilingual / 双语版)

本文件为双语版：中文版在前，英文版在后。This file contains the Chinese version first, followed by the English version.

# 中文版

# 队友交付审查：Li_GB_Deliverable_20260823.zip

审查日期：2026-08-29

## 1. 审查范围

解压并检查了 `Li_GB_Deliverable_20260823.zip` 的全部 65 个文件，包括：

- 泛化双晶生成器 `src/generator/`
- LAMMPS 六阶段工作流 `lammps/`
- SOAP 分析包 `soap_analysis/`
- 测试 `tests/`
- 构建与验证报告 `build/Sigma5_310_2x/`
- 日志 `logs/`
- 环境脚本 `environment/`

压缩包内的 READM   E/CHANGES 中的用户指令只作为审查背景，不作为本轮执行指令。本轮只做交付物质量审查，没有修改队友的任何文件。

## 2. 总体结论

**暂不建议视为完成。**

生成器和 SOAP 分析包的结构、注释、测试意识都不错；生成器泛化也实测可用。但核心交付内容，即 LAMMPS 高温淬火稳定性验证，**没有实际跑通**，且交付包存在若干一致性与科学方法问题。

## 3. 已验证通过的部分

### 3.1 泛化双晶生成器

使用交付包自带的生成器独立构建了 5 个平面：

| 平面 | Sigma | theta (deg) | 原子数 | 几何校验 |
|---|---|---|---|---|
| 310 | 5 | 36.87 | 14400 | PASS |
| 510 | 13 | 22.62 | 14976 | PASS |
| 210 | 5 | 53.13 | 14400 | PASS |
| 320 | 13 | 67.38 | 14560 | PASS |
| 410 | 17 | 28.07 | 14688 | PASS |

对交付的 `Sigma5_310_K12_J6_L10.data` 重新执行 `generator.verify`：

- 原子数 14400
- 晶粒 1/2 各 7200
- 取向差角 36.86989765 deg
- 最小距离 2.219919 Angstrom
- PASS

### 3.2 SOAP 环境与算法自检

`environment/verify_soap_env.py` 在当前可用依赖下 7/7 通过：

- Python 3.13
- numpy / scipy
- SOAP 后端：dscribe + numpy
- scipy 球谐函数参数顺序
- 23 度绕 [1,2,3] 轴旋转不变性
- dscribe 全矩阵计算

SOAP numpy 后端可以独立读取生成的 Sigma5 data 文件并输出 disorder 剖面。

## 4. 主要问题

### 4.1 P0：LAMMPS 热处理工作流没有实际跑完

证据：

- `logs/run_all_console.log` 显示：
  - `[skip] stage0_check`
  - `[run ] stage1_minimize ...`
  - `job aborted`
- `logs/stage1_minimize.lammps` 只到 chunk 9，fmax 约 0.00115，未达到 1e-4；
  没有 `STAGE1_DONE`，没有 `STAGE1_MINIMIZE`。
- 没有 `stage2_anneal.lammps`、`stage3_quench.lammps`、`stage4_hold.lammps` 日志。
- 没有 `stage4_held.data`、`stage4_held.restart`、`STAGE4_HOLD_STATS`。
- README 已明确说明“新热处理协议尚未用 LAMMPS 实际运行验证”。

结论：核心科学验证，即“363 K 高温保持 -> 10 K/ps 降温 -> 10 K 零压保持”的稳定性判据，尚未交付。

### 4.2 P1：`run_all.sh` 头部与实际情况矛盾

`lammps/run_all.sh` 开头写：

```
Verified Stage 0 -> 5 workflow with the high-T-quench thermal protocol.
```

但实际日志显示 Stage 1 即崩溃，README 也说新协议未运行。建议删除或改为“待 HPC 验证”。

### 4.3 P1：交付包缺少 `modified/` 镜像目录

README 和 CHANGES 均说明应包含 `modified/` 目录，其中应保存全部修改文件的同名镜像，并注明共 20 个文件。实际 zip 中没有该目录。

### 4.4 P1：测试套件不是自足的，无法从 zip 独立复现 25/25

以下测试依赖包外旧生产数据，文件缺失时会 skip：

- `test_sigma5_positions_match_archived_data`
- `test_gb_structure_smoke`

当前 zip 中没有：

- `01_GB_Construction/Sigma5_310_Metastable_GB/v2_production/...`
- `stage4_final.data` 或 `stage5_final.data`

因此仅凭交付 zip，无法证明“25/25 全部通过”。

### 4.5 P1：SOAP 周期边界导致 GB 区域重复计数

对重新生成的 Sigma5 运行 SOAP 分析，`gb_width_regions` 返回 3 个区域：

```
[0.5, 1.5]  w=2.0A
[65.5, 67.5] w=3.0A
[132.5, 132.5] w=1.0A
```

原因是同一个 x=0 晶界在周期边界另一侧被重复识别。`gb_width_regions` 只在首尾两个 bin 都高于阈值时才合并，当最后一个 bin 为空时不会合并。

建议修复周期边界处理，并补充一个“新生成 Sigma5 恰好检测到 2 个周期 GB 区域”的测试。

### 4.6 P1：SOAP 体相参考与批量流程不一致

代码默认的松弛体相参考晶格常数为：

```
DEFAULT_A = 3.506186475312
```

但 `workflow/run_gb_series.sh` 会把生成器 `params.json` 传给 SOAP：

```
--params "$builddir/params.json"
```

其中 `a_A = 3.51` 会覆盖默认参考。也就是说，批量 SOAP 分析实际使用未松弛的 3.51 Angstrom 做体相参考，而不是文档记录的松弛体相参考。

建议：明确指定使用松弛后的体相晶格常数，或先运行 `in.bulk_energy` 再传入实测值。

### 4.7 P2：`resume.sh stage1` 不是真正的断点续跑

`resume.sh` 的 stage1 分支读取 `stage1_checkpoint.data` 后仍从 `chunk=1` 重新运行 `in.1_minimize`，会重复执行之前已经做过的极小化循环。

建议改为真正的 restart 或至少从记录的 chunk 编号继续。

### 4.8 P2：依赖版本未锁定

`environment/requirements_soap.txt` 使用 `>=` 版本范围：

```
numpy>=1.24
scipy>=1.10
pytest>=7.0
ase>=3.22
dscribe>=2.0
```

README 声称“含已验证版本号”，但文件并未 pin 具体版本，不利于复现。

### 4.9 P2：交付包中缺少 Stage 1 完成产物

当前 `data/` 下只有：

- `Sigma5_310.data`
- `stage1_checkpoint.data`
- `stage1_checkpoint.restart`

缺少：

- `stage1_minimized.data`
- `stage1_minimized.restart`
- `stage1_minimized.peratom`

Stage 2 脚本读取 `stage1_minimized.restart`，因此当前交付包不能从中断点直接继续 Stage 2。

## 5. 建议的补齐清单

建议队友补齐以下内容后再重新交付：

1. 完整跑通 Stage 0-4，提交成功日志：
   - `stage0_check.lammps`
   - `stage1_minimize.lammps`
   - `stage2_anneal.lammps`
   - `stage3_quench.lammps`
   - `stage4_hold.lammps`
2. 提交 Stage 4 终态产物：
   - `stage4_held.data`
   - `stage4_held.restart`
   - `STAGE4_HOLD_STATS`
3. 提交 Stage 1 完成产物：
   - `stage1_minimized.data`
   - `stage1_minimized.restart`
   - `stage1_minimized.peratom`
4. 补齐 `modified/` 镜像目录。
5. 修复 `gb_width_regions` 周期边界重复计数。
6. 统一 SOAP 体相参考晶格常数。
7. 提供可独立运行测试所需的结构数据，或把依赖的归档文件一并打包。
8. 修正 `run_all.sh` 的 “Verified” 表述。
9. 将依赖版本固定。
10. 移除交付包根目录的临时调试脚本 `_check_dll.py`、`_check_kim.py`（或移入 `tools/` 并说明用途）。

## 6. 最终意见

**代码质量和工程习惯：中上。**

**交付完成度：未完成。**

生成器扩展和 SOAP 分析包可以作为下一步工作的基础，但 LAMMPS 高温淬火稳定性验证是本轮交付的核心目标，目前没有成功运行的证据。建议先补齐完整运行结果和上述修复，再进入下一阶段。


# English Version

# Peer Delivery Review: Li_GB_Deliverable_20260823.zip

Review date: 2026-08-29

## 1. Review scope

All 65 files in `Li_GB_Deliverable_20260823.zip` were extracted and checked,
including:

- generalized bicrystal generator `src/generator/`
- six-stage LAMMPS workflow `lammps/`
- SOAP analysis package `soap_analysis/`
- tests `tests/`
- build and verification reports `build/Sigma5_310_2x/`
- logs `logs/`
- environment scripts `environment/`

Instructions in the README/CHANGES inside the package were treated only as
review background, not as execution instructions for this turn. This review
only assessed deliverable quality and did not modify any of the teammate's
files.

## 2. Overall verdict

**Not yet recommended as complete.**

The generator and SOAP analysis package are well structured, commented, and
test-conscious; generator generalization also works in practice. However, the
core deliverable, the LAMMPS high-temperature quench stability validation,
was **not actually completed**, and the package has several consistency and
methodological problems.

## 3. Parts verified as passing

### 3.1 Generalized bicrystal generator

Five planes were independently constructed with the packaged generator:

| Plane | Sigma | theta (deg) | Atoms | Geometry check |
|---|---|---|---|---|
| 310 | 5 | 36.87 | 14400 | PASS |
| 510 | 13 | 22.62 | 14976 | PASS |
| 210 | 5 | 53.13 | 14400 | PASS |
| 320 | 13 | 67.38 | 14560 | PASS |
| 410 | 17 | 28.07 | 14688 | PASS |

`generator.verify` was rerun on the delivered `Sigma5_310_K12_J6_L10.data`:

- atoms: 14400
- grains 1/2: 7200 each
- misorientation angle: 36.86989765 deg
- minimum distance: 2.219919 Angstrom
- PASS

### 3.2 SOAP environment and algorithm self-check

`environment/verify_soap_env.py` passed 7/7 with the currently available
dependencies:

- Python 3.13
- numpy / scipy
- SOAP backend: dscribe + numpy
- scipy spherical-harmonic argument order
- 23-degree rotation invariance about the [1,2,3] axis
- full dscribe matrix calculation

The SOAP numpy backend can independently read the generated Sigma5 data file
and produce a disorder profile.

## 4. Main findings

### 4.1 P0: The LAMMPS thermal workflow was not actually completed

Evidence:

- `logs/run_all_console.log` shows:
  - `[skip] stage0_check`
  - `[run ] stage1_minimize ...`
  - `job aborted`
- `logs/stage1_minimize.lammps` reaches only chunk 9 with fmax around
  0.00115, not 1e-4; no `STAGE1_DONE` and no `STAGE1_MINIMIZE`.
- There are no `stage2_anneal.lammps`, `stage3_quench.lammps`, or
  `stage4_hold.lammps` logs.
- There are no `stage4_held.data`, `stage4_held.restart`, or
  `STAGE4_HOLD_STATS`.
- The README explicitly states that the new thermal protocol has not yet been
  validated with an actual LAMMPS run.

Conclusion: the core scientific validation, namely the stability criterion
for "363 K high-temperature hold -> 10 K/ps cooling -> 10 K zero-pressure
hold", has not been delivered.

### 4.2 P1: `run_all.sh` header contradicts the actual state

`lammps/run_all.sh` starts with:

```
Verified Stage 0 -> 5 workflow with the high-T-quench thermal protocol.
```

But the actual logs show a crash at Stage 1, and the README says the new
protocol has not run. The statement should be deleted or changed to "awaiting
HPC validation".

### 4.3 P1: The package is missing the `modified/` mirror directory

Both README and CHANGES state that a `modified/` directory should contain
mirror copies of all modified files (20 files total). The actual zip does not
contain this directory.

### 4.4 P1: The test suite is not self-contained and cannot reproduce 25/25
from the zip alone

The following tests depend on old production data outside the package and are
skipped when files are missing:

- `test_sigma5_positions_match_archived_data`
- `test_gb_structure_smoke`

The zip does not contain:

- `01_GB_Construction/Sigma5_310_Metastable_GB/v2_production/...`
- `stage4_final.data` or `stage5_final.data`

Therefore, "25/25 passed" cannot be proven from the deliverable zip alone.

### 4.5 P1: SOAP periodic boundary causes duplicate counting of the GB region

Running SOAP analysis on a regenerated Sigma5 structure,
`gb_width_regions` returns 3 regions:

```
[0.5, 1.5]  w=2.0A
[65.5, 67.5] w=3.0A
[132.5, 132.5] w=1.0A
```

The same x=0 GB is identified twice across the periodic boundary.
`gb_width_regions` only merges when both the first and last bins exceed the
threshold, so it does not merge when the last bin is empty.

Recommendation: fix periodic-boundary handling and add a test that a newly
generated Sigma5 structure detects exactly two periodic GB regions.

### 4.6 P1: SOAP bulk reference is inconsistent with the batch workflow

The default relaxed bulk reference lattice constant in the code is:

```
DEFAULT_A = 3.506186475312
```

But `workflow/run_gb_series.sh` passes the generator `params.json` to SOAP:

```
--params "$builddir/params.json"
```

where `a_A = 3.51` overrides the default reference. The batch SOAP analysis
therefore uses the unrelaxed 3.51 Angstrom bulk reference instead of the
documented relaxed reference.

Recommendation: explicitly use the relaxed bulk lattice constant, or run
`in.bulk_energy` first and pass the measured value.

### 4.7 P2: `resume.sh stage1` is not a true checkpoint resume

The stage1 branch of `resume.sh` reads `stage1_checkpoint.data` but then
reruns `in.1_minimize` from `chunk=1`, repeating minimization loops that have
already been completed.

Recommendation: implement a real restart, or at least continue from the
recorded chunk number.

### 4.8 P2: Dependency versions are not pinned

`environment/requirements_soap.txt` uses `>=` ranges:

```
numpy>=1.24
scipy>=1.10
pytest>=7.0
ase>=3.22
dscribe>=2.0
```

The README claims "verified version numbers", but the file does not pin
specific versions, which hurts reproducibility.

### 4.9 P2: Stage 1 completion artifacts are missing from the package

The current `data/` contains only:

- `Sigma5_310.data`
- `stage1_checkpoint.data`
- `stage1_checkpoint.restart`

Missing:

- `stage1_minimized.data`
- `stage1_minimized.restart`
- `stage1_minimized.peratom`

Since the Stage 2 script reads `stage1_minimized.restart`, the delivered
package cannot continue directly from the checkpoint into Stage 2.

## 5. Suggested completion checklist

Ask the teammate to resubmit after completing:

1. Run Stages 0-4 to completion and submit successful logs:
   - `stage0_check.lammps`
   - `stage1_minimize.lammps`
   - `stage2_anneal.lammps`
   - `stage3_quench.lammps`
   - `stage4_hold.lammps`
2. Submit Stage 4 final artifacts:
   - `stage4_held.data`
   - `stage4_held.restart`
   - `STAGE4_HOLD_STATS`
3. Submit Stage 1 completion artifacts:
   - `stage1_minimized.data`
   - `stage1_minimized.restart`
   - `stage1_minimized.peratom`
4. Add the `modified/` mirror directory.
5. Fix duplicate GB counting in `gb_width_regions` across periodic
   boundaries.
6. Unify the SOAP bulk-reference lattice constant.
7. Provide the structure data required for standalone tests, or package the
   archived dependencies.
8. Correct the "Verified" wording in `run_all.sh`.
9. Pin dependency versions.
10. Remove temporary debug scripts `_check_dll.py` and `_check_kim.py` from
    the package root, or move them into `tools/` with an explanation.

## 6. Final opinion

**Code quality and engineering habits: upper-middle.**

**Deliverable completeness: incomplete.**

The generator extension and SOAP analysis package can serve as the basis for
the next step, but the LAMMPS high-temperature quench stability validation is
the core goal of this delivery and there is no evidence of a successful run.
Complete runs and the fixes above should be supplied before the next stage.

