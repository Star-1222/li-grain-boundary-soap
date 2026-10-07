# MIGRATION_MAP — 旧文件到新目录的映射

旧权威项目：

```text
D:\Research\Lithium_GrainBoundary_Project\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production
```

| 旧位置 | 新位置 | 用途 |
|---|---|---|
| `docs/README.md` | `legacy/README_v2_final.md` | 旧版权威说明 |
| `src/generator/` | `generator/` 与 `legacy/generator_old/` | 生成器 v2.0.0 副本，源码未改 |
| `lammps/` | `legacy/lammps_v2/` | 旧版 LAMMPS 脚本参考 |
| `tests/` | `legacy/tests_v2/` | 旧版测试参考 |
| `docs/reports/` | `legacy/reports_v2/` | 旧版 FINAL_STATUS、ENVIRONMENT 等 |
| `requirements.txt` | `legacy/requirements_v2.txt` | Python 依赖参考 |

新实验目录根级 `README.md`、`PROTOCOL.md`、`PARAMETERS.json`、
`RUN_STATUS.md`、`PROVENANCE.md` 为本次新写文件；`lammps/` 下为本次
新协议脚本；`stage0_build_validate/` 为本次新生成的结构与校验结果。

规则：

1. 旧文件只读，不修改；
2. 旧输出不得冒充新实验结果；
3. 所有新 MD 输出必须写在本实验自己的 `stage*_*/` 目录。
