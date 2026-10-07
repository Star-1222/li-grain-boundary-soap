# 研究产物台账工具

这套工具把 `D:\Research` 下的每一个产物与它的说明、层次、校验和、内容分析
绑定在一起，生成一份会自我更新的台账。

## 产出物

| 文件 | 内容 |
|---|---|
| `ARTIFACT_LEDGER.md` | 人读版台账：层次分布、自动识别的问题、内容分析、逐文件清单 |
| `ARTIFACT_LEDGER.csv` | 机器可读版：每行一个产物，含 SHA256，可用 Excel 筛选 |

文件名刻意使用 ASCII，避免中文路径在脚本与 git 钩子中出问题；内容仍是中文。

## 每个产物绑定了什么

1. **层次** —— 结论 / 续算 / 中间 / 追溯 / 代码 / 输入 / 文献 / 文档，并给出处置建议
2. **用途** —— 这个文件是什么、在项目里扮演什么角色
3. **内容分析** —— 直接从文件里读出来的事实，例如：
   - `.data` → 原子数、盒子尺寸、timestep
   - `.dump` → 帧数、timestep 范围、原子数
   - `.stats` → 压力分量的时均值与极差（会标注各向异性）
   - `.npy` / `.npz` → 数组形状、类型、是否全为有限值
   - `.csv` → 行列数、前若干列名
   - CNA 报告 → 各阶段晶界/晶内 BCC 比例与无序带宽
   - 能量模型表 → 各模型的 R² 均值
4. **SHA256** —— 完整性校验；文件一变就能看出来

另外台账会自动做**交叉核对**：把报告里声称的数字与原始数据实测值对照。
例如 `stage6_energy_models.csv` 的实测 R² 与最终报告声称的数值逐项一致。

## 三种更新方式

### 1. 手动（随时）

双击 `D:\Research\refresh_ledger.cmd`，或在命令行运行：

```cmd
D:\Research\refresh_ledger.cmd
```

### 2. 每次提交自动更新（推荐）

在仓库根目录执行一次：

```powershell
powershell -ExecutionPolicy Bypass -File tools\install_git_hook.ps1
```

之后每次 `git commit`，钩子会先刷新台账并把结果一并纳入本次提交。
若台账生成失败，钩子会放行（`exit 0`），不会阻塞你的提交。

### 3. 定时刷新

用任务计划程序定期调用 `refresh_ledger.cmd` 即可，例如每天一次。

## 换目录或用别的项目

通过环境变量指定扫描根：

```cmd
set RESEARCH_ROOT=E:\MyProject
D:\Research\refresh_ledger.cmd
```

`LEDGER_MIRROR` 为空字符串时不产生镜像副本。
如需登记新的项目目录，编辑 `tools\build_artifact_catalog.py` 里的 `PROJECTS` 列表。

## 注意

- 台账是**派生物**，不要手工编辑，改动会在下次刷新时丢失。
- 新增产物后无需登记，重新生成即可自动收录。
- `.cmd` 必须保持 CRLF、`.ps1` 必须保持纯 ASCII —— 这两点在本机都踩过坑。
