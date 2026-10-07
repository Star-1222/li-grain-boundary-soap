# 锂晶界与 SOAP 方法研究

BCC 锂 Σ5(310)[001] 对称倾转晶界的亚稳态构建，以及一套可迁移的 SOAP
局域结构分析方法工作流。

> 私密仓库，内容为未发表的研究工作，请勿外传。

---

## 项目组成

| 目录 | 内容 | 状态 |
|---|---|---|
| `Lithium_GrainBoundary_Project` | 锂晶界主项目：生成器、生产工作流、文献 | 基准模型已定稿（3600 原子，63 项测试通过） |
| `Li_Sigma5_310_2x2_阶段一验收_2026-09-05_new` | 2×2 面内扩大实验（14400 原子）阶段一验收 | 通过，可进入阶段二 |
| `Ni_SOAP_Method_Exploration_by_Stage` | Ni SOAP 方法探索，六阶段归档 | 流程可执行且可审计；参数与距离定义尚未冻结 |
| `_archive` | 历史版本与迁移记录 | 仅供溯源 |
| `tools` | 产物台账生成器、晶界能计算器、git 钩子安装器 | 可用 |

---

## 从这里看起

| 文档 | 内容 |
|---|---|
| [`STUDY_GUIDE.md`](STUDY_GUIDE.md) | ★ **从零理解本项目的学习路线**（7 章，含自检问题） |
| [`PROJECT_REVIEW.md`](PROJECT_REVIEW.md) | ★ **项目评估与改进路线**：缺口、优先级、行动清单 |
| [`ARTIFACT_LEDGER.md`](ARTIFACT_LEDGER.md) | 产物台账：每个文件的层次、用途、内容分析、校验和 |
| [`DATA_MANIFEST.md`](DATA_MANIFEST.md) | 未入库的大文件清单（约 466 MB），数据还原核对表 |
| [`tools/README.md`](tools/README.md) | 工具用法说明 |

---

## 核心结论（截至 2026-10-07）

**亚稳态晶界构建**

热重构后的无序晶界能量比理想 CSL 晶界低约 **12 mJ/m²**（153.2 vs 165.0
mJ/m²），说明理想 CSL 结构并非该体系在 0 K 下的基态。2×2 实验呈现同样
方向（187.0 → 160.8）。

> 该结果由现有日志复算得到，尚未用完全一致的弛豫协议复核，见
> `PROJECT_REVIEW.md` §1.1。

**阶段一验收**

原目标通过：363 K 退火 25 ps → 10 K/ps 淬火 → 10 K 保持 5 ps。
晶界 BCC 占比降至 38%，晶粒内部保持 98.7%，无序带宽 6.5 Å。

> 「0 bar」按标量时均压力口径成立（+1.7 bar），但三个应力分量为
> −187 / +386 / −196 bar，尚需 aniso NPT 敏感性复核。

**SOAP 方法线**

控制物理描述符后，相关系数由 0.83 降至 0.21；仅用无序度 D 的能量模型
R² 为 0.0393。另发现 D 与位移范数相关达 **0.99998** —— 若成立，则
SOAP 派生坐标相对简单几何量未提供额外信息，方法前提需优先核实。

---

## 环境

| 项目 | 值 |
|---|---|
| 计算平台 | WSL / Linux（Windows 上的 WSL2） |
| 分子动力学 | LAMMPS（2025-12-10 构建）+ OpenKIM，OpenMPI 12 ranks |
| 势函数 | `SNAP_ZuoChenLi_2019_Li__MO_732106099012_000` |
| 分析环境 | Windows，Python 3.12+（ase / numpy / scipy / matplotlib） |

---

## 版本控制范围

**已入库**：生成器与分析代码、LAMMPS 输入脚本、报告文档、运行日志、
参数记录、结构文件（`.data`/`.xyz`/`.cif`）、**LAMMPS 检查点（`.restart`）**、
历史归档。

**未入库**：文献 PDF、LAMMPS 轨迹（`.dump`）、SOAP 向量矩阵（`.npy`）、
部分大型数据表，合计约 466 MB。清单、用途与 SHA256 见
[`DATA_MANIFEST.md`](DATA_MANIFEST.md) —— 全新克隆时据此核对缺失的数据。

`.restart` 之所以入库：它保存速度场与热力学状态，而 `.data` 只写坐标，
丢失后必须从上一阶段重跑。

---

## 日常操作

```cmd
refresh_ledger.cmd          :: 刷新产物台账
git commit                  :: 提交前钩子会自动刷新台账
```

```bash
python tools/compute_gb_energy.py    # 从现有日志复算晶界过剩能
```

克隆到新机器后，需运行一次钩子安装脚本才会自动刷新台账：

```powershell
powershell -ExecutionPolicy Bypass -File tools\install_git_hook.ps1
```

---

## 行尾约定

本项目在 Windows 上编辑、在 WSL 下执行计算，行尾由 `.gitattributes` 固定：

- 文本文件与 `*.sh` 一律 LF —— 否则 WSL 下 `sh` 会报 `\r': command not found`
- `*.cmd` / `*.bat` 为 CRLF —— cmd.exe 解析 LF 行尾的批处理会出错
