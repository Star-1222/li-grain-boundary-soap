# -*- coding: utf-8 -*-
"""Build an artifact catalog for the repository.

Walks the research tree, classifies every file by role and layer, attaches a
plain-language purpose, and emits both a Markdown ledger and a CSV manifest
with SHA256 checksums.
"""

import csv
import hashlib
import os
import sys
from collections import defaultdict
from pathlib import Path

def _default_root() -> Path:
    """Repository root, derived from this file's location.

    The script normally lives at <repo>/maintenance/tools/build_artifact_catalog.py, so the
    repository root is three levels up. Override with RESEARCH_ROOT when the
    tree lives somewhere else.
    """
    env = os.environ.get("RESEARCH_ROOT")
    if env:
        return Path(env)
    here = Path(__file__).resolve().parent
    return here.parents[1]


ROOT = _default_root()

# Optional mirror copy of the generated ledger; empty disables it.
_mirror_env = os.environ.get("LEDGER_MIRROR")
OUT_DIR = Path(_mirror_env) if _mirror_env else None

PROJECTS = [
    ("grain_boundary", "主线一：晶界构造"),
    ("soap", "主线二：SOAP 方法"),
    ("maintenance", "仓库维护：文档、工具与记录"),
    ("", "仓库首页与其他文件"),
]


def iter_project_files(proj):
    """Yield files belonging to one project group.

    The empty key includes all files outside the explicitly named projects.
    New experiment and documentation directories must not silently disappear.
    artifact_ledger.py additionally filters this inventory through git ls-files.
    """
    if proj == "":
        named = {name for name, _ in PROJECTS if name}
        excluded = {".git", "__pycache__", ".venv", "venv", ".pytest_cache", "runs"}
        for p in sorted(ROOT.iterdir()):
            if p.is_file():
                yield p
            elif p.is_dir() and p.name not in named | excluded:
                for child in sorted(p.rglob("*")):
                    if child.is_file() and not excluded.intersection(child.relative_to(ROOT).parts):
                        yield child
        return
    base = ROOT / proj
    if not base.is_dir():
        return
    for p in sorted(base.rglob("*")):
        if p.is_file():
            yield p

LAYERS = {
    "结论": "结论性 — 必须保留，丢了要重跑",
    "续算": "续算检查点 — 删了要从上游重跑才能恢复",
    "中间": "中间/可再生 — 有上游输入和脚本就能重算，代价是算力",
    "追溯": "可追溯性 — 体积小但关键，记录“结果怎么来的”",
    "代码": "代码",
    "输入": "输入/配置",
    "文献": "文献",
    "文档": "报告与文档",
}

# --- purpose rules -----------------------------------------------------------
# Exact relative-path suffix matches take priority, then basename, then rules.

PATH_PURPOSE = {
    r"grain_boundary\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\initial\Sigma5_310.data":
        "初始双晶模型的 LAMMPS 格式（3600 原子）。整个项目的结构基准，后续每一级计算都从这里出发。",
    r"grain_boundary\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\initial\Sigma5_310.xyz":
        "同一模型的扩展 XYZ 格式，比 .data 更通用，便于其他软件读取。",
    r"grain_boundary\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\initial\Sigma5_310_grain_id.xyz":
        "每原子标注属于哪个晶粒。用来验证两条晶界在晶体学上等价——这是“对称倾转晶界”成立的前提。",
    r"grain_boundary\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\initial\Sigma5_310.cif":
        "晶体学标准格式。给外部软件或合作者读，通用性最好。",
    r"grain_boundary\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\initial\Sigma5_310.grain_map":
        "晶粒归属图，用于快速核对两条晶界的位置和对称性。",
    r"grain_boundary\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\initial\params.json":
        "构建参数：超胞尺寸 K/J/L、晶格常数、势函数名。复现这个模型所需的最小信息集。",
    r"grain_boundary\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\initial\removal_log.csv":
        "重合位点删除记录：每个晶界删了哪些原子、删了多少。保证删除对称、可审计，是可追溯性的核心。",
    r"grain_boundary\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\validation\rebuilt_stage2.data":
        "从 dump 轨迹反推重建的结构。用来证明 dump 没有丢信息——是“模拟结果可复现”的直接证据。",
    r"grain_boundary\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\docs\reports\SHA256SUMS.txt":
        "本包内文件的校验清单。日后任何人可据此确认文件未被改动。",
    r"grain_boundary\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\docs\reports\ENVIRONMENT.txt":
        "运行环境快照（LAMMPS 版本、MPI、势函数等），复现的必要条件。",
    r"grain_boundary\experiments\phase1_2026-09\Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\5_保持10K零压\p1_held.data":
        "★ 验收的核心产物：10 K 保持结束时的结构。此时晶界无序、晶粒内部仍为 BCC。",
    r"grain_boundary\experiments\phase1_2026-09\Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\5_保持10K零压\p1_held.restart":
        "★ 阶段二的起点。验收报告明确指定用这个检查点继续后续计算。",
    r"grain_boundary\experiments\phase1_2026-09\Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\5_保持10K零压\p1_hold_pressure.stats":
        "★ 保持段的压力时间序列。“10 K、0 bar 下维持”这个结论的直接证据来源。",
    r"grain_boundary\experiments\phase1_2026-09\Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\6_结构分析\cna_phase1_report.txt":
        "★ CNA 结构分析结论：各阶段 BCC 分数（晶界约 38%、晶粒内部 98.7%）。验收判据“晶界真实无序、内部保持 BCC”靠它成立。",
    r"grain_boundary\experiments\phase1_2026-09\阶段一验收汇报_2026-09-04.md":
        "★ 阶段一验收结论：原目标通过，可进入阶段二；并记录了 0 bar 口径的保留意见。",
    r"soap\Stage02_SOAP表示与审计\03_输出\stage2_75_soap_vectors.npy":
        "★ SOAP 主矩阵（4659 原子 × 1014 维，36 MB）。后续全部分析都建立在这个矩阵上；重算耗时较长。",
    r"soap\Stage02_SOAP表示与审计\03_输出\stage2_75_pairwise_sample.csv":
        "原子对距离样本表（48.8 MB）。支撑“原始距离里 42–46% 来自范数差”这一结论，属中间表。",
    r"soap\Stage02_SOAP表示与审计\03_输出\stage2_75_distance_decomposition.png":
        "★ 距离分解图：把 SOAP 距离拆成范数贡献与形状贡献。这是“原始距离含义不纯”的关键证据图。",
    r"soap\Stage03_LAE阈值与LAE_LER构建\03_输出\stage3a_epsilon_sweep.csv":
        "epsilon 阈值扫描结果。证明了阈值轻微变动就会改变分类数——离散 LAE 不稳。",
    r"soap\Stage03_LAE阈值与LAE_LER构建\03_输出\stage3b_partition_ari.csv":
        "★ 不同随机顺序下分区的一致性（ARI）。这是判定贪心 LAE “有顺序依赖、不可作为物理分类”的量化依据。",
    r"soap\Stage04_无序度与能量关联\03_输出\stage4a_atomwise_disorder.csv":
        "★ 逐原子无序度 D。本项目提出的连续结构坐标，是后续所有关联分析的自变量。",
    r"soap\Stage04_无序度与能量关联\03_输出\stage4b_atomwise_energy_disorder.csv":
        "逐原子的无序度与能量配对表，用于考察 D 与能量的关系。",
    r"soap\Stage05_物理空间与方向\03_输出\stage5a_atomwise_physical_descriptors.csv":
        "★ 常规物理描述符（配位数、键应变、体积、CSP）。最终证明“解释力主要来自这些量而非 SOAP 方向”的就是它。",
    r"soap\Stage05_物理空间与方向\03_输出\stage5b_matched_direction_energy_pairs.csv":
        "★ 配对重采样中间表（142.9 MB，占本项目近一半体积）。最终结论只用到由它算出的相关系数，本体属可再生的中间产物。",
    r"soap\Stage06_最终评估\03_输出\stage6_correlations.csv":
        "★ 全项目最终相关系数：原始 0.8317，控制物理描述符后降到 0.2058。这是方法线最重要的产出。",
    r"soap\Stage06_最终评估\03_输出\stage6_energy_models.csv":
        "★ 能量模型 R² 对比：仅 D 为 0.0393，加入物理描述符后 0.9869。量化了 SOAP 坐标的增量信息有限。",
    r"soap\00_汇总报告\Ni_SOAP_最终报告.md":
        "★ 方法探索总报告：项目定位、六个 Stage 的结论、十条下一版建议、以及“不能确定”清单。",

    # --- 锂晶界主项目：热历史各级产物 ---
    r"grain_boundary\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\bulk\bulk_Li_minimized.data":
        "★ 体相参考的零压弛豫结果。给出该势函数下的平衡晶格常数 3.506186 Å —— 全部双晶都用它构建，避免把体相应变强加给晶界周期。",
    r"grain_boundary\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\minimized\stage1_minimized.data":
        "Stage 1 极小化后的结构：初始双晶在 0 K 下消除残余力的结果。",
    r"grain_boundary\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\annealed\stage2_annealed.data":
        "★ Stage 2 退火末态：363 K、0 bar 下 25 ps 后的结构。晶界在此阶段发生热重构。",
    r"grain_boundary\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\annealed\stage2_anneal.dump":
        "退火过程轨迹。用来观察晶界重构是怎样逐步发生的，而不只是看首末两帧。",
    r"grain_boundary\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\quenched\stage3_quenched.data":
        "Stage 3 淬火末态：以 10 K/ps 降到 10 K 后的结构。无序被快速“冻”住。",
    r"grain_boundary\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\quenched\stage3_quench.dump":
        "淬火过程轨迹。用于核对降温速率与系统响应是否正常。",
    r"grain_boundary\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\final\stage4_final.data":
        "★ 全项目主产物：最终弛豫结构。单条退火-淬火-极小化轨迹收敛到的一个局域极小候选。",
    r"grain_boundary\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\final\stage4_final.restart":
        "最终结构的续算检查点，含速度与热力学状态。",

    # --- 阶段一验收包：逐步产物 ---
    r"grain_boundary\experiments\phase1_2026-09\0_初始XYZ与检测模块\Sigma5_310_K6_J6_L20.xyz":
        "14400 原子 2×2 初始模型（K=6, J=6, L=20）。本验收包所有计算的起点。",
    r"grain_boundary\experiments\phase1_2026-09\0_初始XYZ与检测模块\Sigma5_310_K6_J6_L20_grain_id.xyz":
        "初始模型带晶粒标注的版本，用于确认两条晶界等价、晶粒内部无缺陷。",
    r"grain_boundary\experiments\phase1_2026-09\0_初始XYZ与检测模块\README.md":
        "初始结构与检测模块的说明文档。",
    r"grain_boundary\experiments\phase1_2026-09\Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\1_零应力弛豫0K\p1_min0.data":
        "验收步骤 A：0 K 全维度零应力弛豫的末态。判据是三个方向压力约 0、最大力 8.5e-6 eV/Å。",
    r"grain_boundary\experiments\phase1_2026-09\Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\2_退火363K零压\p1_annealed.data":
        "★ 验收步骤 B：363 K 等压退火 25 ps 后的末态。晶界重构就是在这 25 ps 里发生的。",
    r"grain_boundary\experiments\phase1_2026-09\Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\2_退火363K零压\p1_anneal.dump":
        "退火 25 ps 的过程轨迹。报告里“末 16 ps 时均 T=362.8 K”就是从这里统计出来的。",
    r"grain_boundary\experiments\phase1_2026-09\Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\3_末帧0K检查\p1_casecheck_0K.data":
        "★ 验收步骤 C：取退火末帧做 0 K 检查。这是判定“Case B”（晶界真实无序、晶粒内部仍为 BCC）的关键结构。",
    r"grain_boundary\experiments\phase1_2026-09\Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\4_淬火10K\p1_quenched.data":
        "验收步骤 D：以 10 K/ps 从 363 K 淬火到 10 K 后的末态，历时 35.3 ps。",
    r"grain_boundary\experiments\phase1_2026-09\Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\4_淬火10K\p1_quench.dump":
        "淬火过程轨迹（35.3 ps）。用于核实降温速率是否严格为 10 K/ps。",
    r"grain_boundary\experiments\phase1_2026-09\Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\5_保持10K零压\p1_hold.dump":
        "10 K 保持段的过程轨迹。用于确认无序结构在整个保持期内持续存在，而非中途塌陷。",
    r"grain_boundary\experiments\phase1_2026-09\Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\5_保持10K零压\p1_held_0K.data":
        "终态在 0 K 下再做一次复核的结构，对应验收步骤 F（晶界 38.1%/38.0% BCC，内部 98.7% BCC）。",
    r"grain_boundary\experiments\phase1_2026-09\Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\6_结构分析\cna_phase1.py":
        "CNA 结构分析脚本，统计各阶段晶界与晶粒内部的 BCC 比例。",
    r"grain_boundary\experiments\phase1_2026-09\Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\6_结构分析\run_phase1.sh":
        "阶段一验收的完整编排脚本（依次跑五个步骤并收集结果）。",
    r"grain_boundary\experiments\phase1_2026-09\Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\总结报告.md":
        "阶段一零压运行的中文总结报告。",
    r"grain_boundary\experiments\phase1_2026-09\Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\RUN_META.md":
        "本次运行的元信息：时间、机器、参数、势函数标识。",
    r"grain_boundary\experiments\phase1_2026-09\Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\RUN_STATUS.md":
        "本次运行的实时执行记录与各步骤状态。",
    r"grain_boundary\experiments\phase1_2026-09\Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\README.md":
        "阶段一结果包的入口说明。",

    # --- 2×2 构建校验（主项目内） ---
    r"grain_boundary\Li_Sigma5_310_metastable_2x2\stage0_build_validate\Sigma5_310_K6_J6_L20.data":
        "2×2 模型的 Stage 0 构建产物（14400 原子），是验收包的前身。",
    r"grain_boundary\Li_Sigma5_310_metastable_2x2\stage0_build_validate\Sigma5_310_K6_J6_L20.cif":
        "2×2 模型的晶体学标准格式。",
    r"grain_boundary\Li_Sigma5_310_metastable_2x2\stage0_build_validate\Sigma5_310_K6_J6_L20.xyz":
        "2×2 模型的通用结构格式。",
    r"grain_boundary\Li_Sigma5_310_metastable_2x2\stage0_build_validate\removal_log.csv":
        "2×2 模型的重合位点删除记录。",
    r"grain_boundary\Li_Sigma5_310_metastable_2x2\stage0_build_validate\grain_map":
        "2×2 模型的晶粒归属图。",
    r"grain_boundary\Li_Sigma5_310_metastable_2x2\PARAMETERS.json":
        "2×2 扩展的参数记录：K 保持 6，J 由 3 加到 6，L 由 10 加到 20。",
}

BASENAME_PURPOSE = {
    "README.md": "仓库总入口：项目组成、核心结论、文档索引、环境与日常操作。",
    "bulk_Li.data": "体相参考结构的 LAMMPS 格式（128 原子 BCC）。",
    "bulk_Li.xyz": "体相参考结构的通用格式。",
    "bulk_Li.cif": "体相参考结构的晶体学格式。",
    "bulk_Li_minimized.data": "零压弛豫后的体相结构。由此得到该势函数下的平衡晶格常数 3.506186 Å，全部双晶都用它构建，避免强加体相应变。",
    "bulk_Li_minimized.restart": "体相弛豫的续算检查点。",
    "dump2data_validation.lammps": "dump→data 重建验证的运行日志，用于确认轨迹文件未丢信息。",
    "bulk_energy.lammps": "体相能量计算的运行日志，给出每原子势能参考值。",
    "kim.log": "OpenKIM 势函数调用日志，记录实际加载的模型标识。",
    "requirements.txt": "Python 依赖清单，用于重建运行环境。",
    "pytest.ini": "测试框架配置。",
    "conftest.py": "测试公共夹具。",
    "run_all.sh": "完整热历史编排脚本；带 .done 标记实现断点续算。",
    "resume.sh": "从中断处续算的脚本。",
    "run_phase1.sh": "阶段一验收的完整编排脚本。",
    "run_windows.ps1": "Windows 侧的辅助运行脚本。",
    "PARAMETERS.json": "2×2 扩展的参数记录（K=6, J=6, L=20）。",
    "params.json": "构建参数记录。",
    "parameters.json": "构建参数记录（副本）。",
    "removal_log.csv": "重合位点删除记录，保证删除对称可审计。",
    "structure.data": "初始结构的 LAMMPS 格式。",
    "structure.xyz": "初始结构的通用格式。",
    "grain_map": "晶粒归属图。",
    "cna_phase1.py": "CNA（共近邻分析）结构分析脚本，用于统计 BCC 分数。",
    "requirements_soap.txt": "SOAP 分析的 Python 依赖清单。",
    "setup_soap_env.sh": "SOAP 分析环境的一键配置脚本。",
    "verify_soap_env.py": "环境自检脚本，确认 SOAP 依赖可用。",
    "ni_stable.cfg": "输入构型之一：稳定的 Ni 晶界构型（源数据）。",
    "ni_metastable.cfg": "输入构型之一：亚稳的 Ni 晶界构型（源数据）。",
    "cfg_reader.py": "CFG 格式解析脚本。",
    "inspect_structures.py": "结构概览脚本，检查原子数与字段完整性。",
    "soap_test.py": "SOAP 初次可计算性测试脚本。",
    "stage1_qc.py": "Stage01 数据质量检查脚本。",
}

BASENAME_PREFIX_PURPOSE = [
    ("test_", "自动化回归测试：验证生成器几何、周期性、导出格式、工作流脚本等是否仍然正确。"),
    ("stage5b_matched_direction", "配对重采样中间表（体积最大的单项产物）。"),
    ("stage3a_", "Stage03A 阈值扫描产物：衡量 epsilon 取值对结果的影响。"),
    ("stage3b_", "Stage03B 类别构建产物：LAE 分类与 LER 人口指纹。"),
    ("stage4a_", "Stage04A 连续无序度产物。"),
    ("stage4b_", "Stage04B 无序度-能量关联产物。"),
    ("stage5a_", "Stage05A 物理描述符与空间分布产物。"),
    ("stage5b_", "Stage05B SOAP 方向相关产物。"),
    ("stage6_", "Stage06 最终评估产物：相关系数、模型对比与稳健性检验。"),
    ("stage2_75_", "Stage02.75 全人口 SOAP 与距离审计产物。"),
    ("stage2_5_", "Stage02.5 表示审计产物。"),
    ("stage2_", "Stage02 SOAP 表示产物。"),
    ("stage1_", "Stage01 数据准备与初测产物。"),
    ("TEAM_REVIEW_", "队友交付审查报告，含外部质疑与处理意见。"),
]

IN_PURPOSE = {
    "in.0_check": "Stage 0 输入脚本：核对模型几何、原子数与势函数是否被正确加载。",
    "in.1_minimize": "Stage 1 输入脚本：0 K 极小化，仅沿晶界法向做零压弛豫（y/z 固定在零压 CSL 周期），消除初始残余力。",
    "in.2_anneal": "Stage 2 输入脚本：363 K 等压（NPT）退火，是诱发晶界热重构的关键一步。",
    "in.2_anneal_resume": "Stage 2 的续算版本：从中断的检查点继续退火，不重跑已完成的时段。",
    "in.2_highT_npt": "2×2 实验的高温 NPT 输入脚本（对应 in.2_anneal 在旧版命名下的同一角色）。",
    "in.3_quench": "Stage 3 输入脚本：以 10 K/ps 从 363 K 淬火到 10 K，把高温无序结构冻结下来。",
    "in.3_quench_resume": "Stage 3 的续算版本。",
    "in.4_final": "Stage 4 输入脚本：最终弛豫与保持，产出用于结构分析的末态。",
    "in.4_10K_hold": "2×2 实验的 10 K 保持输入脚本（对应 in.4_final 在旧版命名下的同一角色）。",
    "in.bulk_energy": "体相能量参考计算：给出每原子势能基准，用于扣除体相贡献、得到晶界过剩量。",
    "in.dump2data": "dump→data 重建脚本：把轨迹转回结构文件，用来验证轨迹未丢失信息。",
    "in.test_li": "早期环境连通性测试：确认 LAMMPS + OpenKIM 能正确调用 Li 势函数。",
    "in.bicrystal_relax": "早期版本的双晶弛豫脚本（已被 v2 工作流取代）。",
    "in.quench": "早期版本的淬火脚本。",
    "in.resume": "早期版本的续算脚本。注意：旧版此文件曾直接复制主脚本，是已废弃流程中的一个已知缺陷。",
    "anneal_save.in": "早期退火脚本，含中途保存设置。",
    "continue.in": "早期续算脚本。",
    "dump2data.in": "早期 dump→data 转换脚本。",
}

EXT_PURPOSE = {
    ".data": "LAMMPS 原子构型：某个阶段结束时的静态结构快照。",
    ".dump": "LAMMPS 轨迹：该阶段的过程演化，多帧记录。",
    ".restart": "LAMMPS 续算检查点：保存速度场与热力学状态，只能在 LAMMPS 内读取。",
    ".xyz": "扩展 XYZ 结构文件，可携带每原子附加列。",
    ".cif": "晶体学标准结构文件。",
    ".grain_map": "晶粒归属图。",
    ".csv": "数据表。",
    ".npy": "NumPy 数组（中间计算量）。",
    ".npz": "NumPy 压缩数组（中间计算量）。",
    ".npz ": "NumPy 压缩数组。",
    ".png": "图表，可直接用于报告或论文。",
    ".txt": "纯文本结果或检查报告。",
    ".stats": "时间序列统计输出（温度/压力等）——结论的直接证据。",
    ".lammps": "LAMMPS 运行日志，含各阶段数值历史与断言结果。",
    ".screen": "进程标准输出记录。",
    ".done": "阶段完成标记，供编排脚本判断是否跳过。",
    ".log": "运行日志。",
    ".md": "Markdown 文档。",
    ".pdf": "PDF 文档。",
    ".html": "HTML 报告（含公式渲染）。",
    ".pptx": "汇报幻灯片。",
    ".py": "Python 脚本。",
    ".sh": "Shell 脚本（在 WSL/Linux 下运行）。",
    ".ps1": "PowerShell 脚本。",
    ".cfg": "原子构型输入文件。",
    ".json": "参数与元数据。",
    ".css": "网页样式表（随 HTML 报告）。",
    ".lmp": "LAMMPS 结构文件。",
}


def classify_layer(rel: str, ext: str, name: str) -> str:
    low = rel.lower()
    if name in ("ARTIFACT_LEDGER.md", "ARTIFACT_LEDGER.csv"):
        return "文档"
    if ext == ".pdf" or "05_文献" in rel or "_literature" in low:
        return "文献"
    if ext in (".py", ".sh", ".ps1", ".css"):
        return "代码"
    if (ext in (".cfg", ".cmd") or name == ".gitignore"
            or name.startswith("in.") or name.startswith("requirements")):
        return "输入"
    if ext == ".md" or ext in (".pptx", ".html"):
        return "文档"
    if ext in (".done", ".screen", ".log") or name in (
        "SHA256SUMS.txt", "ENVIRONMENT.txt", "params.json",
        "PARAMETERS.json", "parameters.json", "removal_log.csv",
    ):
        return "追溯"
    if ext in (".lammps",):
        return "追溯"
    if ext == ".restart":
        return "续算"
    if ext in (".dump", ".npy", ".npz"):
        return "中间"
    if ext == ".csv":
        size = 0
        try:
            size = os.path.getsize(ROOT / rel)
        except OSError:
            pass
        return "中间" if size > 5 * 1024 * 1024 else "结论"
    if ext in (".data", ".xyz", ".cif", ".grain_map", ".stats", ".txt", ".png", ".lmp"):
        return "结论"
    return "文档"


def describe(rel: str, ext: str, name: str) -> str:
    if rel in PATH_PURPOSE:
        return PATH_PURPOSE[rel]
    if name in BASENAME_PURPOSE:
        return BASENAME_PURPOSE[name]
    if name in IN_PURPOSE:
        return IN_PURPOSE[name]
    if name == ".gitignore":
        return "Git 忽略规则，指定哪些文件不纳入版本控制。"
    if "barofail" in name:
        return "柱压器失败的运行日志。保留了故障现场——本次运行修正了 Pdamp 参数，是问题定位的痕迹。"
    if name.endswith("_BI.md"):
        return "双语版报告（中文版与英文版首尾拼接在同一文件内）。"
    if name.startswith("README_BI"):
        return "双语版总览文档。"
    # Stage-specific shorthand only applies inside the SOAP project; the Li
    # project also has files named stage1_/stage2_ with unrelated meanings.
    if rel.startswith("soap\\"):
        for prefix, text in BASENAME_PREFIX_PURPOSE:
            if name.startswith(prefix):
                return text
    if name.startswith("p1_") and ext == ".lammps":
        return "该步骤的 LAMMPS 运行日志，含温度/压力/能量时间序列与断言检查。"
    if name.startswith("in.p1_"):
        return "该步骤的 LAMMPS 输入脚本。"
    if name.startswith("in."):
        return "LAMMPS 输入脚本。"
    if ext == ".lammps" and name.startswith("stage"):
        return "该阶段的 LAMMPS 运行日志，含数值历史与断言检查。"
    return EXT_PURPOSE.get(ext, "（未分类）")


def human(size: int) -> str:
    if size >= 1024 * 1024:
        return f"{size / 1024 / 1024:.2f} MB"
    if size >= 1024:
        return f"{size / 1024:.1f} KB"
    return f"{size} B"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    if not ROOT.is_dir():
        print(f"missing root: {ROOT}", file=sys.stderr)
        return 1

    rows = []
    for proj, proj_title in PROJECTS:
        for path in iter_project_files(proj):
            rel = str(path.relative_to(ROOT))
            name = path.name
            ext = path.suffix.lower()
            size = path.stat().st_size
            rows.append({
                "project": proj_title,
                "rel": rel,
                "name": name,
                "ext": ext or "(无)",
                "layer": classify_layer(rel, ext, name),
                "size": size,
                "purpose": describe(rel, ext, name),
                "sha256": sha256(path),
            })

    out_dir = OUT_DIR if OUT_DIR is not None else ROOT / "maintenance/records"
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = out_dir / "ARTIFACT_LEDGER.csv"
    with csv_path.open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["项目", "相对路径", "文件名", "扩展名", "层次", "大小(字节)", "用途说明", "SHA256"])
        for r in rows:
            w.writerow([r["project"], r["rel"], r["name"], r["ext"], r["layer"],
                        r["size"], r["purpose"], r["sha256"]])

    md = []
    md.append("# 研究产物台账")
    md.append("")
    md.append(f"扫描根目录：`{ROOT}`")
    md.append("")
    md.append(f"共 {len(rows)} 个文件，合计 {human(sum(r['size'] for r in rows))}。")
    md.append("")

    md.append("## 层次约定")
    md.append("")
    md.append("| 层次 | 含义 | 处置建议 |")
    md.append("|---|---|---|")
    advice = {
        "结论": "结果本身：末态结构、汇总数据、图、运行统计。",
        "续算": "检查点，保存速度场与热力学状态，无法仅凭结构文件还原。",
        "中间": "可从上游输入与脚本重算，代价是算力；可压缩存放或只留清单。",
        "追溯": "记录“结果怎么来的”，体积小、价值高，应随代码一起版本控制。",
        "代码": "生成器与分析脚本。",
        "输入": "构型文件与运行脚本。",
        "文献": "第三方论文，不纳入版本控制。",
        "文档": "报告与演示材料。",
    }
    layer_advice = {
        "结论": "必须保留",
        "续算": "保留至后续阶段完成",
        "中间": "可重算，建议只留清单+摘要",
        "追溯": "随代码进版本控制",
        "代码": "进版本控制",
        "输入": "进版本控制",
        "文献": "不进版本控制",
        "文档": "进版本控制",
    }
    for key, desc in LAYERS.items():
        md.append(f"| **{key}** | {advice[key]} | {layer_advice[key]} |")
    md.append("")

    for proj, proj_title in PROJECTS:
        proj_rows = [r for r in rows if r["project"] == proj_title]
        if not proj_rows:
            continue
        md.append(f"## {proj_title}")
        md.append("")
        md.append(f"`{proj}` — {len(proj_rows)} 个文件，{human(sum(r['size'] for r in proj_rows))}")
        md.append("")

        by_layer = defaultdict(list)
        for r in proj_rows:
            by_layer[r["layer"]].append(r)

        for key in LAYERS:
            group = by_layer.get(key)
            if not group:
                continue
            total = human(sum(r["size"] for r in group))
            md.append(f"### {key}（{len(group)} 个，{total}）")
            md.append("")
            md.append("| 文件 | 大小 | 是什么 / 有什么用 |")
            md.append("|---|---|---|")
            for r in sorted(group, key=lambda x: x["rel"]):
                shown = r["rel"].split("\\", 1)[1] if "\\" in r["rel"] else r["rel"]
                md.append(f"| `{shown}` | {human(r['size'])} | {r['purpose']} |")
            md.append("")

    md_path = out_dir / "ARTIFACT_LEDGER.md"
    md_path.write_text("\n".join(md), encoding="utf-8")

    print(f"files={len(rows)}")
    print(f"md={md_path}")
    print(f"csv={csv_path}")
    counts = defaultdict(int)
    sizes = defaultdict(int)
    for r in rows:
        counts[r["layer"]] += 1
        sizes[r["layer"]] += r["size"]
    for k in LAYERS:
        if counts[k]:
            print(f"  {k}: {counts[k]} files, {human(sizes[k])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
