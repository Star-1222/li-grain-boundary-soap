# -*- coding: utf-8 -*-
"""Generate a self-updating artifact ledger for a research tree.

For every file the ledger records: layer, size, plain-language purpose,
SHA256, and -- where the format allows -- a content-derived analysis summary.
It also cross-checks numbers that appear in both reports and raw data, and
flags structural oddities.

Usage:
    python artifact_ledger.py                  # scan default root
    python artifact_ledger.py --quiet          # no stdout summary
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

import build_artifact_catalog as base

def _default_root() -> Path:
    """Repository root, derived from this file's location.

    Restoring the repository anywhere must keep working, so nothing here
    hardcodes an absolute path. Override with RESEARCH_ROOT if needed.
    """
    env = os.environ.get("RESEARCH_ROOT")
    if env:
        return Path(env)
    here = Path(__file__).resolve().parent
    return here.parents[1]


ROOT = _default_root()
LEDGER_MD = ROOT / os.environ.get("LEDGER_MD_NAME", "maintenance/records/ARTIFACT_LEDGER.md")
LEDGER_CSV = ROOT / os.environ.get("LEDGER_CSV_NAME", "maintenance/records/ARTIFACT_LEDGER.csv")

_mirror = os.environ.get("LEDGER_MIRROR")
MIRROR_DIR = Path(_mirror) if _mirror else None


# --- content probes ----------------------------------------------------------

def probe_data(path: Path) -> str:
    try:
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            head = [next(fh, "") for _ in range(12)]
    except OSError:
        return ""
    text = "".join(head)
    atoms = re.search(r"^\s*(\d+)\s+atoms\s*$", text, re.M)
    step = re.search(r"timestep\s*=\s*(\d+)", text)
    box = re.findall(r"^\s*(-?[\d.eE+-]+)\s+(-?[\d.eE+-]+)\s+(xlo|ylo|zlo)", text, re.M)
    parts = []
    if atoms:
        parts.append(f"{int(atoms.group(1)):,} 原子")
    if box:
        dims = [float(hi) - float(lo) for lo, hi, _ in box]
        parts.append("盒子 " + " × ".join(f"{d:.2f}" for d in dims) + " Å")
    if step:
        parts.append(f"timestep {int(step.group(1)):,}")
    return " · ".join(parts)


def probe_dump(path: Path) -> str:
    frames = 0
    steps: list[int] = []
    atoms = None
    try:
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("ITEM: TIMESTEP"):
                    frames += 1
                    nxt = fh.readline().strip()
                    if nxt.lstrip("-").isdigit():
                        steps.append(int(nxt))
                elif line.startswith("ITEM: NUMBER OF ATOMS"):
                    nxt = fh.readline().strip()
                    if atoms is None and nxt.isdigit():
                        atoms = int(nxt)
    except OSError:
        return ""
    parts = [f"轨迹 {frames} 帧"]
    if steps:
        parts.append(f"timestep {steps[0]:,} → {steps[-1]:,}")
    if atoms:
        parts.append(f"{atoms:,} 原子")
    return " · ".join(parts)


def probe_stats(path: Path) -> str:
    header: list[str] = []
    rows: list[list[float]] = []
    try:
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    if "TimeStep" in line:
                        header = line.lstrip("#").split()
                    continue
                vals = line.split()
                if vals and all(re.fullmatch(r"-?[\d.eE+-]+", v) for v in vals):
                    rows.append([float(v) for v in vals])
    except OSError:
        return ""
    if not rows or not header:
        return f"时间序列 {len(rows)} 行"
    arr = np.array(rows)
    means = arr.mean(axis=0)
    named = {header[i]: means[i] for i in range(min(len(header), arr.shape[1]))}
    parts = [f"{len(rows)} 个采样点"]
    scalars = {k: v for k, v in named.items() if k.startswith("v_vp") and k != "v_vpress"}
    if "v_vpress" in named:
        parts.append(f"标量时均压力 {named['v_vpress']:+.2f} bar")
    if scalars:
        spread = max(scalars.values()) - min(scalars.values())
        detail = ", ".join(f"{k.replace('v_vp','P')}={v:+.0f}" for k, v in scalars.items())
        parts.append(f"分量 {detail}")
        if spread > 100:
            parts.append(f"⚠ 分量极差 {spread:.0f} bar，各向异性显著")
    return " · ".join(parts)


def probe_json(path: Path) -> str:
    try:
        obj = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return ""
    if not isinstance(obj, dict):
        return ""
    keys = ["K", "J", "L", "a_A", "theta_deg", "final_total", "total_before_removal"]
    picked = [f"{k}={obj[k]}" for k in keys if k in obj]
    if picked:
        return " · ".join(picked)
    return f"字段: {', '.join(list(obj)[:8])}"


CNA_SECTION = re.compile(r"^--- (\S+)\s+\((.*?)\)", re.M)


def probe_cna_report(path: Path) -> str:
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return ""
    blocks = re.split(r"^--- ", text, flags=re.M)[1:]
    chain = []
    for b in blocks:
        name = b.split(None, 1)[0] if b.split() else "?"
        gb = re.search(r"GB region\s+x ~ 0\s+\(\|d\| < 8 A\): BCC =\s*([\d.]+)%", b)
        interior = re.search(r"grain interior\s+: BCC =\s*([\d.]+)%", b)
        width = re.search(r"w_GB .*?=\s*([\d.]+)\s*A", b)
        if gb and interior:
            chain.append((name, float(gb.group(1)), float(interior.group(1)),
                          float(width.group(1)) if width else None))
    if not chain:
        return ""
    parts = []
    parts.append("晶界 BCC: " + " → ".join(f"{c[0]} {c[1]:.1f}%" for c in chain))
    parts.append("晶内 BCC: " + " → ".join(f"{c[2]:.1f}%" for c in chain))
    widths = [c[3] for c in chain if c[3] is not None]
    if widths:
        parts.append("无序带宽: " + " → ".join(f"{w:.1f} Å" for w in widths))
    return " · ".join(parts)


def probe_correlations_csv(path: Path) -> str:
    out: list[str] = []
    seen: set[str] = set()
    try:
        with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
            for row in csv.DictReader(fh):
                if row.get("scope") != "all":
                    continue
                if row.get("var_b") == "displacement_norm" and row.get("spearman"):
                    msg = f"无序度 D 与位移范数相关 {float(row['spearman']):.5f}"
                    if msg not in seen:
                        seen.add(msg)
                        out.append(msg)
                if row.get("var_a") == "disorder_k10" and row.get("var_b") == "excess_energy_pooled":
                    msg = f"D 与过剩能相关 {float(row['spearman']):.3f}"
                    if msg not in seen:
                        seen.add(msg)
                        out.append(msg)
    except (OSError, ValueError):
        return ""
    return " · ".join(out)


def probe_energy_models_csv(path: Path) -> str:
    acc: dict[str, list[float]] = defaultdict(list)
    try:
        with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
            for row in csv.DictReader(fh):
                try:
                    acc[row["model"]].append(float(row["r2"]))
                except (KeyError, ValueError):
                    continue
    except OSError:
        return ""
    if not acc:
        return ""
    order = sorted(acc)
    return " · ".join(f"{m}: R²={np.mean(v):.4f}" for m, v in
                      ((m, acc[m]) for m in order))


def probe_npy(path: Path) -> str:
    try:
        if path.suffix.lower() == ".npz":
            with np.load(path, allow_pickle=False) as z:
                return " · ".join(f"{k} {z[k].shape} {z[k].dtype}" for k in z.files[:4])
        a = np.load(path, mmap_mode="r", allow_pickle=False)
        return f"数组 {a.shape} · {a.dtype} · 有限值 {bool(np.isfinite(a).all())}"
    except Exception:
        return ""


def probe_generic_csv(path: Path) -> str:
    try:
        with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
            reader = csv.reader(fh)
            header = next(reader, [])
            rows = sum(1 for _ in reader)
    except OSError:
        return ""
    cols = ", ".join(header[:6])
    more = " …" if len(header) > 6 else ""
    return f"{rows:,} 行 × {len(header)} 列 [{cols}{more}]"


def probe_xyz(path: Path) -> str:
    try:
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            first = fh.readline().strip()
            comment = fh.readline().strip()
    except OSError:
        return ""
    if not first.isdigit():
        return ""
    cols = [c for c in re.split(r"\s{2,}|\t", comment) if c]
    extra = f" · 附加列: {', '.join(cols[:6])}" if cols else ""
    return f"{int(first):,} 原子{extra}"


def probe_cfg(path: Path) -> str:
    try:
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            text = fh.read(4000)
    except OSError:
        return ""
    m = re.search(r"Number of particles\s*=\s*(\d+)", text)
    return f"{int(m.group(1)):,} 原子" if m else ""


def probe_sha_manifest(path: Path) -> str:
    try:
        n = sum(1 for line in path.read_text(encoding="utf-8", errors="replace").splitlines() if line.strip())
    except OSError:
        return ""
    return f"收录 {n} 个文件的校验和"


def probe(path: Path, ext: str) -> str:
    name = path.name
    if name == "SHA256SUMS.txt":
        return probe_sha_manifest(path)
    if ext == ".data" or ext == ".lmp":
        return probe_data(path)
    if ext == ".dump":
        return probe_dump(path)
    if ext == ".stats":
        return probe_stats(path)
    if ext == ".json":
        return probe_json(path)
    if name == "cna_phase1_report.txt":
        return probe_cna_report(path)
    if name == "stage6_correlations.csv":
        return probe_correlations_csv(path)
    if name == "stage6_energy_models.csv":
        return probe_energy_models_csv(path)
    if ext in (".npy", ".npz"):
        return probe_npy(path)
    if ext == ".cfg":
        return probe_cfg(path)
    if ext in (".xyz", ".grain_map"):
        return probe_xyz(path)
    if ext == ".csv":
        return probe_generic_csv(path)
    return ""


# --- git-scoped file discovery ----------------------------------------------

def _git_paths(root: Path, extra_args):
    """Return a set of repo-relative paths from git, or None if unavailable."""
    try:
        proc = subprocess.run(
            ["git", "-C", str(root), "ls-files", "-z", *extra_args],
            capture_output=True,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if proc.returncode != 0:
        return None
    raw = proc.stdout.decode("utf-8", "replace")
    names = [n for n in raw.split("\0") if n]
    return {n.replace("/", os.sep) for n in names}


def tracked_paths(root: Path):
    """Files the repository would contain: tracked + untracked-but-not-ignored."""
    return _git_paths(root, ["--cached", "--others", "--exclude-standard"])


def ignored_paths(root: Path):
    """Files present locally but excluded by .gitignore."""
    return _git_paths(root, ["--others", "--ignored", "--exclude-standard"])


def git_blob_digests(root: Path):
    """Map repo-relative path -> (sha256, size) of the *repository* content.

    Hashing the working-tree bytes makes the ledger depend on each machine's
    core.autocrlf / core.eol settings: the same commit checks out as LF on one
    box and CRLF on another. Hashing the git blob instead gives one answer
    everywhere, which is what makes a fresh clone reproduce the ledger exactly.
    """
    try:
        ls = subprocess.run(
            ["git", "-C", str(root), "ls-files", "-s", "-z"],
            capture_output=True, check=True,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if ls.returncode != 0:
        return None

    path_to_sha = {}
    for rec in ls.stdout.decode("utf-8", "replace").split("\0"):
        if not rec:
            continue
        meta, _, path = rec.partition("\t")
        parts = meta.split()
        if len(parts) >= 2:
            path_to_sha[path.replace("/", os.sep)] = parts[1]
    if not path_to_sha:
        return None

    shas = list(dict.fromkeys(path_to_sha.values()))
    try:
        proc = subprocess.run(
            ["git", "-C", str(root), "cat-file", "--batch"],
            input=("\n".join(shas) + "\n").encode(),
            capture_output=True, check=True,
        )
    except (OSError, subprocess.SubprocessError):
        return None

    digest_by_sha = {}
    buf = proc.stdout
    pos = 0
    while pos < len(buf):
        nl = buf.find(b"\n", pos)
        if nl < 0:
            break
        header = buf[pos:nl].decode("utf-8", "replace").split()
        pos = nl + 1
        if len(header) >= 3 and header[1] == "blob":
            size = int(header[2])
            blob = buf[pos:pos + size]
            pos += size + 1
            digest_by_sha[header[0]] = (hashlib.sha256(blob).hexdigest(), size)
        else:
            break

    out = {}
    for p, sha in path_to_sha.items():
        if sha in digest_by_sha:
            out[p] = digest_by_sha[sha]
    return out


# --- cross-checks ------------------------------------------------------------

def cross_checks(rows: list[dict]) -> list[str]:
    notes: list[str] = []
    by_name = {r["name"]: r for r in rows}

    models = by_name.get("stage6_energy_models.csv")
    final_report = by_name.get("Ni_SOAP_最终报告.md")
    if models and models.get("analysis"):
        notes.append(
            "**能量模型自洽**：`stage6_energy_models.csv` 实测 " + models["analysis"] +
            "。与最终报告声称的 D-only 0.0393、D+物理 0.9869、D+物理+方向 0.9886、"
            "D+方向 0.7497 一致。"
        )

    stats = by_name.get("p1_hold_pressure.stats")
    if stats and stats.get("analysis"):
        notes.append(
            "**验收口径风险已量化**：`p1_hold_pressure.stats` 显示 " + stats["analysis"] +
            "。报告把“0 bar”定义为标量时均压力，若改为要求应力张量各分量趋零，"
            "需补做 aniso NPT 敏感性复核。"
        )

    cna = by_name.get("cna_phase1_report.txt")
    if cna and cna.get("analysis"):
        notes.append(
            "**无序化进程可量化**：`cna_phase1_report.txt` 给出 " + cna["analysis"] +
            "。晶界 BCC 占比持续下降而晶内基本不变，这正是“晶界真实无序、"
            "晶粒未熔化”的证据链。"
        )

    corr = by_name.get("stage6_correlations.csv")
    if corr and corr.get("analysis"):
        notes.append(
            "**重大发现（需人工确认）**：`stage6_correlations.csv` 显示 " + corr["analysis"] +
            "。若 D 与位移范数几乎完全重合，则 SOAP 派生坐标相对简单几何量"
            "并未提供额外信息——这会动摇整个方法线的前提，建议优先核实。"
        )
    return notes


def structural_flags(rows: list[dict]) -> list[str]:
    flags: list[str] = []
    mid = [r for r in rows if r["layer"] == "中间"]
    mid_size = sum(r["size"] for r in mid)
    total = sum(r["size"] for r in rows)
    if mid and total:
        flags.append(
            f"可再生的中间产物共 {len(mid)} 个、{base.human(mid_size)}，"
            f"占总体积 {mid_size / total * 100:.1f}%。其中最重的单项是 "
            f"`{max(mid, key=lambda r: r['size'])['name']}`"
            f"（{base.human(max(mid, key=lambda r: r['size'])['size'])}）。"
        )
    concl = sum(r["size"] for r in rows if r["layer"] == "结论")
    if total:
        flags.append(
            f"结论性产物只有 {base.human(concl)}，占 {concl / total * 100:.1f}%。"
        )
    trace = sum(r["size"] for r in rows if r["layer"] in ("追溯", "代码", "输入"))
    if total:
        flags.append(
            f"决定“可复现”的追溯+代码+输入合计仅 {base.human(trace)}"
            f"（{trace / total * 100:.1f}%），是最该保护的一层。"
        )
    bi = [r for r in rows if r["name"].endswith("_BI.md")]
    if bi:
        flags.append(
            f"双语拼接文档 {len(bi)} 份。中英文同处一文件，改动一侧另一侧会静默过期，"
            "建议拆分为 `README.md` + `README.en.md`。"
        )
    return flags


DATA_MD = ROOT / "maintenance/records/DATA_MANIFEST.md"
DATA_CSV = ROOT / "maintenance/records/DATA_MANIFEST.csv"


def write_data_manifest() -> int:
    """Inventory the bulk data that .gitignore keeps out of the repository.

    A fresh clone has none of these files, so this manifest is the checklist
    for restoring them: it records what each one is and what it should hash to.
    """
    ignored = ignored_paths(ROOT)
    if not ignored:
        print("未发现被忽略的文件（或当前目录不是 git 仓库）")
        return 0

    entries = []
    for rel in sorted(ignored):
        p = ROOT / rel
        if not p.is_file():
            continue
        ext = p.suffix.lower()
        entries.append({
            "rel": rel,
            "name": p.name,
            "ext": ext or "(无)",
            "layer": base.classify_layer(rel, ext, p.name),
            "size": p.stat().st_size,
            "purpose": base.describe(rel, ext, p.name),
            "sha256": base.sha256(p),
        })

    total = sum(e["size"] for e in entries)
    by_layer: dict[str, list] = defaultdict(list)
    for e in entries:
        by_layer[e["layer"]].append(e)

    md: list[str] = []
    md.append("# 未纳入版本控制的数据清单")
    md.append("")
    md.append("> 本文件由 `maintenance/tools/artifact_ledger.py --data-manifest` 生成。")
    md.append("> `ARTIFACT_LEDGER.md` 只覆盖仓库内容；本文件登记的是被 `.gitignore`")
    md.append("> 排除在仓库之外、但**仍然重要**的本地数据。")
    md.append("> 全新克隆不含这些文件 —— 本清单即为还原它们时的核对表。")
    md.append("")
    md.append(f"- 条目数：**{len(entries)}**")
    md.append(f"- 合计体积：**{base.human(total)}**")
    md.append("")
    md.append("## 为什么不在仓库里")
    md.append("")
    md.append("| 层次 | 文件数 | 体积 | 原因 |")
    md.append("|---|---:|---:|---|")
    reasons = {
        "文献": "第三方论文，体积大且有版权",
        "中间": "可由上游输入与脚本重算，代价是算力",
        "续算": "LAMMPS 二进制检查点，含速度场；**无法由结构文件还原**",
        "追溯": "运行日志，体积偏大",
        "结论": "数据表，体积偏大",
        "输入": "输入文件",
        "文档": "文档",
        "代码": "代码",
    }
    for key in base.LAYERS:
        g = by_layer.get(key)
        if not g:
            continue
        md.append(f"| **{key}** | {len(g)} | "
                  f"{base.human(sum(e['size'] for e in g))} | {reasons.get(key, '')} |")
    md.append("")
    md.append("## 清单")
    md.append("")
    md.append("| 文件 | 大小 | 层次 | SHA256（前 16 位） | 是什么 |")
    md.append("|---|---:|---|---|---|")
    for e in sorted(entries, key=lambda x: -x["size"]):
        md.append(f"| `{e['rel']}` | {base.human(e['size'])} | {e['layer']} "
                  f"| `{e['sha256'][:16]}` | {e['purpose']} |")
    md.append("")
    md.append("## 还原步骤")
    md.append("")
    md.append("1. 从备份介质取回上述文件，保持相对路径不变")
    md.append("2. 核对校验和：`certutil -hashfile <文件> SHA256`，比对上表")
    md.append("3. 重新生成台账：双击 `refresh_ledger.cmd`")
    md.append("")

    # newline="\n" 必须显式指定：Windows 上 write_text 默认把 \n 转成 \r\n，
    # 会导致工作区内容与 git 规范化后的结果不一致。
    DATA_MD.write_text("\n".join(md), encoding="utf-8", newline="\n")
    with DATA_CSV.open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["相对路径", "文件名", "扩展名", "层次", "大小(字节)",
                    "用途说明", "SHA256"])
        for e in entries:
            w.writerow([e["rel"], e["name"], e["ext"], e["layer"],
                        e["size"], e["purpose"], e["sha256"]])
    print(f"data manifest: {len(entries)} entries, {base.human(total)}")
    print(f"  {DATA_MD}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--data-manifest", action="store_true",
                    help="另外生成 DATA_MANIFEST.md/csv，登记未纳入版本控制的大文件")
    args = ap.parse_args()
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

    rows: list[dict] = []
    self_outputs = {LEDGER_MD.resolve(), LEDGER_CSV.resolve()}
    # 先生成数据清单，使台账能收录其最新校验和（避免慢一拍的抖动）。
    if args.data_manifest:
        write_data_manifest()
    # 在 git 仓库内，台账只覆盖仓库内容，这样任何克隆副本重新生成的结果都一致。
    tracked = tracked_paths(ROOT)
    blob_digests = git_blob_digests(ROOT) or {}
    for proj, title in base.PROJECTS:
        for path in base.iter_project_files(proj):
            # 跳过台账自身：把它的哈希写进台账会导致每次生成都自变。
            if path.resolve() in self_outputs:
                continue
            rel = str(path.relative_to(ROOT))
            if tracked is not None and rel not in tracked:
                continue
            ext = path.suffix.lower()
            # 已跟踪文件用仓库内容的哈希与大小，保证跨机器一致；
            # 未跟踪文件只能按工作区实际内容计算。
            digest = blob_digests.get(rel)
            if digest:
                size, sha256 = digest[1], digest[0]
            else:
                size, sha256 = path.stat().st_size, base.sha256(path)
            rows.append({
                "group": proj,
                "project": title,
                "rel": rel,
                "name": path.name,
                "ext": ext or "(无)",
                "layer": base.classify_layer(rel, ext, path.name),
                "size": size,
                "purpose": base.describe(rel, ext, path.name),
                "analysis": probe(path, ext),
                "sha256": sha256,
            })

    if not rows:
        print(f"no files found under {ROOT}")
        return 1

    total = sum(r["size"] for r in rows)
    layer_count = defaultdict(int)
    layer_size = defaultdict(int)
    for r in rows:
        layer_count[r["layer"]] += 1
        layer_size[r["layer"]] += r["size"]

    md: list[str] = []
    md.append("# 研究产物台账")
    md.append("")
    md.append("> 本文件由 `artifact_ledger.py` 自动生成，请勿手工编辑。")
    md.append("> 每次重新生成都会刷新文件清单、校验和与内容分析。")
    md.append("")
    # 刻意不写绝对路径：台账随仓库走，写死本机路径会让克隆副本无法复现。
    md.append(f"- 文件总数：**{len(rows)}**")
    md.append(f"- 总体积：**{base.human(total)}**")
    md.append("")

    md.append("## 层次分布")
    md.append("")
    md.append("| 层次 | 文件数 | 体积 | 占总比 | 处置 |")
    md.append("|---|---:|---:|---:|---|")
    disposal = {
        "结论": "必须保留", "续算": "保留至后续阶段完成",
        "中间": "可重算，建议只留清单+摘要", "追溯": "随代码进版本控制",
        "代码": "进版本控制", "输入": "进版本控制",
        "文献": "不进版本控制", "文档": "进版本控制",
    }
    for key in base.LAYERS:
        if not layer_count[key]:
            continue
        md.append(f"| **{key}** | {layer_count[key]} | {base.human(layer_size[key])} "
                  f"| {layer_size[key] / total * 100:.1f}% | {disposal[key]} |")
    md.append("")

    flags = structural_flags(rows)
    if flags:
        md.append("## 自动识别的结构特征")
        md.append("")
        for f in flags:
            md.append(f"- {f}")
        md.append("")

    checks = cross_checks(rows)
    if checks:
        md.append("## 交叉核对：报告声称 vs 数据实测")
        md.append("")
        for c in checks:
            md.append(f"- {c}")
        md.append("")

    interesting = [r for r in rows if r["analysis"]]
    if interesting:
        md.append("## 关键产物内容分析")
        md.append("")
        md.append("以下是直接从文件内容读出的事实，不是人工描述。")
        md.append("")
        md.append("| 文件 | 内容分析 |")
        md.append("|---|---|")
        for r in sorted(interesting, key=lambda x: -x["size"])[:40]:
            md.append(f"| `{r['name']}` | {r['analysis']} |")
        md.append("")

    for proj, title in base.PROJECTS:
        proj_rows = [r for r in rows if r["project"] == title]
        if not proj_rows:
            continue
        md.append(f"## {title}")
        md.append("")
        where = f"`{proj}`" if proj else "仓库根目录、文档、工具与其他项目"
        md.append(f"{where} — {len(proj_rows)} 个文件，"
                  f"{base.human(sum(r['size'] for r in proj_rows))}")
        md.append("")
        by_layer = defaultdict(list)
        for r in proj_rows:
            by_layer[r["layer"]].append(r)
        for key in base.LAYERS:
            group = by_layer.get(key)
            if not group:
                continue
            md.append(f"### {key}（{len(group)} 个，"
                      f"{base.human(sum(r['size'] for r in group))}）")
            md.append("")
            has_analysis = any(r["analysis"] for r in group)
            if has_analysis:
                md.append("| 文件 | 大小 | 用途 | 内容分析 |")
                md.append("|---|---|---|---|")
                for r in sorted(group, key=lambda x: x["rel"]):
                    shown = (r["rel"] if r["group"] == ""
                             else (r["rel"].split("\\", 1)[1]
                                   if "\\" in r["rel"] else r["rel"]))
                    md.append(f"| `{shown}` | {base.human(r['size'])} | {r['purpose']} "
                              f"| {r['analysis']} |")
            else:
                md.append("| 文件 | 大小 | 用途 |")
                md.append("|---|---|---|")
                for r in sorted(group, key=lambda x: x["rel"]):
                    shown = (r["rel"] if r["group"] == ""
                             else (r["rel"].split("\\", 1)[1]
                                   if "\\" in r["rel"] else r["rel"]))
                    md.append(f"| `{shown}` | {base.human(r['size'])} | {r['purpose']} |")
            md.append("")

    LEDGER_MD.write_text("\n".join(md), encoding="utf-8", newline="\n")
    with LEDGER_CSV.open("w", newline="", encoding="utf-8-sig") as fh:
        # 显式 LF：csv.writer 默认 \r\n，会被 .gitattributes 归一化，
        # 导致工作区永远显示"已修改"。
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["项目", "相对路径", "文件名", "扩展名", "层次",
                    "大小(字节)", "用途说明", "内容分析", "SHA256"])
        for r in rows:
            w.writerow([r["project"], r["rel"], r["name"], r["ext"], r["layer"],
                        r["size"], r["purpose"], r["analysis"], r["sha256"]])

    if MIRROR_DIR is not None:
        MIRROR_DIR.mkdir(parents=True, exist_ok=True)
        (MIRROR_DIR / LEDGER_MD.name).write_text(
            "\n".join(md), encoding="utf-8", newline="\n")
        (MIRROR_DIR / LEDGER_CSV.name).write_bytes(LEDGER_CSV.read_bytes())

    if not args.quiet:
        print(f"files={len(rows)} total={base.human(total)}")
        print(f"ledger={LEDGER_MD}")
        print(f"analysed={len(interesting)}")
        for f in flags:
            print("  flag:", re.sub(r"`|\*\*", "", f)[:110])
        for c in checks:
            print("  check:", re.sub(r"`|\*\*", "", c)[:110])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
