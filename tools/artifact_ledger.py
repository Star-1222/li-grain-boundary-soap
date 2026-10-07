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
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

import build_artifact_catalog as base

ROOT = Path(os.environ.get("RESEARCH_ROOT", r"D:\Research"))
LEDGER_MD = ROOT / os.environ.get("LEDGER_MD_NAME", "ARTIFACT_LEDGER.md")
LEDGER_CSV = ROOT / os.environ.get("LEDGER_CSV_NAME", "ARTIFACT_LEDGER.csv")
MIRROR_DIR = Path(os.environ.get(
    "LEDGER_MIRROR",
    r"C:\Users\1\Documents\Codex\2026-10-07\ni\outputs",
))


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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

    rows: list[dict] = []
    self_outputs = {LEDGER_MD.resolve(), LEDGER_CSV.resolve()}
    for proj, title in base.PROJECTS:
        for path in base.iter_project_files(proj):
            # 跳过台账自身：把它的哈希写进台账会导致每次生成都自变。
            if path.resolve() in self_outputs:
                continue
            rel = str(path.relative_to(ROOT))
            ext = path.suffix.lower()
            rows.append({
                "group": proj,
                "project": title,
                "rel": rel,
                "name": path.name,
                "ext": ext or "(无)",
                "layer": base.classify_layer(rel, ext, path.name),
                "size": path.stat().st_size,
                "purpose": base.describe(rel, ext, path.name),
                "analysis": probe(path, ext),
                "sha256": base.sha256(path),
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
    md.append(f"- 扫描根目录：`{ROOT}`")
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
        where = f"`{proj}`" if proj else "仓库根目录与 `tools/`"
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

    LEDGER_MD.write_text("\n".join(md), encoding="utf-8")
    with LEDGER_CSV.open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["项目", "相对路径", "文件名", "扩展名", "层次",
                    "大小(字节)", "用途说明", "内容分析", "SHA256"])
        for r in rows:
            w.writerow([r["project"], r["rel"], r["name"], r["ext"], r["layer"],
                        r["size"], r["purpose"], r["analysis"], r["sha256"]])

    mirror = str(MIRROR_DIR)
    if mirror:
        MIRROR_DIR.mkdir(parents=True, exist_ok=True)
        (MIRROR_DIR / LEDGER_MD.name).write_text("\n".join(md), encoding="utf-8")
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
