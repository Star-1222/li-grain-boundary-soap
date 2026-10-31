# -*- coding: utf-8 -*-
"""Compute grain-boundary excess energy from existing LAMMPS logs.

gamma_GB = (E_bicrystal - N * E_bulk_per_atom) / (2 * Ly * Lz)

The factor 2 counts the two crystallographically equivalent GBs that a
periodic bicrystal necessarily contains along the GB normal.

Reads thermo output by column NAME, because the header layout differs
between stages (some include Temp/KinEng, some include Lx/Ly/Lz).
"""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path


def _default_root() -> Path:
    """Repository root, derived from this file's location (maintenance/tools/../..)."""
    env = os.environ.get("RESEARCH_ROOT")
    if env:
        return Path(env)
    here = Path(__file__).resolve().parent
    return here.parents[1]


ROOT = _default_root()
PROD = (ROOT / "grain_boundary/history/construction_baseline_3600"
        / "Sigma5_310_Metastable_GB" / "v2_production")
ACC = (ROOT / "grain_boundary/history/experiments/phase1_2026-09"
       / "Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03")

# 1 eV/A^2 -> mJ/m^2
EV_A2_TO_MJ_M2 = 16021.76634
NUM = re.compile(r"-?[\d.]+(?:[eE][+-]?\d+)?")


def parse_thermo(path):
    """Return (header, last_row_as_dict, n_atoms)."""
    header = None
    last = None
    n_atoms = None
    try:
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if n_atoms is None:
                    m = re.match(r"^\s*(\d+)\s+atoms\s*$", line)
                    if m:
                        n_atoms = int(m.group(1))
                if re.match(r"^\s*Step\s", line):
                    header = line.split()
                    last = None
                    continue
                if header:
                    vals = line.split()
                    if len(vals) == len(header) and NUM.fullmatch(vals[0]):
                        last = vals
    except OSError:
        return None, None, None
    if not header or not last:
        return header, None, n_atoms
    return header, dict(zip(header, last)), n_atoms


def bulk_ref_per_atom(path):
    _header, row, n = parse_thermo(path)
    if not row or not n:
        return None
    return float(row["PotEng"]) / n, n


def gamma(pot_eng, n, e_bulk, ly, lz):
    excess = pot_eng - n * e_bulk
    area = 2.0 * ly * lz
    g = excess / area
    return excess, area, g, g * EV_A2_TO_MJ_M2


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

    bulk_log = PROD / "logs" / "bulk_energy.lammps"
    ref = bulk_ref_per_atom(bulk_log)
    if not ref:
        print("无法从 %s 读取体相参考能量" % bulk_log)
        return 1
    e_bulk, n_bulk = ref
    print("体相参考: %s" % bulk_log.name)
    print("  %d 原子, PotEng/n = %.10f eV/atom" % (n_bulk, e_bulk))
    print()

    # (标签, 日志, 是否为 0 K 弛豫结构)
    # 只有 0 K 弛豫结构的势能才能直接当作晶界能；有限温度快照的势能
    # 含热贡献，比出来的是"能量"而不是"晶界过剩能"。
    targets = [
        ("理想 CSL 晶界", PROD / "logs" / "stage1_minimize_checkpoint_mpi12.lammps", True),
        ("退火末帧 (300 K NVT)", PROD / "logs" / "stage2_anneal_mpi12.lammps", False),
        ("淬火末帧 (10 K)", PROD / "logs" / "stage3_quench_mpi12.lammps", False),
        ("最终弛豫 (亚稳无序)", PROD / "logs" / "stage4_final_mpi12.lammps", True),
    ]

    # 面内尺寸全流程不变（只有 x 方向弛豫），取 stage1 的 Ly/Lz。
    _h, s1, n1 = parse_thermo(targets[0][1])
    if not s1:
        print("无法读取 stage1 日志以确定面内尺寸")
        return 1
    ly, lz = float(s1["Ly"]), float(s1["Lz"])
    print("单条晶界面积 = Ly x Lz = %.6f x %.6f = %.3f A^2" % (ly, lz, ly * lz))
    print("盒子含两条等价晶界，总面积 = 2 x Ly x Lz")
    print()
    print("%-28s %6s %12s %9s %10s %9s %9s %s" % (
        "状态", "N", "PotEng(eV)", "dE(eV)", "g(eV/A2)", "g(mJ/m2)",
        "Press(bar)", "可比性"))
    print("-" * 112)

    results = []
    for label, log, is_zero_k in targets:
        _h, row, n = parse_thermo(log)
        if not row:
            print("%-30s 解析失败" % label)
            continue
        n_use = n or n1
        pe = float(row["PotEng"])
        excess, _area, g, g_mj = gamma(pe, n_use, e_bulk, ly, lz)
        if is_zero_k:
            results.append((label, g_mj))
        press = float(row.get("Press", "nan"))
        marks = []
        marks.append("0K" if is_zero_k else "含热贡献")
        if press == press and abs(press) > 50:
            marks.append("残余应力")
        mark = ", ".join(marks)
        print("%-28s %6d %12.4f %9.4f %10.6f %9.1f %9.1f %s" % (
            label, n_use, pe, excess, g, g_mj, press, mark))

    print()
    if len(results) >= 2:
        ideal = results[0][1]
        final = results[-1][1]
        print("【仅比较 0 K 弛豫结构】")
        print("理想 CSL 晶界  gamma = %.1f mJ/m^2" % ideal)
        print("亚稳无序晶界    gamma = %.1f mJ/m^2" % final)
        print("差值            dgamma = %+.1f mJ/m^2" % (final - ideal))
        if final < ideal:
            print("-> 热重构后能量更低：理想 CSL 结构并非该体系基态。")
        else:
            print("-> 热重构使能量升高，该差值即亚稳态的能量代价。")

    print()
    acc_logs = [p for p in sorted((ACC / "7_运行日志").glob("p1_*.lammps"))
                if "barofail" not in p.name]
    if acc_logs:
        print("2x2 验收包 (14400 原子，面积不同，仅供对照):")
        for p in acc_logs:
            _h, row, n = parse_thermo(p)
            if not row or "PotEng" not in row:
                continue
            if "Ly" in row and "Lz" in row:
                ly2, lz2 = float(row["Ly"]), float(row["Lz"])
                g_mj = gamma(float(row["PotEng"]), n or 14400,
                             e_bulk, ly2, lz2)[3]
                print("  %-22s N=%s  gamma = %8.1f mJ/m2  (Ly x Lz = %.2f x %.2f)"
                      % (p.stem, n, g_mj, ly2, lz2))
            else:
                print("  %-22s N=%s  PotEng = %.4f eV  (日志无盒子尺寸)"
                      % (p.stem, n, float(row["PotEng"])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
