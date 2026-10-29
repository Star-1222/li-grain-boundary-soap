#!/usr/bin/env python3
"""Figures for the barostat comparison at Stage 4 (10 K / 0 bar).

Inputs  : exp_barostat/analysis/thermo_<tag>.npz   (from compare_barostat.py)
          exp_barostat/analysis/<name>_profiles.npz (from analyze_structure.py)
Outputs : exp_barostat/fig_stress.png
          exp_barostat/fig_strain.png
          exp_barostat/fig_profile.png
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
ANA = ROOT / "analysis"

COL = {"xonly": "#c0392b", "aniso": "#2471a3",
       "xonly_fast": "#c0392b", "aniso_fast": "#2471a3"}
LBL = {"xonly": "x-only (Pdamp 1000 ps)", "aniso": "aniso (Pdamp 1000 ps)",
       "xonly_fast": "x-only (Pdamp 0.1 ps)", "aniso_fast": "aniso (Pdamp 0.1 ps)"}


def load(tag):
    p = ANA / f"thermo_{tag}.npz"
    return np.load(p) if p.is_file() else None


def main():
    tags = [t for t in ("xonly_fast", "aniso_fast", "xonly", "aniso")
            if (ANA / f"thermo_{t}.npz").is_file()]
    if not tags:
        raise SystemExit("[error] no thermo npz found; run compare_barostat.py first")

    # ------------------------------------------------ fig 1: stress -----
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2), sharex=True)
    for tag in tags:
        d = load(tag)
        t = d["Step"] * 1e-3                    # ps
        for ax, (k, lab) in zip(
                axes, [("Pxx", r"$P_{xx}$"), ("Pyy", r"$P_{yy}$"), ("Pzz", r"$P_{zz}$")]):
            ax.plot(t, d[k], color=COL[tag], lw=1.0,
                    label=LBL[tag] if k == "Pxx" else None)
            ax.axhline(0.0, color="grey", lw=0.6, ls=":")
    for ax, (k, lab) in zip(axes,
                            [("Pxx", r"$P_{xx}$"), ("Pyy", r"$P_{yy}$"), ("Pzz", r"$P_{zz}$")]):
        ax.set_title(lab, fontsize=13)
        ax.set_xlabel("time (ps)", fontsize=11)
        ax.set_ylabel("stress (bar)", fontsize=11)
        ax.grid(alpha=0.25)
    axes[0].legend(fontsize=9, loc="best")
    fig.suptitle("normal stress components during the 10 K / 0 bar hold "
                 "(target for every piston = 0 bar)", fontsize=12)
    fig.tight_layout()
    fig.savefig(ROOT / "fig_stress.png", dpi=200)
    plt.close(fig)

    # --------------------------------- fig 2: deviatoric + strain -------
    fig, axes = plt.subplots(1, 4, figsize=(17, 3.8), sharex=True)
    for tag in tags:
        d = load(tag)
        t = d["Step"] * 1e-3
        axes[0].plot(t, d["v_devP"], color=COL[tag], lw=1.2, label=LBL[tag])
        for ax, k, lab in zip(axes[1:], ["v_eXX", "v_eYY", "v_eZZ"],
                              [r"$\varepsilon_{xx}$", r"$\varepsilon_{yy}$",
                               r"$\varepsilon_{zz}$"]):
            ax.plot(t, 100.0 * d[k], color=COL[tag], lw=1.2)
    axes[0].set_title("deviatoric stress  $\\sqrt{\\frac{1}{2}\\sum_i(P_{ii}-\\bar P)^2}$",
                      fontsize=11)
    axes[0].set_ylabel("bar", fontsize=11)
    axes[0].legend(fontsize=9)
    for ax, lab in zip(axes[1:], [r"$\varepsilon_{xx}$", r"$\varepsilon_{yy}$",
                                  r"$\varepsilon_{zz}$"]):
        ax.axhline(0.0, color="grey", lw=0.6, ls=":")
        ax.set_title(lab, fontsize=13)
        ax.set_ylabel("%", fontsize=11)
    for ax in axes:
        ax.set_xlabel("time (ps)", fontsize=11)
        ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(ROOT / "fig_strain.png", dpi=200)
    plt.close(fig)

    # --------------------------------------------- fig 3: Q(x) ----------
    names = ["start_quenched", "final_xonly_fast", "final_aniso_fast"]
    available = [n for n in names if (ANA / f"{n}_profiles.npz").is_file()]
    if available:
        fig, axes = plt.subplots(1, 2, figsize=(13, 4.0))
        style = {"start_quenched": ("grey", "--", "start (quenched, shared)"),
                 "final_xonly_fast": (COL["xonly"], "-", "x-only, Pdamp 0.1 ps (final)"),
                 "final_aniso_fast": (COL["aniso"], "-", "aniso, Pdamp 0.1 ps (final)")}
        for n in available:
            z = np.load(ANA / f"{n}_profiles.npz")
            c, ls, lab = style.get(n, ("k", "-", n))
            for ax, key, ttl in zip(axes, ["qC", "qB"],
                                    ["coordination  CN(r<3.25 A) = 8",
                                     "bond-orientational  near-perfect-BCC in (q4,q6)"]):
                q = z[key]
                ax.plot(z["x"], q, color=c, ls=ls, lw=1.3, label=lab)
        for ax, ttl in zip(axes, ["coordination  CN(r < 3.25 A) = 8",
                                  "bond-orientational  near-perfect-BCC in (q4, q6)"]):
            ax.set_title(ttl, fontsize=11)
            ax.set_xlabel("x, GB normal (A)", fontsize=11)
            ax.set_ylabel("ordered fraction", fontsize=11)
            ax.set_ylim(-0.05, 1.05)
            ax.grid(alpha=0.25)
        axes[0].legend(fontsize=9)
        fig.tight_layout()
        fig.savefig(ROOT / "fig_profile.png", dpi=200)
        plt.close(fig)

    print("[ok] figures written to", ROOT)


if __name__ == "__main__":
    main()
