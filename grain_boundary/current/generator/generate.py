"""CLI entry point: build, verify and export the Sigma5 (310)[001] bicrystal.

Usage:
    python -m generator.generate [--a 3.506186475312] [--k 6] [--j 6] [--l 20]
                                 [--outdir outputs] [--model-size full|small]

The full model (default) reproduces the reference geometry: 14880 sites,
480 coincident-atom removals on the two GB planes, 14400 final atoms.
The small model (--model-size small, K=3 J=3 L=5) is for CI/tests.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from .builder import build_bicrystal
from .config import (DEFAULT_A, DEFAULT_J, DEFAULT_K, DEFAULT_L, ELEMENT,
                     KIM_MODEL, MASS_LI, THETA_DEG)
from .exporters import (write_cif, write_extxyz, write_grain_id_map,
                        write_lammps_data, write_params_json,
                        write_removal_log_csv, write_xyz_with_grain_id)
from .verify import run_all_checks, VerificationError

PRESETS = {"full": (DEFAULT_K, DEFAULT_J, DEFAULT_L),
           "small": (3, 3, 5)}


def render_build_report(bc, outdir: Path, prefix: str) -> str:
    m = bc.meta
    Lx, Ly, Lz = bc.box
    rows = [
        f"# Build report: Sigma5 (310)[001] BCC Li bicrystal ({prefix})",
        "",
        "## Construction method",
        "",
        "The corrected generator uses the **mirror (CSL) construction**:",
        "",
        "1. Grain 1 is generated from its own rotated BCC lattice with orientation",
        f"   matrix M1 = Rz(+theta/2), theta = 2*atan(1/3) = {THETA_DEG:.8f} deg.",
        "   Conventional cubic BCC basis (0,0,0)+(1/2,1/2,1/2) is used; the",
        "   primitive-cell mistake of the legacy generator is not repeated.",
        "2. Grain 2 is the exact mirror image of grain 1 across x = Lx/2.",
        "   No rotation inside a fixed periodic box is ever performed.",
        "3. Overlap removal deletes ONLY cross-grain coincident pairs sitting on",
        "   the two GB planes (x=0 and x=Lx/2). Nothing else is deleted.",
        f"4. Removal: {m['removed_total']} atoms ({m['removal_fraction_pct']:.4f}%),",
        f"   g1={m['removed_grain1']}, g2={m['removed_grain2']},",
        f"   GB1={sum(1 for r in bc.removal_log if r['gb']=='GB1')},",
        f"   GB2={sum(1 for r in bc.removal_log if r['gb']=='GB2')}.",
        "",
        "## Parameters",
        "",
        f"- lattice constant a = {bc.a} A",
        f"- misorientation theta = {THETA_DEG:.8f} deg about [001] (Sigma5)",
        f"- GB plane (310); box periods K={bc.K}, J={bc.J}, L={bc.L}",
        f"- box: {Lx:.6f} x {Ly:.6f} x {Lz:.6f} A (orthogonal, p p p)",
        f"- in-plane periods: y = sqrt(10)*a = {math.sqrt(10)*bc.a:.6f} A, z = a = {bc.a} A",
        f"- atoms: {m['total_before_removal']} -> {m['final_total']} (per grain {m['final_grain1']}/{m['final_grain2']})",
        f"- density: {m['final_total']/(Lx*Ly*Lz):.12f} atoms/A^3 = 2/a^3",
        "",
        "## Orientation matrices (rows = crystal axes in the lab frame)",
        "",
        "M1 (grain 1) = Rz(+theta/2):",
        "    [[0.9486832981, -0.3162277660, 0],",
        "     [0.3162277660,  0.9486832981, 0],",
        "     [0, 0, 1]]",
        "M2 (grain 2) = Rz(-theta/2):",
        "    [[0.9486832981,  0.3162277660, 0],",
        "     [-0.3162277660, 0.9486832981, 0],",
        "     [0, 0, 1]]",
        f"M2 . M1^T = Rz(-theta), angle {THETA_DEG:.8f} deg about [001].",
        "",
        "## Files written",
        "",
        f"- `{outdir}/{prefix}.data`       LAMMPS data (orthogonal, Velocities=0)",
        f"- `{outdir}/{prefix}.xyz`        extended XYZ with grain_id",
        f"- `{outdir}/{prefix}.cif`        CIF (P1)",
        f"- `{outdir}/{prefix}_grain_id.xyz`  plain XYZ with grain_id column",
        f"- `{outdir}/{prefix}.grain_map`  atom_id -> grain_id map",
        f"- `{outdir}/removal_log.csv`     full record of deleted atoms",
        f"- `{outdir}/params.json`         construction parameters",
        f"- `{outdir}/verification_report.md`  all geometric checks",
        "",
        "## Next steps (LAMMPS)",
        "",
        "Run the staged workflow in v2/lammps/ (see v2/README.md).",
        f"Potential: OpenKIM {KIM_MODEL} (tentative).",
    ]
    return "\n".join(rows) + "\n"


def render_verification_report(res: dict) -> str:
    """Render the verification results as Markdown."""
    rows: list = []
    def h(level, text):
        rows.append("#" * level + " " + text)

    h(1, "Geometric verification report (v2 Sigma5 (310)[001] BCC Li)")
    rows.append("")
    p = res["parameters"]
    rows.append(f"- a = {p['a_A']} A, theta = {p['theta_deg']:.8f} deg, axis [001], plane {p['gb_plane']}")
    rows.append(f"- box = {p['box_A'][0]:.6f} x {p['box_A'][1]:.6f} x {p['box_A'][2]:.6f} A (K={p['K']} J={p['J']} L={p['L']})")
    rows.append(f"- in-plane periods y = {p['in_plane_periods_A'][0]:.6f} A, z = {p['in_plane_periods_A'][1]:.6f} A")
    rows.append("")
    h(2, "1. Misorientation and rotation axis")
    o = res["orientation"]
    rows.append(f"- angle = {o['misorientation_angle_deg']:.8f} deg (target 36.86989765)")
    rows.append(f"- axis = {o['rotation_axis']} (target [001])")
    rows.append(f"- cubic-symmetry minimum = {o['cubic_symmetry_min_angle_deg']:.8f} deg")
    rows.append("")
    h(2, "2. GB plane")
    gp = res["gb_plane"]
    rows.append(f"- normal in crystal 1 = {gp['gb_normal_crystal1']} ~ {gp['gb_normal_crystal1_miller']}")
    rows.append(f"- normal in crystal 2 = {gp['gb_normal_crystal2']} ~ {gp['gb_normal_crystal2_miller']}")
    rows.append("")
    h(2, "3. Two GBs, explicit and equivalent")
    g = res["two_gbs"]
    rows.append(f"- GB1 at x = {g['gb1_x']}, GB2 at x = {g['gb2_x']} A")
    rows.append(f"- cross-GB pairs < 4.5 A: GB1 = {g['cross_gb_pairs_below_4p5A_gb1']}, GB2 = {g['cross_gb_pairs_below_4p5A_gb2']}")
    rows.append(f"- spectra identical: {g['spectra_identical']}")
    rows.append(f"- min cross-GB distance: GB1 = {g['min_gb1_A']:.6f} A, GB2 = {g['min_gb2_A']:.6f} A")
    rows.append(f"- pairs at 2.2199 A: GB1 = {g['pairs_at_2p22A_gb1']}, GB2 = {g['pairs_at_2p22A_gb2']}")
    rows.append("")
    h(2, "4. Bulk BCC away from the GBs")
    b = res["bulk_bcc"]
    rows.append(f"- margin = {b['margin_A']} A; sampled = {b['sampled']}, perfect = {b['perfect']}, failed = {b['failed']}")
    rows.append("")
    h(2, "5. Distances")
    d = res["distances"]
    rows.append(f"- global min = {d['global_min_A']:.6f} A (analytic 2d = {d['analytic_2d_A']:.6f} A)")
    rows.append(f"- BCC NN = {d['nn_bcc_A']:.6f} A, NNN = {d['nnn_bcc_A']:.6f} A")
    rows.append("")
    h(2, "6. Voids")
    v = res["voids"]
    rows.append(
        f"- occupied x planes: {v['occupied_x_planes']}; "
        f"expected gap = {v['expected_plane_gap_A']} A; "
        f"gap range = [{v['min_plane_gap_A']}, {v['max_plane_gap_A']}] A; "
        f"missing planes = {v['missing_x_planes']}"
    )
    rows.append("")
    h(2, "7. Removal log")
    r = res["removal"]
    rows.append(f"- total removed = {r['removed_total']} (g1 = {r['removed_grain1']}, g2 = {r['removed_grain2']}; GB1 = {r['gb1_removed']}, GB2 = {r['gb2_removed']})")
    rows.append(f"- all deleted atoms on GB planes: {r['all_on_gb_plane']}")
    rows.append("")
    h(2, "8. Counts and density")
    c = res["counts"]
    rows.append(f"- atoms = {c['final_total']} (g1 = {c['final_grain1']}, g2 = {c['final_grain2']})")
    rows.append(f"- density = {c['density_atoms_per_A3']:.12f} atoms/A^3 (2/a^3 = {c['density_target_2_over_a3']:.12f})")
    rows.append("")
    h(2, "9. In-plane periodicity (y, z)")
    pp = res["periodicity"]
    rows.append(f"- lattice-index misses: grain1 = {pp['grain1_misses']}, grain2 = {pp['grain2_misses']}")
    rows.append(f"- y period = {pp['y_period_A']:.6f} A = sqrt(10)*a; z period = {pp['z_period_A']:.6f} A = a")
    rows.append("")
    h(2, "10. Legacy regression")
    lc = res["legacy_comparison"]
    rows.append(f"- this model: {lc['before']} -> {lc['after']} atoms, removal {lc['removal_pct']:.4f}%")
    rows.append("- legacy model: 3600 -> 3224, 10.44% (must NOT be reproduced)")
    rows.append(f"- legacy-bad scenario present: {lc['legacy_bad']}")
    rows.append("")
    rows.append(f"**Overall: PASSED = {res['passed']}**")
    rows.append("")
    return "\n".join(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Build + verify + export the Sigma5 GB")
    ap.add_argument("--a", type=float, default=DEFAULT_A)
    ap.add_argument("--k", type=int, default=None)
    ap.add_argument("--j", type=int, default=None)
    ap.add_argument("--l", type=int, default=None)
    ap.add_argument("--outdir", type=str, default="outputs")
    ap.add_argument("--model-size", choices=["full", "small"], default=None)
    args = ap.parse_args(argv)

    if args.model_size:
        K, J, L = PRESETS[args.model_size]
        prefix = "Sigma5_310_small" if args.model_size == "small" else "Sigma5_310"
    else:
        K = args.k or DEFAULT_K
        J = args.j or DEFAULT_J
        L = args.l or DEFAULT_L
        prefix = "Sigma5_310"
    if args.k or args.j or args.l:
        prefix = f"Sigma5_310_K{K}_J{J}_L{L}"

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    bc = build_bicrystal(a=args.a, K=K, J=J, L=L)
    try:
        res = run_all_checks(bc)
    except VerificationError as exc:
        print(f"[verify] FAILED: {exc}")
        return 2

    write_lammps_data(bc, outdir / f"{prefix}.data")
    write_extxyz(bc, outdir / f"{prefix}.xyz")
    write_cif(bc, outdir / f"{prefix}.cif")
    write_xyz_with_grain_id(bc, outdir / f"{prefix}_grain_id.xyz")
    write_grain_id_map(bc, outdir / f"{prefix}.grain_map")
    write_removal_log_csv(bc, outdir / "removal_log.csv")
    write_params_json(bc, outdir / "params.json")
    (outdir / "build_report.md").write_text(render_build_report(bc, outdir, prefix), encoding="utf-8")
    (outdir / "verification_report.md").write_text(render_verification_report(res), encoding="utf-8")

    print(f"[ok] built {prefix}: {len(bc.atoms)} atoms, box {bc.box[0]:.4f} x {bc.box[1]:.4f} x {bc.box[2]:.4f} A")
    print("[ok] all geometric checks passed; outputs written to", outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
