# Bulk reference energy (same SNAP potential)

A perfect BCC Li single crystal is used as the bulk reference so the GB excess
energy can be computed with the same OpenKIM potential
`SNAP_ZuoChenLi_2019_Li__MO_732106099012_000`.

## Build the reference cell

```bash
cd /mnt/c/Users/1/Documents/Codex/2026-08-01/bcc-li-5-001-sigma5-310/v2
python -m generator.bulk --a 3.51 --n 4 --outdir outputs/bulk
# -> outputs/bulk/bulk_Li.data (128 atoms, 14.04 A cell), .xyz, .cif
```

## Compute the reference energy (LAMMPS)

```bash
lmp -in lammps/in.bulk_energy -var outdir ../outputs \
    -var kim_model SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
```

The script relaxes all three box directions at P = 0 (a bulk crystal has no GB
normal that must be held fixed), minimizes with the same force criterion
(`minimize 0.0 1.0e-6 10000 100000`), and prints `pe_per_atom`, `fmax`, `fnorm`
plus the full pressure tensor.  Outputs: `bulk_Li_minimized.data` and
`bulk_Li_minimized.restart`.

## GB excess energy (formula)

```
E_gb = (E_bicrystal - N_bicrystal * E_bulk_per_atom) / (2 * Ly * Lz)
     # two equivalent GBs, each of area Ly*Lz
```

The 2 in the denominator accounts for the two periodic GBs in the box.
