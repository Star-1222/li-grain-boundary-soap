#!/usr/bin/env python3
"""Command-line entry point for bicrystal generation.

Usage
-----
    python main.py --name Sigma5_310
    python main.py --name Sigma13_510 --lattice 3.52 --output-dir my_output
    python main.py --list

All physical parameters can be overridden via command-line flags.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from .builder import BiCrystalBuilder, GrainBoundaryDescriptor
from .config import (
    get_default_config,
    LATTICE_CONSTANT,
    OVERLAP_CUTOFF,
    OUTPUT_DIR,
    PERIODIC,
    SUPERCELL,
    VACUUM,
    GRAIN_BOUNDARIES,
)
from .exporter import SUPPORTED_FORMATS, export
from .utils import list_available_gbs, resolve_gb


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate BCC lithium bicrystal models for LAMMPS.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python -m Li_Bicrystal_Generator --name Sigma5_310\n"
            "  python -m Li_Bicrystal_Generator --name Sigma13_510 --lattice 3.52\n"
            "  python -m Li_Bicrystal_Generator --list\n"
        ),
    )

    # Required (or use --list)
    parser.add_argument(
        "--name", "-n", type=str, default=None,
        help="Grain boundary name from the database (use --list to see all).",
    )
    parser.add_argument(
        "--list", "-l", action="store_true",
        help="List all available grain boundaries and exit.",
    )

    # Optional overrides
    parser.add_argument(
        "--lattice", type=float, default=None,
        help=f"Lattice constant in Angstrom (default: {LATTICE_CONSTANT}).",
    )
    parser.add_argument(
        "--supercell", type=int, nargs=3, default=None,
        metavar=("NX", "NY", "NZ"),
        help=f"Supercell dimensions (default: {SUPERCELL[0]} {SUPERCELL[1]} {SUPERCELL[2]}).",
    )
    parser.add_argument(
        "--overlap-cutoff", type=float, default=None,
        help=f"Overlap removal cutoff in Angstrom (default: {OVERLAP_CUTOFF}).",
    )
    parser.add_argument(
        "--vacuum", type=float, default=None,
        help=f"Vacuum padding in Angstrom (default: {VACUUM}).",
    )
    parser.add_argument(
        "--periodic", type=str, nargs=3, default=None,
        metavar=("X", "Y", "Z"),
        help="Periodicity flags: True/False (default: False True True).",
    )

    # Output
    parser.add_argument(
        "--output-dir", "-o", type=str, default=None,
        help=f"Output directory (default: {OUTPUT_DIR}).",
    )
    parser.add_argument(
        "--formats", type=str, nargs="+", default=None,
        choices=SUPPORTED_FORMATS,
        help=f"Output formats (default: all).",
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Enable verbose logging.",
    )

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    # Logging
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=level, format="%(message)s")

    # --list mode
    if args.list:
        list_available_gbs()
        return

    if args.name is None:
        parser.print_help()
        sys.exit(1)

    # Resolve GB
    gb_dict = resolve_gb(args.name)

    # Build config
    config = get_default_config()
    if args.lattice is not None:
        config["lattice_constant"] = args.lattice
    if args.supercell is not None:
        config["supercell"] = tuple(args.supercell)
    if args.overlap_cutoff is not None:
        config["overlap_cutoff"] = args.overlap_cutoff
    if args.vacuum is not None:
        config["vacuum"] = args.vacuum
    if args.periodic is not None:
        config["periodic"] = tuple(
            p.lower() in ("true", "1", "yes") for p in args.periodic
        )
    output_dir = args.output_dir or config.get("output_dir", "output")
    formats = args.formats or config.get("output_formats", ("lammps", "xyz", "cif"))

    # Build
    builder = BiCrystalBuilder(config)
    gb_desc = GrainBoundaryDescriptor.from_dict(gb_dict)
    atoms = builder.build(gb_desc)

    # Export
    export(atoms, gb_desc.name, output_dir=output_dir, formats=formats)

    # Summary
    cell = atoms.cell.lengths()
    print(f"\nDone.  Bicrystal: {len(atoms)} atoms, "
          f"cell = ({cell[0]:.3f}, {cell[1]:.3f}, {cell[2]:.3f}) A\n")


if __name__ == "__main__":
    main()