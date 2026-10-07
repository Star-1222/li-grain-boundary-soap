"""Example: generate Sigma5 (310) [001] tilt GB in BCC Li.

Run from the project root:
    python examples/example.py
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from Li_Bicrystal_Generator.builder import BiCrystalBuilder, GrainBoundaryDescriptor
from Li_Bicrystal_Generator.config import get_default_config
from Li_Bicrystal_Generator.exporter import export


def main():
    config = get_default_config()
    config["lattice_constant"] = 3.51
    config["supercell"] = (20, 12, 12)
    config["overlap_cutoff"] = 2.3
    config["vacuum"] = 0.0
    config["periodic"] = (False, True, True)

    gb = GrainBoundaryDescriptor(
        name="Sigma5_310",
        sigma=5,
        axis=(0, 0, 1),
        angle=36.87,
        plane=(3, 1, 0),
        description="Sigma5 [001] tilt GB, (310) plane",
    )

    builder = BiCrystalBuilder(config)
    atoms = builder.build(gb)

    export(atoms, gb.name, output_dir="output", formats=("lammps", "xyz", "cif"))

    cell = atoms.cell.lengths()
    print(f"\nSigma5 (310) bicrystal generated: {len(atoms)} atoms")
    print(f"Cell: ({cell[0]:.3f}, {cell[1]:.3f}, {cell[2]:.3f}) A\n")
    print("Files written to output/")


if __name__ == "__main__":
    main()