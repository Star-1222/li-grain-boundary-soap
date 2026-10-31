#!/bin/bash
# ────────────────────────────────────────────────────────────────────────
#  Run script for bicrystal relaxation with LAMMPS + OpenKIM
# ────────────────────────────────────────────────────────────────────────
#  Usage:
#    chmod +x run_relax.sh
#    ./run_relax.sh
#
#  Environment variables:
#    LMP  — LAMMPS executable (default: lmp)
# ────────────────────────────────────────────────────────────────────────

LMP_EXE="${LMP:-lmp}"

# Copy the bicrystal data file into the working directory
# (Adjust the source path to match your output location)
DATA_SRC="../output/Sigma5_310.lmp"
DATA_DST="Sigma5_310.lmp"

if [ ! -f "$DATA_DST" ]; then
    if [ -f "$DATA_SRC" ]; then
        cp "$DATA_SRC" "$DATA_DST"
        echo "Copied data file: $DATA_SRC → $DATA_DST"
    else
        echo "Error: cannot find $DATA_SRC"
        echo "Run the bicrystal generator first, or copy the .lmp file manually."
        exit 1
    fi
fi

echo "Running LAMMPS relaxation..."
echo "  Executable: ${LMP_EXE}"
echo "  Input:      in.bicrystal_relax"
echo ""

${LMP_EXE} -in in.bicrystal_relax

if [ $? -eq 0 ]; then
    echo ""
    echo "Relaxation completed successfully."
    echo "Output files:"
    ls -lh Sigma5_310_*.data
else
    echo ""
    echo "LAMMPS run failed.  Check log.lammps for details."
    exit 1
fi