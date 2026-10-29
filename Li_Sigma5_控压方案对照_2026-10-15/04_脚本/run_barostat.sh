#!/bin/bash
# =====================================================================
# Head-to-head barostat test at Stage 4 (10 K / 0 bar hold).
# Same start (stage3b_quenched.restart), same cell, same length.
#
#   xonly       : existing protocol - GB normal (x) only, Pdamp = 1000 ps
#   aniso       : x, y, z independent pistons,        Pdamp = 1000 ps
#   xonly_fast  : GB normal (x) only,                 Pdamp = 0.1 ps
#   aniso_fast  : x, y, z independent pistons,        Pdamp = 0.1 ps
#
# NOTE: in `units metal` Pdamp is in picoseconds, so the classic
# "... 0.0 0.0 1000.0" used by the Stage 2-4 inputs means 1000 ps, i.e. a
# piston three orders of magnitude too heavy for a 5 ps hold.
#
# The branches listed in TAGS run concurrently, NCPU MPI ranks each.
# Keep NCPU * (#branches) <= number of PHYSICAL cores; the box reports 24
# logical cores (hyperthreads) but only ~12 are physical.
# =====================================================================
set -euo pipefail
cd /home/lenovo/Li_meta_2x2

MODEL=SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
REST=/home/lenovo/Li_meta_2x2/stage3b_quench/stage3b_quenched.restart
IN=/home/lenovo/Li_meta_2x2/exp_barostat/lammps/in.baro_hold
NSTEPS=${NSTEPS:-5000}
NCPU=${NCPU:-6}
TAGS=${TAGS:-xonly_fast,aniso_fast}

if [ ! -f "$REST" ]; then
  echo "[fail] restart not found: $REST" >&2
  exit 1
fi

echo "[run ] start $(date)  nsteps=$NSTEPS  ncpu/branch=$NCPU  tags=$TAGS"

IFS=',' read -ra TLIST <<< "$TAGS"
for TAG in "${TLIST[@]}"; do
  OUT=/home/lenovo/Li_meta_2x2/exp_barostat/$TAG
  mkdir -p "$OUT"
  OMP_NUM_THREADS=1 mpirun -np "$NCPU" --oversubscribe /usr/local/bin/lmp-kim \
      -in "$IN" \
      -var outdir    "$OUT" \
      -var baro      "$TAG" \
      -var nsteps    "$NSTEPS" \
      -var restart   "$REST" \
      -var kim_model "$MODEL" \
      -log /home/lenovo/Li_meta_2x2/logs/baro_${TAG}.lammps \
      > /home/lenovo/Li_meta_2x2/logs/baro_${TAG}.out 2>&1 &
  echo "[spawn] $TAG (pid $!) -> $OUT"
done

wait
for TAG in "${TLIST[@]}"; do
  if grep -q "BARO_HOLD_DONE" "/home/lenovo/Li_meta_2x2/logs/baro_${TAG}.lammps"; then
    touch "/home/lenovo/Li_meta_2x2/logs/baro_${TAG}.done"
    echo "[done] $TAG $(date)"
  else
    echo "[WARN] $TAG did not finish cleanly" >&2
  fi
done
echo "[ok] ALL_BARO_DONE $(date)"
