#!/usr/bin/env bash
# Resume an interrupted stage for the Li Sigma5 2x2 workflow.
# Usage: bash lammps/resume.sh stage1|stage2|stage3|stage4
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
STAGE="${1:?usage: resume.sh stage1|stage2|stage3|stage4}"

OUTDIR="${OUTDIR:-$ROOT}"
LOGDIR="${LOGDIR:-$ROOT/logs}"
KIM="${KIM:-SNAP_ZuoChenLi_2019_Li__MO_732106099012_000}"
SEED="${SEED:-12345}"
MPI_NP="${MPI_NP:-12}"
OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
MPIEXEC="${MPIEXEC:-mpirun}"
LMP="${LMP:-lmp}"

export OMP_NUM_THREADS
mkdir -p "$LOGDIR"

case "$STAGE" in
  stage1)
    SCRIPT="in.1_minimize"
    EXTRA="-var input_data $ROOT/stage0_build_validate/Sigma5_310_K6_J6_L20.data"
    ;;
  stage2)
    SCRIPT="in.2_highT_npt"
    EXTRA=""
    ;;
  stage3)
    SCRIPT="in.3_quench"
    EXTRA=""
    ;;
  stage4)
    SCRIPT="in.4_10K_hold"
    EXTRA=""
    ;;
  *)
    echo "unknown stage: $STAGE"
    exit 1
    ;;
esac

"$MPIEXEC" -np "$MPI_NP" "$LMP" \
  -in "$HERE/$SCRIPT" \
  -var outdir "$ROOT/$STAGE" \
  -var seed "$SEED" \
  -var kim_model "$KIM" \
  $EXTRA \
  -log "$LOGDIR/$STAGE.resume.lammps"

echo "[ok] $STAGE resume finished"
