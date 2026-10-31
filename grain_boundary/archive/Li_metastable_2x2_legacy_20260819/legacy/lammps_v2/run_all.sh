#!/usr/bin/env bash
# Verified fresh Stage 0 -> 4 workflow.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"

OUTDIR="${OUTDIR:-$ROOT/data}"
LOGDIR="${LOGDIR:-$ROOT/logs}"
KIM="${KIM:-SNAP_ZuoChenLi_2019_Li__MO_732106099012_000}"
SEED="${SEED:-12345}"
MPI_NP="${MPI_NP:-12}"
OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
MPIEXEC="${MPIEXEC:-mpirun}"
LMP="${LMP:-lmp}"

export OMP_NUM_THREADS
mkdir -p "$OUTDIR" "$LOGDIR"

run_stage() {
  local stage="$1"
  local script="$2"
  local marker="$LOGDIR/$stage.done"

  if [[ -f "$marker" ]]; then
    echo "[skip] $stage"
    return 0
  fi

  echo "[run ] $stage with $MPI_NP MPI tasks"

  "$MPIEXEC" -np "$MPI_NP" "$LMP" \
    -in "$HERE/$script" \
    -var outdir "$OUTDIR" \
    -var seed "$SEED" \
    -var kim_model "$KIM" \
    -var input_data "$OUTDIR/Sigma5_310.data" \
    -var input_restart "$OUTDIR/stage3_quenched.restart" \
    -log "$LOGDIR/$stage.lammps"

  touch "$marker"
  echo "[done] $stage"
}

run_stage stage0_check stage0_test/in.0_check
run_stage stage1_minimize stage1_minimize/in.1_minimize
run_stage stage2_anneal stage2_anneal/in.2_anneal
run_stage stage3_quench stage3_quench/in.3_quench
run_stage stage4_final stage4_relax/in.4_final

echo "[ok] all stages completed"
