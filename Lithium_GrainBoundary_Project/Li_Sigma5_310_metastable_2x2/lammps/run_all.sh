#!/usr/bin/env bash
# Li Sigma5 (310)[001] 2x2 -- Stage 0 -> 4 runner (requires LAMMPS + OpenKIM).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"

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

run_stage() {
  local stage="$1"
  local script="$2"
  local stage_outdir="$3"
  local marker="$LOGDIR/$stage.done"

  if [[ -f "$marker" ]]; then
    echo "[skip] $stage"
    return 0
  fi

  echo "[run ] $stage"
  "$MPIEXEC" -np "$MPI_NP" "$LMP" \
    -in "$HERE/$script" \
    -var outdir "$ROOT/$stage_outdir" \
    -var seed "$SEED" \
    -var kim_model "$KIM" \
    -var input_data "$ROOT/stage0_build_validate/Sigma5_310_K6_J6_L20.data" \
    -log "$LOGDIR/$stage.lammps"

  touch "$marker"
  echo "[done] $stage"
}

run_stage stage0_check in.0_check stage0_build_validate
run_stage stage1_relax in.1_minimize stage1_relax
run_stage stage2_highT_NPT in.2_highT_npt stage2_highT_NPT
run_stage stage3_quench in.3_quench stage3_quench
run_stage stage4_10K_hold in.4_10K_hold stage4_10K_hold

echo "[ok] all stages completed"
