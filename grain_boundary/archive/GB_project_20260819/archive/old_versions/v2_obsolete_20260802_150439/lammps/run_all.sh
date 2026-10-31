#!/usr/bin/env bash
# v2 Sigma5 (310)[001] BCC Li -- run the full staged LAMMPS workflow.
#
# Checkpoint-safe: each stage writes logs/<stage>.done ONLY after it
# finishes; re-running this script skips finished stages and continues
# from the first unfinished one.  Every stage uses its own log file
# (logs/<stage>.lammps).  Interrupt with Ctrl-C; then re-run to resume.
#
# Usage:
#   bash run_all.sh                    # full workflow
#   OUTDIR=/path/to/outputs bash run_all.sh
#   SEED=42 bash run_all.sh            # different RNG seed
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"
OUTDIR="${OUTDIR:-$HERE/../outputs}"
LOGDIR="${LOGDIR:-$HERE/../logs}"
KIM="${KIM:-SNAP_ZuoChenLi_2019_Li__MO_732106099012_000}"
SEED="${SEED:-12345}"
LMP="${LMP:-lmp}"
mkdir -p "$OUTDIR" "$LOGDIR"

run_stage() {
  local stage="$1" script="$2"
  if [[ -f "$LOGDIR/$stage.done" ]]; then
    echo "[skip] $stage (already done)"
    return 0
  fi
  echo "[run ] $stage  ->  logs/$stage.lammps"
  "$LMP" -in "$script" -var outdir "$OUTDIR" -var seed "$SEED" -var kim_model "$KIM" \
         -log "$LOGDIR/$stage.lammps"
  touch "$LOGDIR/$stage.done"
  echo "[done] $stage"
}

run_stage stage0_check in.0_check
run_stage stage1_minimize in.1_minimize
run_stage stage2_anneal in.2_anneal
run_stage stage3_quench in.3_quench
run_stage stage4_final in.4_final

echo "[ok] all stages completed"
