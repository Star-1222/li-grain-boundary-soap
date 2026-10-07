#!/usr/bin/env bash
# Resume one interrupted stage from its latest validated checkpoint.
# Usage: bash resume.sh stage1|stage2|stage3|stage4
set -euo pipefail

STAGE="${1:-}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"

OUTDIR="${OUTDIR:-$ROOT/data}"
LOGDIR="${LOGDIR:-$ROOT/logs}"
KIM="${KIM:-SNAP_ZuoChenLi_2019_Li__MO_732106099012_000}"
MPI_NP="${MPI_NP:-12}"
OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
MPIEXEC="${MPIEXEC:-mpirun}"
LMP="${LMP:-lmp}"

export OMP_NUM_THREADS
mkdir -p "$OUTDIR" "$LOGDIR"

case "$STAGE" in
  stage1)
    INPUT="$OUTDIR/stage1_checkpoint.data"
    [[ -f "$INPUT" ]] || {
      echo "Missing $INPUT"
      exit 2
    }
    SCRIPT="$HERE/stage1_minimize/in.1_minimize"
    EXTRA=(-var input_data "$INPUT")
    MARKER="stage1_minimize.done"
    ;;

  stage2)
    compgen -G "$OUTDIR/stage2_anneal.checkpoint.*" >/dev/null || {
      echo "No Stage 2 checkpoint found"
      exit 2
    }
    SCRIPT="$HERE/stage2_anneal/in.2_anneal_resume"
    EXTRA=(-var restart "$OUTDIR/stage2_anneal.checkpoint.*")
    MARKER="stage2_anneal.done"
    ;;

  stage3)
    compgen -G "$OUTDIR/stage3_quench.checkpoint.*" >/dev/null || {
      echo "No Stage 3 checkpoint found"
      exit 2
    }
    SCRIPT="$HERE/stage3_quench/in.3_quench_resume"
    EXTRA=(-var restart "$OUTDIR/stage3_quench.checkpoint.*")
    MARKER="stage3_quench.done"
    ;;

  stage4)
    INPUT="$OUTDIR/stage4_checkpoint.restart"
    [[ -f "$INPUT" ]] || {
      echo "Missing $INPUT"
      exit 2
    }
    SCRIPT="$HERE/stage4_relax/in.4_final"
    EXTRA=(-var input_restart "$INPUT")
    MARKER="stage4_final.done"
    ;;

  *)
    echo "Usage: bash resume.sh stage1|stage2|stage3|stage4"
    exit 2
    ;;
esac

echo "[resume] $STAGE with $MPI_NP MPI tasks"

"$MPIEXEC" -np "$MPI_NP" "$LMP" \
  -in "$SCRIPT" \
  "${EXTRA[@]}" \
  -var outdir "$OUTDIR" \
  -var kim_model "$KIM" \
  -log "$LOGDIR/${STAGE}_resume.lammps"

touch "$LOGDIR/$MARKER"
echo "[done] $STAGE"
