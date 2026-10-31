#!/usr/bin/env bash
# v2 Sigma5 (310)[001] BCC Li -- resume an interrupted workflow.
#
# run_all.sh already skips finished stages (logs/<stage>.done markers), so
# "resume" is simply re-running it; it continues from the first unfinished
# stage using the restart files written by the previous stage.
#
# For mid-stage continuation (e.g. a quench that was killed after a
# checkpoint), use in.resume with the checkpoint restart, for example:
#   lmp -in in.resume -var restart ../outputs/stage3_quench.checkpoint \
#       -var stage quench -var nsteps 30000 \
#       -var tstart 300.0 -var tstop 1.0 -var outdir ../outputs
set -euo pipefail
echo "[resume] re-running staged workflow (finished stages are skipped)..."
exec "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_all.sh" "$@"
