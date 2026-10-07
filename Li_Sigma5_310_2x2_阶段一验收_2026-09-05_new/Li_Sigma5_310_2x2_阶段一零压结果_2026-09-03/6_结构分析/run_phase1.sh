#!/bin/bash
# Phase-1 runner: true zero-pressure NPT full pipeline (guide 2026-09-02)
# Idempotent: each stage guarded by logs/p1_*.done marker.
# Resume after interruption: this script restarts any unfinished stage from
# its own starting point; if a stage crashed mid-run, delete its .done (none
# written) and use the latest phase1/*.checkpoint.* restart manually.
cd ~/Li_meta_2x2
NP=12
LMP=/usr/local/bin/lmp-kim
KIM=SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
IN=phase1_inputs

# ---- Stage A: 0 K full aniso relaxation (ordered GB, zero stress) ----
if [ ! -f logs/p1_min0.done ]; then
  echo "=== P1 stage A: 0K aniso relax ==="
  mpirun -np $NP $LMP -in $IN/in.p1_min0 -var kim_model $KIM -log logs/p1_min0.lammps > logs/p1_min0.out 2>&1
  grep -q P1_MIN0_DONE logs/p1_min0.lammps && touch logs/p1_min0.done && echo "min0 OK"
fi

# ---- Stage B: 363 K iso 0 bar NPT anneal, 25 ps ----
if [ ! -f logs/p1_anneal.done ]; then
  echo "=== P1 stage B: 363K iso zero-P anneal ==="
  mpirun -np $NP $LMP -in $IN/in.p1_anneal -var kim_model $KIM -var seed 12345 -log logs/p1_anneal.lammps > logs/p1_anneal.out 2>&1
  grep -q P1_ANNEAL_DONE logs/p1_anneal.lammps && touch logs/p1_anneal.done && echo "anneal OK"
fi

# ---- Stage C: 0 K Case-B diagnostic (fixed box) ----
if [ ! -f logs/p1_casecheck.done ]; then
  echo "=== P1 stage C: 0K Case-B check ==="
  mpirun -np $NP $LMP -in $IN/in.p1_casecheck -var kim_model $KIM -log logs/p1_casecheck.lammps > logs/p1_casecheck.out 2>&1
  grep -q P1_CASECHECK_DONE logs/p1_casecheck.lammps && touch logs/p1_casecheck.done && echo "casecheck OK"
fi

# ---- Stage D: quench 10 K/ps, iso zero-P NPT ----
if [ ! -f logs/p1_quench.done ]; then
  echo "=== P1 stage D: quench 363->10K ==="
  mpirun -np $NP $LMP -in $IN/in.p1_quench -var kim_model $KIM -log logs/p1_quench.lammps > logs/p1_quench.out 2>&1
  grep -q P1_QUENCH_DONE logs/p1_quench.lammps && touch logs/p1_quench.done && echo "quench OK"
fi

# ---- Stage E: 10 K hold 5 ps, iso zero-P NPT ----
if [ ! -f logs/p1_hold.done ]; then
  echo "=== P1 stage E: 10K hold ==="
  mpirun -np $NP $LMP -in $IN/in.p1_hold -var kim_model $KIM -log logs/p1_hold.lammps > logs/p1_hold.out 2>&1
  grep -q P1_HOLD_DONE logs/p1_hold.lammps && touch logs/p1_hold.done && echo "hold OK"
fi

echo "=== PHASE1_RUNNER_COMPLETE ==="
ls logs/p1_*.done 2>/dev/null
