#!/bin/bash
# meshMedium Phase 2 (user-approved plan, 2026-10-28): m1 = CFL cap 0.5 for 5000 iterations from the interpolated s6
# field (0_interp), then m2 = cap 1.0 for 45000 iterations from m1. Floor 0.1, fully second order. Each stage runs its
# full length: the watcher only logs check_run.py flags and redraws convergence.png (STOP_FLAGS=NONE).
# Run on the desktop node: tmux new -d -s srp_med "bash <this folder>/run_medium.sh"   (log: run_medium.log)
C=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs/SA_p0j8.015_p0inf3.247_meshCoarse   # tools live here
M=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs/SA_p0j8.015_p0inf3.247_meshMedium
exec > $M/run_medium.log 2>&1
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh
[ -f $M/0_interp/restart.dat ] || { echo "no interpolated restart.dat"; exit 1; }
for s in m1_cap0.5 m2_cap1; do
  cd $M/2_jet_staged/$s || exit 1
  echo "$s start $(date)"
  { echo "host $(hostname)  ranks $(nproc)"; date; } > $s.log
  (sleep 30; STOP_FLAGS=NONE bash $C/watch_stop.sh ../SA_p0j8.015_p0inf3.247_meshMedium/2_jet_staged/$s) &
  mpirun -n $(nproc) SU2_CFD -t 1 $s.cfg >> $s.log 2>&1
  rc=$?
  echo "SU2_CFD exit=$rc  $(date)" >> $s.log
  wait
  echo "$s end exit=$rc $(date): $(tail -1 watch.log 2>/dev/null)"
  [ $rc -eq 0 ] && [ -f restart.dat ] || { echo "stage $s failed: chain stopped"; exit 1; }
done
echo "all stages done $(date)"
