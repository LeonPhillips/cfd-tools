#!/bin/bash
# Mesh v4 Phase 2 with a staged CFL cap (user, 2026-10-27): s1 cap 0.5 (5000 it), s2 cap 1.0 (5000), s3 cap 2.0 (10000),
# s4 cap 5.0 (90000), floor 0.1; each stage restarts from the previous one. Each stage runs its full length: the watcher
# only logs check_run.py flags (STOP_FLAGS=NONE) and redraws convergence.png every 5 minutes. The chain stops if a
# stage exits with an error or leaves no restart.dat.
# Run on the desktop node: tmux new -d -s srp_staged "bash <this folder>/run_staged.sh"   (log: run_staged.log)
V4=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs/SA_p0j8.015_p0inf3.247_meshCoarse; R=$V4   # tools moved here from ../SA_p0j8.015_p0inf3.247/ on 2026-10-28
exec > $V4/run_staged.log 2>&1
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh
[ -f $V4/1_freestream/restart.dat ] || { echo "no v4 Phase 1 restart.dat"; exit 1; }
for s in s1_cap0.5 s2_cap1; do   # s3_cap2 and s4_cap5 were deleted on 2026-10-28
  cd $V4/2_jet_staged/$s || exit 1
  echo "$s start $(date)"
  { echo "host $(hostname)  ranks $(nproc)"; date; } > $s.log
  (sleep 30; STOP_FLAGS=NONE bash $R/watch_stop.sh 2_jet_staged/$s) &
  mpirun -n $(nproc) SU2_CFD -t 1 $s.cfg >> $s.log 2>&1
  rc=$?
  echo "SU2_CFD exit=$rc  $(date)" >> $s.log
  wait
  echo "$s end exit=$rc $(date): $(tail -1 watch.log 2>/dev/null)"
  [ $rc -eq 0 ] && [ -f restart.dat ] || { echo "stage $s failed: chain stopped"; exit 1; }
done
echo "all stages done $(date)"
