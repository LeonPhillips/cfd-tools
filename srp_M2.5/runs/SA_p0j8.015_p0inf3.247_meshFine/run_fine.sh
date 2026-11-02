#!/bin/bash
# meshFine Phase 2 (user-approved plan, 2026-10-28): f1 = CFL cap 0.5 for 5000 iterations from the interpolated medium m3
# field (0_interp), f2 = cap 1.0 for 25000 from f1, f3 = cap 3.0 for 20000 from f2. Floor 0.1, fully second order. Each
# stage runs its full length: the watcher only logs check_run.py flags and redraws convergence.png (STOP_FLAGS=NONE).
# Run on the desktop node: tmux new -d -s srp_fine "bash <this folder>/run_fine.sh"   (log: run_fine.log)
C=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs/SA_p0j8.015_p0inf3.247_meshCoarse   # tools live here
F=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs/SA_p0j8.015_p0inf3.247_meshFine
exec > $F/run_fine.log 2>&1
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh
[ -f $F/0_interp/restart.dat ] || { echo "no interpolated restart.dat"; exit 1; }
for s in f1_cap0.5 f2_cap1 f3_cap3; do
  cd $F/2_jet_staged/$s || exit 1
  echo "$s start $(date)"
  { echo "host $(hostname)  ranks $(nproc)"; date; } > $s.log
  (sleep 30; STOP_FLAGS=NONE bash $C/watch_stop.sh ../SA_p0j8.015_p0inf3.247_meshFine/2_jet_staged/$s) &
  mpirun -n $(nproc) SU2_CFD -t 1 $s.cfg >> $s.log 2>&1
  rc=$?
  echo "SU2_CFD exit=$rc  $(date)" >> $s.log
  wait
  echo "$s end exit=$rc $(date): $(tail -1 watch.log 2>/dev/null)"
  [ $rc -eq 0 ] && [ -f restart.dat ] || { echo "stage $s failed: chain stopped"; exit 1; }
done
echo "all stages done $(date)"
