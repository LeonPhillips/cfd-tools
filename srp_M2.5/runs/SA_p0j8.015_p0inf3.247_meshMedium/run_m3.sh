#!/bin/bash
# meshMedium stage m3 (user, 2026-10-28): restart from m2_cap1 (stopped at iteration 26673), CFL cap 3.0, start 1.0,
# floor 0.1, 20000 iterations. Log-only watcher (STOP_FLAGS=NONE): it records check_run.py flags and redraws convergence.png.
# Run on the desktop node: tmux new -d -s srp_m3 "bash <this folder>/run_m3.sh"   (log: run_m3.log)
C=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs/SA_p0j8.015_p0inf3.247_meshCoarse   # tools live here
M=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs/SA_p0j8.015_p0inf3.247_meshMedium
s=m3_cap3
exec > $M/run_m3.log 2>&1
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh
cd $M/2_jet_staged/$s || exit 1
echo "$s start $(date)"
{ echo "host $(hostname)  ranks $(nproc)"; date; } > $s.log
(sleep 30; STOP_FLAGS=NONE bash $C/watch_stop.sh ../SA_p0j8.015_p0inf3.247_meshMedium/2_jet_staged/$s) &
mpirun -n $(nproc) SU2_CFD -t 1 $s.cfg >> $s.log 2>&1
rc=$?
echo "SU2_CFD exit=$rc  $(date)" >> $s.log
wait
echo "$s end exit=$rc $(date): $(tail -1 watch.log 2>/dev/null)"
