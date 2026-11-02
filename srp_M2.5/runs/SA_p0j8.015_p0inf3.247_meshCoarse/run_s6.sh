#!/bin/bash
# v4 Phase 2 continuation (user, 2026-10-28): restart from the end of s5_cap1_long, CFL cap 3.0, floor 0.1, start 1.0,
# up to 200000 iterations. The watcher stops SU2 on CONVERGED (CD band < 1 % over 10000 it and imbalance < 0.5 %) or
# DIVERGE; CFLCYCLE is only logged.  Run: tmux new -d -s srp_s6 "bash <this folder>/run_s6.sh"
V4=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs/SA_p0j8.015_p0inf3.247_meshCoarse; R=$V4   # tools moved here from ../SA_p0j8.015_p0inf3.247/ on 2026-10-28
s=s6_cap3
exec > $V4/run_s6.log 2>&1
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh
cd $V4/2_jet_staged/$s || exit 1
echo "$s start $(date)"
{ echo "host $(hostname)  ranks $(nproc)"; date; } > $s.log
(sleep 30; STOP_FLAGS="CONVERGED DIVERGE" bash $R/watch_stop.sh 2_jet_staged/$s) &
mpirun -n $(nproc) SU2_CFD -t 1 $s.cfg >> $s.log 2>&1
rc=$?
echo "SU2_CFD exit=$rc  $(date)" >> $s.log
wait
echo "$s end exit=$rc $(date): $(grep -v 'watcher done' watch.log 2>/dev/null | tail -1)"
