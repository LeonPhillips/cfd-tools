#!/bin/bash
# v4 Phase 2 continuation (user, 2026-10-27): restart from the end of s2_cap1, same settings as s2 (cap 1.0, floor 0.1,
# CFL start 1.0), 200000 iterations. Log only: the watcher records check_run.py flags and redraws convergence.png
# every 5 minutes but stops nothing.  Run: tmux new -d -s srp_s5 "bash <this folder>/run_s5.sh"
V4=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs/SA_p0j8.015_p0inf3.247_meshCoarse; R=$V4   # tools moved here from ../SA_p0j8.015_p0inf3.247/ on 2026-10-28
s=s5_cap1_long
exec > $V4/run_s5.log 2>&1
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh
cd $V4/2_jet_staged/$s || exit 1
echo "$s start $(date)"
{ echo "host $(hostname)  ranks $(nproc)"; date; } > $s.log
(sleep 30; STOP_FLAGS=NONE bash $R/watch_stop.sh 2_jet_staged/$s) &
mpirun -n $(nproc) SU2_CFD -t 1 $s.cfg >> $s.log 2>&1
rc=$?
echo "SU2_CFD exit=$rc  $(date)" >> $s.log
wait
echo "$s end exit=$rc $(date): $(tail -1 watch.log 2>/dev/null)"
