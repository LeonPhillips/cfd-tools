#!/bin/bash
# Mesh v4: Phase 1, then check_phase1.py, then (only if it prints RESULT PASS) Phase 2 with the watcher in tmux
# session srp_watch. Run on the desktop node: tmux new -d -s srp_v4 "bash <this folder>/run_v4.sh"
# Log of this script: run_v4.log
V4=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs/SA_p0j8.015_p0inf3.247_meshCoarse; R=$V4   # tools moved here from ../SA_p0j8.015_p0inf3.247/ on 2026-10-28
exec > $V4/run_v4.log 2>&1
echo "Phase 1 start $(date)"
bash $R/run_su2.sh v4_phase1
tail -1 $V4/1_freestream/phase1.log
(unset LD_PRELOAD; module load paraview/5.13 >/dev/null 2>&1
 /apps/paraview/5.13/bin/pvbatch --mesa $V4/check_phase1.py $V4/1_freestream 2>&1 | grep -v Warning) > $V4/1_freestream/phase1_check.txt
cat $V4/1_freestream/phase1_check.txt
if ! grep -q "RESULT PASS" $V4/1_freestream/phase1_check.txt; then
  echo "Phase 1 check failed: Phase 2 not started $(date)"; exit 1
fi
echo "Phase 1 checked $(date). The original Phase 2 here (2_jet, CFL 1 -> 50) was deleted on 2026-10-28; Phase 2 is now 2_jet_staged/."; exit 0
echo "Phase 2 start $(date)"
tmux new -d -s srp_watch "sleep 60; bash $R/watch_stop.sh 2_jet"
bash $R/run_su2.sh v4_phase2
tail -1 $V4/2_jet/phase2.log
echo "done $(date)"
