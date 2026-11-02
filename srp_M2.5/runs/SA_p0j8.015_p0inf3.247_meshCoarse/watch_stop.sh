#!/bin/bash
# Stop a running phase cleanly when check_run.py flags CONVERGED, OSCILLATE, CORRUPT or DIVERGE (checked
# every 2 minutes), and redraw <phase folder>/convergence.png every 5 minutes (overwritten, so only
# the latest plot is kept).
# Usage, on the desktop node: tmux new -d -s srp_watch "bash watch_stop.sh <phase folder>"
# SU2 treats SIGTERM as "converged": it finishes the iteration and writes restart/flow/surface files.
# Only SU2_CFD processes whose working directory is <phase folder> are watched and stopped, so other SU2 runs
# on the node are left alone (fixed 2026-10-27; it used to pkill every SU2_CFD).
# Log: <phase folder>/watch.log
# STOP_FLAGS (env) = the check_run.py flags that stop the run; default CONVERGED OSCILLATE CORRUPT DIVERGE.
# Other flags are only logged. E.g. STOP_FLAGS=CFLCYCLE bash watch_stop.sh <folder>
R=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs/SA_p0j8.015_p0inf3.247_meshCoarse
D=$R/${1:?phase folder}
STOP_FLAGS=${STOP_FLAGS:-CONVERGED OSCILLATE CORRUPT DIVERGE}
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh
set -u
cd $R
DR=$(readlink -f $D)
mine() { for p in $(pgrep -u $USER -x SU2_CFD); do [ "$(readlink /proc/$p/cwd)" = "$DR" ] && echo $p; done; }
t=0
while [ -n "$(mine)" ]; do
  sleep 60; t=$((t + 1))
  [ -f $D/history.csv ] || continue
  [ $((t % 5)) -eq 0 ] && python plot_history.py $D 10000 > /dev/null 2>&1
  [ $((t % 2)) -eq 0 ] || continue
  s=$(python check_run.py $D 2>&1)
  echo "$(date +%T) $s" >> $D/watch.log
  for f in $STOP_FLAGS; do
    case $s in *"flags "*$f*)
      echo "$(date +%T) stopping SU2 (SIGTERM) on $f" >> $D/watch.log
      kill -TERM $(mine)
      break 2 ;;
    esac
  done
done
python plot_history.py $D 10000 > /dev/null 2>&1   # final plot
echo "$(date +%T) watcher done" >> $D/watch.log
