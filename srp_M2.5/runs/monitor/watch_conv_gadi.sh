#!/bin/bash
# Live convergence monitor for one stage, Gadi version (user, 2026-11-05). Same logic and outputs as watch_conv.sh (M3):
# every minute, each new flow_<iter>.vtu (older than 60 s) -> snap.py -> snaps/flow_<iter>.npz (the vtu is deleted),
# conv_update.py -> box_change.csv / line_change.csv / features.csv, then the three figures (convergence.png,
# conv_lines.png, conv_features.png), each copied to the run folder as *_latest.png. Every 2 minutes check_run.py into
# watch.log. Nothing here stops the run (convergence.md).
# Differences from watch_conv.sh:
#   - ParaView 6.1.1 (module paraview/6.1.1-mesa, an OSMesa build: pvbatch has no --mesa flag). It loads openmpi/5.0.5,
#     so it runs in a subshell with the SU2 environment (Intel MPI, venv, PYTHONPATH) removed.
#   - Gadi SU2 environment; check_run.py is the copy in this folder.
#   - No HOLD option.
# Usage: [PREV=<previous stage folder>] bash watch_conv_gadi.sh <stage folder>      (ONCE=1: one pass, then exit)
# Only SU2_CFD processes whose working directory is <stage folder> are watched. Log: <stage folder>/watch.log
MON=$(dirname $(readlink -f $0))
PV_MODULE=${PV_MODULE:-paraview/6.1.1-mesa}
D=$(readlink -f ${1:?stage folder}); PREV=${PREV:--}
RUN=$(dirname $(dirname $D))
source /scratch/os85/lp3540/su2/SU2-v8.5.0/su2_env.sh >/dev/null
type module >/dev/null 2>&1 || source /opt/Modules/v4.3.0/init/bash
mkdir -p $D/snaps
mine() { for p in $(pgrep -u $USER -x SU2_CFD); do [ "$(readlink /proc/$p/cwd)" = "$D" ] && echo $p; done; }
pv() {  # pv <script> <args...>: pvbatch from ParaView 6.1 in a clean subshell
  ( unset LD_PRELOAD PYTHONPATH PYTHONHOME VIRTUAL_ENV; module purge >/dev/null 2>&1; module load $PV_MODULE >/dev/null 2>&1
    pvbatch "$@" )
}
plots() {
  python $MON/plot_conv_roi.py $D $PREV > /dev/null 2>> $D/watch.log && cp $D/convergence.png $RUN/convergence_latest.png
  python $MON/plot_conv_lines.py $D $PREV > /dev/null 2>> $D/watch.log && cp $D/conv_lines.png $RUN/conv_lines_latest.png
  python $MON/plot_conv_features.py $D > /dev/null 2>> $D/watch.log && cp $D/conv_features.png $RUN/conv_features_latest.png
}
snaps() {
  local new=() f
  for f in $D/flow_*.vtu; do
    [ -e "$f" ] && [ $(( $(date +%s) - $(stat -L -c %Y $f) )) -gt 60 ] && new+=($f)
  done
  if [ ${#new[@]} -gt 0 ]; then
    pv $MON/snap.py $D/snaps "${new[@]}" 2>&1 | grep -vi warning >> $D/watch.log
    for f in "${new[@]}"; do [ -f $D/snaps/$(basename ${f%.vtu}).npz ] && rm -f $f; done
    python $MON/conv_update.py $D $PREV 2>> $D/watch.log
  fi
  [ -f $D/history.csv ] && plots
}
if [ "${ONCE:-0}" = 1 ]; then snaps; echo "$(date +%T) watcher ONCE pass done" >> $D/watch.log; exit 0; fi
t=0
while [ -n "$(mine)" ]; do
  sleep 60; t=$((t + 1))
  [ -f $D/history.csv ] || continue
  [ $((t % 2)) -eq 0 ] && echo "$(date +%T) $(python $MON/check_run.py $D 2>&1)" >> $D/watch.log
  snaps
done
while [ -n "$(mine)" ]; do sleep 10; done   # let SU2 finish its final write
sleep 70; snaps                             # final snapshot and plots
echo "$(date +%T) watcher done" >> $D/watch.log
