#!/bin/bash
# Live convergence monitor for one free-jet stage on Gadi (user, 2026-11-23). As monitor/watch_conv_gadi.sh, with the
# free-jet lines and features: every minute each new flow_<iter>.vtu (older than 60 s) -> snap_fj.py (ParaView 6.1)
# -> snaps/flow_<iter>.npz (the vtu is deleted), conv_fj.py -> line_history.npz, line_change.csv, features.csv, then
# plot_fj.py redraws the chain figures of the whole run (copied as *_latest.png). One status line per pass in
# watch.log. Nothing here stops the run (convergence.md).
# Usage: bash watch_fj_gadi.sh <stage folder> "<title>"      (ONCE=1: one pass, then exit)
T=$(dirname $(readlink -f $0))
PV_MODULE=${PV_MODULE:-paraview/6.1.1-mesa}
D=$(readlink -f ${1:?stage folder}); TITLE=${2:-}
RUN=$(dirname $(dirname $D))
source /scratch/os85/lp3540/su2/SU2-v8.5.0/su2_env.sh >/dev/null
type module >/dev/null 2>&1 || source /opt/Modules/v4.3.0/init/bash
mkdir -p $D/snaps
mine() { for p in $(pgrep -u $USER -x SU2_CFD); do [ "$(readlink /proc/$p/cwd)" = "$D" ] && echo $p; done; }
pv() { ( unset LD_PRELOAD PYTHONPATH PYTHONHOME VIRTUAL_ENV; module purge >/dev/null 2>&1
         module load $PV_MODULE >/dev/null 2>&1; pvbatch "$@" ); }
pass() {
  local new=() f
  for f in $D/flow_*.vtu; do
    [ -e "$f" ] && [ $(( $(date +%s) - $(stat -L -c %Y $f) )) -gt 60 ] && new+=($f)
  done
  if [ ${#new[@]} -gt 0 ]; then
    pv $T/snap_fj.py $D/snaps "${new[@]}" 2>&1 | grep -vi warning >> $D/watch.log
    for f in "${new[@]}"; do
      n=$(basename $f .vtu); printf -v npz 'flow_%06d.npz' $((10#${n#flow_}))
      [ -f $D/snaps/$npz ] && rm -f $f
    done
    python $T/conv_fj.py $D 2>> $D/watch.log
  fi
  if [ -f $D/history.csv ]; then
    echo "$(date +%T) $(python $T/plot_fj.py $RUN "$TITLE" 2>&1 | tail -1)" >> $D/watch.log
    for p in convergence conv_lines conv_features; do [ -f $RUN/${p}_chain.png ] && cp $RUN/${p}_chain.png $RUN/${p}_latest.png; done
  fi
}
if [ "${ONCE:-0}" = 1 ]; then pass; echo "$(date +%T) watcher ONCE pass done" >> $D/watch.log; exit 0; fi
while [ -n "$(mine)" ]; do sleep 60; [ -f $D/history.csv ] && pass; done
sleep 70; pass   # the final snapshot and figures
echo "$(date +%T) watcher done" >> $D/watch.log
