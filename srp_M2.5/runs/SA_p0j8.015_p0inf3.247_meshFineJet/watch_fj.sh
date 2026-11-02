#!/bin/bash
# Watcher for one meshFineJet stage (2026-10-29). Usage: [HOLD=1] [PREV=<previous stage folder>] bash watch_fj.sh <stage folder>
# Every minute:
#   - if HOLD=1, end the stage (SIGTERM: SU2 finishes the iteration and writes restart/flow/surface) once the average CFL
#     has not changed for HOLD_N (5000) iterations (user rule, 2026-10-29). Logs CFLHOLD.
#   - snapshots (user, 2026-10-29): each new flow_<iter>.vtu (one every 1000 iterations) is resampled over the jet box
#     (box_snap.py) and then deleted; box_change.py logs the change and keeps only the latest two box snapshots. Then
#     convergence.png is redrawn (plot_conv_roi.py) and copied to <run folder>/convergence_latest.png, which the PC
#     copies into the Simulation figures folder. SU2's own flow.vtu (latest full field) is kept.
# Every 2 minutes: log check_run.py (information only; nothing else ever stops the run).
# Only SU2_CFD processes whose working directory is <stage folder> are watched. Log: <stage folder>/watch.log
FJ=$(dirname $(readlink -f $0))
C=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs/SA_p0j8.015_p0inf3.247_meshCoarse   # check_run.py
D=$(readlink -f ${1:?stage folder}); HOLD=${HOLD:-0}; HOLD_N=${HOLD_N:-5000}; PREV=${PREV:--}
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh
mkdir -p $D/snaps
mine() { for p in $(pgrep -u $USER -x SU2_CFD); do [ "$(readlink /proc/$p/cwd)" = "$D" ] && echo $p; done; }
# number of trailing history rows with the same average CFL as the last row
held() {
  awk -F, 'NR == 1 { for (i = 1; i <= NF; i++) { h = $i; gsub(/[" ]/, "", h); if (h == "AvgCFL") c = i }; nf = NF; next }
           NF == nf { v[++n] = $c + 0 }
           END { k = 0; for (i = n; i >= 1 && (v[i] - v[n] <= 1e-6 * v[n]) && (v[n] - v[i] <= 1e-6 * v[n]); i--) k++; print k }' \
    <(head -1 $D/history.csv; tail -n $((HOLD_N + 1)) $D/history.csv)
}
# resample, log and plot new snapshots (files older than 60 s, so SU2 has finished writing them)
snaps() {
  local new=() f
  for f in $D/flow_*.vtu; do
    [ -f "$f" ] && [ $(( $(date +%s) - $(stat -c %Y $f) )) -gt 60 ] && new+=($f)
  done
  [ ${#new[@]} -gt 0 ] || return 0
  ( unset LD_PRELOAD; module load paraview/5.13 >/dev/null 2>&1
    pvbatch --mesa $FJ/box_snap.py $D/snaps "${new[@]}" 2>&1 | grep -vi warning >> $D/watch.log )
  for f in "${new[@]}"; do [ -f $D/snaps/$(basename ${f%.vtu}).npz ] && rm -f $f; done
  python $FJ/box_change.py $D $PREV 2>> $D/watch.log
  python $FJ/plot_conv_roi.py $D $PREV > /dev/null 2>> $D/watch.log && cp $D/convergence.png $FJ/convergence_latest.png
}
t=0
while [ -n "$(mine)" ]; do
  sleep 60; t=$((t + 1))
  [ -f $D/history.csv ] || continue
  if [ "$HOLD" = 1 ]; then
    k=$(held)
    if [ "${k:-0}" -ge "$HOLD_N" ]; then
      echo "$(date +%T) CFLHOLD: average CFL unchanged for $k iterations: ending the stage (SIGTERM)" >> $D/watch.log
      kill -TERM $(mine)
      break
    fi
  fi
  [ $((t % 2)) -eq 0 ] && echo "$(date +%T) $(python $C/check_run.py $D 2>&1)" >> $D/watch.log
  snaps
done
while [ -n "$(mine)" ]; do sleep 10; done   # let SU2 finish its final write
sleep 70; snaps                             # final snapshot and plot
echo "$(date +%T) watcher done" >> $D/watch.log
