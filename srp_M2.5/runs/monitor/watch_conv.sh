#!/bin/bash
# Live convergence monitor for one stage (user, 2026-11-04; replaces watch_fj.sh for new runs).
# Usage: [PREV=<previous stage folder>] [HOLD=1] bash watch_conv.sh <stage folder>      (ONCE=1: one pass, then exit)
# Needs in the stage config: OUTPUT_WRT_FREQ= <restart>, 1000, <surface> and WRT_VOLUME_OVERWRITE= NO (numbered
# flow_<iter>.vtu every 1000 iterations), and CUSTOM_OUTPUTS from probes_md.py with CUSTOM in HISTORY_OUTPUT.
# Every minute, for each new flow_<iter>.vtu (older than 60 s, so SU2 has finished writing it):
#   snap.py      -> snaps/flow_<iter>.npz (0.025 mm box x 0-20, r 0-9.5 mm + Mach on the six lines); the vtu is deleted
#   conv_update.py -> box_change.csv, line_change.csv, features.csv; keeps only the latest two snapshots
#   three figures, each copied to the run folder as *_latest.png:
#     convergence.png    (plot_conv_roi.py: residuals, probes, box change, mass flows, CFL, difference map x 0-16, r 0-6)
#     conv_lines.png     (plot_conv_lines.py: Mach change on the lines r 0.01/1.3/2.3/3/5 mm and the shear layer)
#     conv_features.png  (plot_conv_features.py: Mach disk x and r, shock crossing, bow shock, stagnation)
# Every 2 minutes: check_run.py into watch.log. Nothing here stops the run (convergence.md), except HOLD=1: end the stage
# (SIGTERM, clean write) once the average CFL is unchanged for HOLD_N (5000) iterations, logged CFLHOLD.
# Only SU2_CFD processes whose working directory is <stage folder> are watched. Log: <stage folder>/watch.log
MON=$(dirname $(readlink -f $0))
C=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs/SA_p0j8.015_p0inf3.247_meshCoarse   # check_run.py
D=$(readlink -f ${1:?stage folder}); HOLD=${HOLD:-0}; HOLD_N=${HOLD_N:-5000}; PREV=${PREV:--}
RUN=$(dirname $(dirname $D))
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh >/dev/null
mkdir -p $D/snaps
mine() { for p in $(pgrep -u $USER -x SU2_CFD); do [ "$(readlink /proc/$p/cwd)" = "$D" ] && echo $p; done; }
held() {  # number of trailing history rows with the same average CFL as the last row
  awk -F, 'NR == 1 { for (i = 1; i <= NF; i++) { h = $i; gsub(/[" ]/, "", h); if (h == "AvgCFL") c = i }; nf = NF; next }
           NF == nf { v[++n] = $c + 0 }
           END { k = 0; for (i = n; i >= 1 && (v[i] - v[n] <= 1e-6 * v[n]) && (v[n] - v[i] <= 1e-6 * v[n]); i--) k++; print k }' \
    <(head -1 $D/history.csv; tail -n $((HOLD_N + 1)) $D/history.csv)
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
    ( unset LD_PRELOAD; module load paraview/5.13 >/dev/null 2>&1
      pvbatch --mesa $MON/snap.py $D/snaps "${new[@]}" 2>&1 | grep -vi warning >> $D/watch.log )
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
sleep 70; snaps                             # final snapshot and plots
echo "$(date +%T) watcher done" >> $D/watch.log
