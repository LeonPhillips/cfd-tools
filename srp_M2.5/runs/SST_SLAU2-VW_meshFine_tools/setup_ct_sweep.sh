#!/bin/bash
# CT sweep (user, 2026-11-26): Schlieren1 conditions from ct_sweep_matrix.txt (copied from
# cases/srp_M2.5/CT_sweep_conditions.md on Market), run with the 6.6 bar SLAU2-VW meshFine setup.
# The file-name CT/CTa/CTd are labels only; the RANS results are not matched to them (user, 2026-11-25).
# Boundary conditions (confirmed, user 2026-11-25/26):
#   p0inf = tunnel p0 gauge + 1.01325 bar (absolute); T0inf = 288.15 K, M 2.5: T = 128.067 K, U = 567.16 m/s,
#     p = p0inf / 2.25^3.5.
#   p0j = 0.887732 * nozzle target - 0.044750 bar abs (purple fit, "p0 (direct) fit, slope 0.888" in
#     Leon_pitotTube_data/pitot_pdrop_p0_vs_target.fig); T0j = 288.15 K.
# Each case restarts from a converged solution on the same mesh given on the command line (the nearest case).
# Builds, under /scratch/os85/lp3540/su2/cases/srp_M2.5/runs/CT-sweep/:
#   CT<label>_p0j<P>_p0inf<Q>/
#     0_start/                  copies of the start stage's restart.dat and flow.vtu, SOURCE.txt
#     2_jet_staged/t1_cap0.5 ... t1-t3 = 1k/1k/8k (CFL 0.1->0.5, 0.5->1, 1->3); auto-stop chunks are added by the driver
#     title.txt
# Each config is the 6.6 bar t5 config with only: the freestream pressure (FREESTREAM_PRESSURE and
# MARKER_SUPERSONIC_INLET), the jet p0 (MARKER_INLET), ITER, CFL start and cap, SOLUTION_FILENAME, comments.
# Usage (login node; copies ~160 MB per case):
#   bash setup_ct_sweep.sh CT0.1626=<start stage folder> CT0.5202=<start stage folder> ...
#   CT0.1931=@CT0.1808: start from the final solution of sweep case CT0.1808 (copied when the job starts; submit
#   with -W depend=afterany:<CT0.1808 job>; the parent must have AUTOSTOP_PASS.txt or AUTOSTOP_CAP.txt by then).
#   then per case:  qsub -N <name> -l walltime=04:30:00 -v RUN=<case folder>,AUTOSTOP=1 chain_sr.pbs
#   (OUT=<folder> builds elsewhere, e.g. for a test)
set -e
R=/scratch/os85/lp3540/su2/cases/srp_M2.5/runs
OUT=${OUT:-$R/CT-sweep}
TEMPLATE=$R/SST_p0j6.6_p0inf3.247_meshFine_SLAU2-VW/2_jet_staged/t5_cap3/t5_cap3.cfg
TOOLS=$(dirname $(readlink -f $0))
MATRIX=$TOOLS/ct_sweep_matrix.txt
COMMIT=$(git -C $TOOLS rev-parse --short HEAD 2>/dev/null || echo unknown)
# Auto-stop sweep (user, 2026-11-27): only t1-t3 (10k) are built here; run_chain_gadi.sh with AUTOSTOP=1 adds 10k chunks
# (t4_c01, ...) until autostop_check.py passes (from 50k) or 250k is reached. STAGES="name:ITER:CFL0:cap ..." overrides.
STAGES=${STAGES:-"t1_cap0.5:1000:0.1:0.5 t2_cap1:1000:0.5:1.0 t3_cap3:8000:1.0:3.0"}

mkcfg() {  # mkcfg <out> <label> <p0j Pa> <p0inf Pa> <stage> <ITER> <CFL0> <cap> <solution> <start>
  local out=$1 lab=$2 pj=$3 p0=$4
  local T=128.067 U=567.16
  local p=$(awk -v a=$p0 'BEGIN{printf "%.1f", a/2.25^3.5}')
  local rho=$(awk -v p=$p 'BEGIN{printf "%.5f", p/(287.058*128.067)}')
  local pjb=$(awk -v a=$pj 'BEGIN{printf "%.4f", a/1e5}') p0b=$(awk -v a=$p0 'BEGIN{printf "%.5f", a/1e5}')
  {
    echo "% srp_M2.5 CT sweep ($lab, Schlieren1), meshFine, k-omega SST, SLAU2 + Venkatakrishnan-Wang (user, 2026-11-26)."
    echo "% p0j $pjb bar (purple fit), p0inf $p0b bar abs. Stage $5: $6 iterations, CFL $7 -> cap $8 (x1.01), floor 0.1."
    echo "% Made by SST_SLAU2-VW_meshFine_tools/setup_ct_sweep.sh (cfd-tools $COMMIT) from $TEMPLATE"
    echo "% (only the freestream pressure, jet p0, ITER, CFL and solution path changed). Nothing stops on convergence."
    sed -n '/^% ---- physics/,$p' $TEMPLATE
  } | sed -e "s|^% ---- freestream: .*|% ---- freestream: M2.5 from p0inf = $p0b bar abs, T0 = 288.15 K (isentropic; user 2026-11-25) ----|" \
          -e "s|^% T = 288.15/2.25 .*|% T = 288.15/2.25 = $T K; p = $p0/2.25^3.5 = $p Pa; U = $U m/s in -x|" \
          -e "s|^% (AOA = 180 deg); rho = .*|% (AOA = 180 deg); rho = $rho kg/m3. Restarted from ${10}.|" \
          -e "s|^FREESTREAM_PRESSURE=.*|FREESTREAM_PRESSURE= $p|" \
          -e "s|^MARKER_SUPERSONIC_INLET=.*|MARKER_SUPERSONIC_INLET= ( inlet_freestream, $T, $p, -$U, 0.0, 0.0 )|" \
          -e "s|^% Phase 2: inlet_jet .*|% inlet_jet: total conditions, p0j = $pjb bar (purple fit), T0j = 288.15 K, flow in +x (user 2026-11-26).|" \
          -e "s|^% inlet_jet: total conditions.*|% inlet_jet: total conditions, p0j = $pjb bar (purple fit), T0j = 288.15 K, flow in +x (user 2026-11-26).|" \
          -e "s|^MARKER_INLET=.*|MARKER_INLET= ( inlet_jet, 288.15, $pj, 1.0, 0.0, 0.0 )   % p0j $pjb bar (total), T0j 288.15 K (user, 2026-11-26)|" \
          -e "s|^ITER=.*|ITER= $6|" \
          -e "s|^CFL_NUMBER=.*|CFL_NUMBER= $7|" \
          -e "s|^CFL_ADAPT_PARAM=.*|CFL_ADAPT_PARAM= ( 0.5, 1.01, 0.1, $8, 0.001, 0 )|" \
          -e "s|^SOLUTION_FILENAME=.*|SOLUTION_FILENAME= $9|" > $out
  for k in "FREESTREAM_PRESSURE= $p" "MARKER_SUPERSONIC_INLET= ( inlet_freestream, $T, $p, -$U" \
           "MARKER_INLET= ( inlet_jet, 288.15, $pj," "ITER= $6" "CFL_NUMBER= $7" "CONV_RESIDUAL_MINVAL= -20" \
           "CONV_NUM_METHOD_FLOW= SLAU2" "SLOPE_LIMITER_FLOW= VENKATAKRISHNAN_WANG" "RESTART_SOL= YES" \
           "SOLUTION_FILENAME= $9" "WRT_VOLUME_OVERWRITE= NO"; do
    grep -qF "$k" $out || { echo "CHECK FAILED: '$k' missing in $out"; exit 1; }
  done
}

[ $# -ge 1 ] || { echo "usage: bash setup_ct_sweep.sh <label>=<start stage folder> ..."; exit 1; }
[ -s $TEMPLATE ] && [ -s $MATRIX ] || { echo "missing $TEMPLATE or $MATRIX"; exit 1; }
mkdir -p $OUT
for arg in "$@"; do
  lab=${arg%%=*}; src=${arg#*=}; PARENT=
  if [[ $src == @* ]]; then   # deferred start: the final solution of another sweep case, copied when the job starts
    PARENT=$(ls -d $OUT/${src#@}_p0j* 2>/dev/null | head -1)
    [ -n "$PARENT" ] && [ -d $PARENT/2_jet_staged ] || { echo "$lab: parent case ${src#@} not found in $OUT"; exit 1; }
    START="final stage of $PARENT (at job start)"
  else
    START=$(readlink -f $src)
    [ -s $START/restart.dat ] && [ -s $START/flow.vtu ] || { echo "$lab: no restart.dat/flow.vtu in $START"; exit 1; }
  fi
  row=$(awk -v l=$lab '$1==l' $MATRIX)
  [ -n "$row" ] || { echo "$lab not in $MATRIX"; exit 1; }
  read _ g t <<< "$row"
  p0=$(awk -v g=$g 'BEGIN{printf "%.1f", (g+1.01325)*1e5}')
  pj=$(awk -v t=$t 'BEGIN{printf "%.1f", (0.887732*t-0.044750)*1e5}')
  P=$(awk -v a=$pj 'BEGIN{printf "%.3f", a/1e5}'); Q=$(awk -v a=$p0 'BEGIN{printf "%.3f", a/1e5}')
  RUN=$OUT/${lab}_p0j${P}_p0inf${Q}
  [ -e $RUN ] && { echo "$RUN exists: not touched"; exit 1; }
  mkdir -p $RUN/0_start
  if [ -n "$PARENT" ]; then
    echo $PARENT > $RUN/0_start/PARENT.txt   # run_chain_gadi.sh copies its final restart.dat and flow.vtu at job start
  else
    cp $START/restart.dat $START/flow.vtu $RUN/0_start/
    echo "start: copies of $START/{restart.dat,flow.vtu} ($(head -1 $START/../../title.txt 2>/dev/null))" > $RUN/0_start/SOURCE.txt
  fi
  echo "SST $lab: p0j $P bar, p0inf $Q bar, meshFine, SLAU2 + VW" > $RUN/title.txt
  printf "label %s\ntunnel p0 gauge %s bar -> p0inf %s Pa\nnozzle target %s bar -> p0j %s Pa (purple fit)\nstart %s\ncfd-tools %s\n" \
    $lab $g $p0 $t $pj "$START" $COMMIT > $RUN/BC.txt
  prev=../../0_start/restart
  for st in $STAGES; do
    IFS=: read n it c0 cap <<< "$st"
    mkdir -p $RUN/2_jet_staged/$n
    mkcfg $RUN/2_jet_staged/$n/$n.cfg $lab $pj $p0 $n $it $c0 $cap $prev "$START"
    prev=../$n/restart
  done
  echo "== $RUN"; cat $RUN/BC.txt; ls $RUN/2_jet_staged
  echo "-- diff vs the 6.6 bar t5 config (last stage, settings only):"
  diff $TEMPLATE $RUN/2_jet_staged/$n/$n.cfg | grep -E "^[<>] [A-Z]" || true
done
