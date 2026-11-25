#!/bin/bash
# CT sweep pilot (user, 2026-11-25): two Schlieren1 conditions from cases/srp_M2.5/CT_sweep_conditions.md (Market),
# run with the 6.6 bar SLAU2-VW meshFine setup, restarted from that run's converged t5 solution (iteration 99,999).
# The file-name CT/CTa/CTd are labels only; the RANS results are not matched to them (user, 2026-11-25).
# Boundary conditions (confirmed, user 2026-11-25): p0inf = tunnel p0 absolute (gauge + 1.01325 bar);
#   p0j = nozzle p01 after losses = 0.81895 * target + 0.03094 bar (green fit, pitot_pdrop_close_p0_vs_target.fig).
#   T0inf = T0j = 288.15 K, M 2.5: T = 128.067 K, U = 567.16 m/s, p = p0inf / 2.25^3.5.
# Builds, under /scratch/os85/lp3540/su2/cases/srp_M2.5/runs/:
#   SST_p0j<P>_p0inf<Q>_meshFine_SLAU2-VW[_test]/
#     0_start/                  copies of the 6.6 bar SLAU2-VW t5 restart.dat and flow.vtu (same mesh)
#     2_jet_staged/t1_cap0.5 ... the stages (FULL: t1-t5 = 1k/1k/8k/40k/50k; TEST: t1 1k + t2 1k)
#     title.txt
# Each config is the 6.6 bar t5 config with only: the freestream pressure (FREESTREAM_PRESSURE and
# MARKER_SUPERSONIC_INLET), the jet p0 (MARKER_INLET), ITER, CFL start and cap, SOLUTION_FILENAME, comments.
# Usage (login node; copies ~160 MB per case):
#   TEST=1 bash setup_ct_pilot.sh        then, after the tests:   bash setup_ct_pilot.sh
#   CASES="label:p0j_Pa:p0inf_Pa ..." overrides the case list.
set -e
R=/scratch/os85/lp3540/su2/cases/srp_M2.5/runs
SRC=$R/SST_p0j6.6_p0inf3.247_meshFine_SLAU2-VW/2_jet_staged/t5_cap3
TOOLS=$(dirname $(readlink -f $0))
COMMIT=$(git -C $TOOLS rev-parse --short HEAD 2>/dev/null || echo unknown)
# label : p0j [Pa] : p0inf [Pa]
#   CT 0.1626: tunnel 5.90 bar gauge -> 6.91325 bar abs; target 5.50 -> p01 4.535165 bar
#   CT 0.5202: tunnel 2.00 bar gauge -> 3.01325 bar abs; target 9.00 -> p01 7.40149 bar
CASES=${CASES:-"CT0.1626:453516.5:691325 CT0.5202:740149.0:301325"}

mkcfg() {  # mkcfg <out> <label> <p0j Pa> <p0inf Pa> <stage> <ITER> <CFL0> <cap> <solution>
  local out=$1 lab=$2 pj=$3 p0=$4
  local T=128.067 U=567.16
  local p=$(awk -v a=$p0 'BEGIN{printf "%.1f", a/2.25^3.5}')
  local rho=$(awk -v p=$p 'BEGIN{printf "%.5f", p/(287.058*128.067)}')
  local pjb=$(awk -v a=$pj 'BEGIN{printf "%.4f", a/1e5}') p0b=$(awk -v a=$p0 'BEGIN{printf "%.4f", a/1e5}')
  {
    echo "% srp_M2.5 CT sweep pilot ($lab, Schlieren1), meshFine, k-omega SST, SLAU2 + Venkatakrishnan-Wang (user, 2026-11-25)."
    echo "% p0j $pjb bar (p01 after losses), p0inf $p0b bar abs. Stage $5: $6 iterations, CFL $7 -> cap $8 (x1.01), floor 0.1."
    echo "% Made by SST_SLAU2-VW_meshFine_tools/setup_ct_pilot.sh (cfd-tools $COMMIT) from $SRC/t5_cap3.cfg"
    echo "% (only the freestream pressure, jet p0, ITER, CFL and solution path changed). Nothing stops on convergence."
    sed -n '/^% ---- physics/,$p' $SRC/t5_cap3.cfg
  } | sed -e "s|^% ---- freestream: .*|% ---- freestream: M2.5 from p0inf = $p0b bar abs, T0 = 288.15 K (isentropic; user 2026-11-25) ----|" \
          -e "s|^% T = 288.15/2.25 .*|% T = 288.15/2.25 = $T K; p = $p0/2.25^3.5 = $p Pa; U = $U m/s in -x|" \
          -e "s|^% (AOA = 180 deg); rho = .*|% (AOA = 180 deg); rho = $rho kg/m3. Restarted from the 6.6 bar / 3.247 bar SLAU2-VW t5 solution.|" \
          -e "s|^FREESTREAM_PRESSURE=.*|FREESTREAM_PRESSURE= $p|" \
          -e "s|^MARKER_SUPERSONIC_INLET=.*|MARKER_SUPERSONIC_INLET= ( inlet_freestream, $T, $p, -$U, 0.0, 0.0 )|" \
          -e "s|^% Phase 2: inlet_jet .*|% inlet_jet: total conditions, p0j = $pjb bar (p01 after losses), T0j = 288.15 K, flow in +x (user 2026-11-25).|" \
          -e "s|^MARKER_INLET=.*|MARKER_INLET= ( inlet_jet, 288.15, $pj, 1.0, 0.0, 0.0 )   % p0j $pjb bar (total), T0j 288.15 K (user, 2026-11-25)|" \
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

[ -s $SRC/restart.dat ] && [ -s $SRC/flow.vtu ] && [ -s $SRC/t5_cap3.cfg ] || { echo "missing source files in $SRC"; exit 1; }
if [ "${TEST:-0}" = 1 ]; then
  STAGES="t1_cap0.5:1000:0.1:0.5 t2_cap1:1000:0.5:1.0"; SUF=_test
else
  STAGES="t1_cap0.5:1000:0.1:0.5 t2_cap1:1000:0.5:1.0 t3_cap3:8000:1.0:3.0 t4_cap3:40000:1.0:3.0 t5_cap3:50000:1.0:3.0"; SUF=
fi
for c in $CASES; do
  IFS=: read lab pj p0 <<< "$c"
  P=$(awk -v a=$pj 'BEGIN{printf "%.3f", a/1e5}'); Q=$(awk -v a=$p0 'BEGIN{printf "%.3f", a/1e5}')
  RUN=$R/SST_p0j${P}_p0inf${Q}_meshFine_SLAU2-VW$SUF
  [ -e $RUN ] && { echo "$RUN exists: not touched"; exit 1; }
  mkdir -p $RUN/0_start
  cp $SRC/restart.dat $SRC/flow.vtu $RUN/0_start/
  echo "start: copies of $SRC/{restart.dat,flow.vtu} (p0j 6.6 / p0inf 3.247 bar SLAU2-VW t5, cumulative iteration 99999)" > $RUN/0_start/SOURCE.txt
  echo "${SUF:+TEST }SST $lab: p0j $P bar, p0inf $Q bar, meshFine, SLAU2 + VW" > $RUN/title.txt
  prev=../../0_start/restart
  for st in $STAGES; do
    IFS=: read n it c0 cap <<< "$st"
    mkdir -p $RUN/2_jet_staged/$n
    mkcfg $RUN/2_jet_staged/$n/$n.cfg $lab $pj $p0 $n $it $c0 $cap $prev
    prev=../$n/restart
  done
  echo "== $RUN"; ls $RUN/2_jet_staged
  echo "-- diff vs the 6.6 bar t5 config (last stage, settings only):"
  diff $SRC/t5_cap3.cfg $RUN/2_jet_staged/$n/$n.cfg | grep -E "^[<>] [A-Z]" || true
done
