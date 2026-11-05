#!/bin/bash
# SLAU2 + Venkatakrishnan-Wang on meshFine, p0j 7.0 and 7.4945 bar, on Gadi (user, 2026-11-05). Builds the run folders:
#   <R>/SST_p0j<P>_p0inf3.247_meshFine_SLAU2-VW/
#     0_start/                 copies (not symlinks) of the AUSM+UP2 t3 restart.dat and flow.vtu on the same mesh
#     2_jet_staged/t1_cap0.5/  1000 iterations, CFL 0.1 -> cap 0.5
#     2_jet_staged/t2_cap1/    1000 iterations, CFL 0.5 -> cap 1
#     2_jet_staged/t3_cap3/    8000 iterations, CFL 1.0 -> cap 3      (10,000 in total)
# Each config is the AUSM+UP2 t3 config of that case with only these changes:
#   CONV_NUM_METHOD_FLOW= SLAU2, SLOPE_LIMITER_FLOW= VENKATAKRISHNAN_WANG (VENKAT_LIMITER_COEFF 0.05 kept),
#   ITER / CFL_NUMBER / CFL_ADAPT_PARAM cap per stage, SOLUTION_FILENAME, MESH_FILENAME (Gadi path), the header.
#   CUSTOM_OUTPUTS is reset per stage by the chain driver (probes_md.py). Nothing stops on convergence.
# TEST=1: builds <R>/SST_p0j7.4945_p0inf3.247_meshFine_SLAU2-VW_test with t1_cap0.5 only, ITER 2000 (the test job).
# Usage (login node; copies ~160 MB per case): bash setup_slau2vw.sh            or   TEST=1 bash setup_slau2vw.sh
set -e
R=/scratch/os85/lp3540/su2/cases/srp_M2.5/runs
MESH=/scratch/os85/lp3540/su2/cases/srp_M2.5/mesh/srp_M2.5_meshFine.su2
TOOLS=$(dirname $(readlink -f $0))
COMMIT=$(git -C $TOOLS rev-parse --short HEAD 2>/dev/null || echo unknown)

mkcfg() {  # mkcfg <source cfg> <out cfg> <P> <stage> <ITER> <CFL start> <CFL cap> <solution file>
  local src=$1 out=$2
  {
    echo "% srp_M2.5, case SST_p0j$3_p0inf3.247, meshFine, k-omega SST, SLAU2 + Venkatakrishnan-Wang (user, 2026-11-05)."
    echo "% Stage $4: $5 iterations, CFL starts at $6, ramps x1.01 per iteration to a cap of $7, halved on trouble, floor 0.1."
    echo "% Made by SST_SLAU2-VW_meshFine_tools/setup_slau2vw.sh (cfd-tools $COMMIT) from $src"
    echo "% (only the scheme, limiter, ITER, CFL, solution and mesh paths changed). Nothing stops on convergence."
    sed -n '/^% ---- physics/,$p' $src
  } | sed -e "s|^CONV_NUM_METHOD_FLOW=.*|CONV_NUM_METHOD_FLOW= SLAU2   % was AUSMPLUSUP2 (user, 2026-11-05)|" \
          -e "s|^SLOPE_LIMITER_FLOW=.*|SLOPE_LIMITER_FLOW= VENKATAKRISHNAN_WANG   % was VENKATAKRISHNAN (user, 2026-11-05)|" \
          -e "s|^ITER=.*|ITER= $5|" \
          -e "s|^CFL_NUMBER=.*|CFL_NUMBER= $6|" \
          -e "s|^CFL_ADAPT_PARAM=.*|CFL_ADAPT_PARAM= ( 0.5, 1.01, 0.1, $7, 0.001, 0 )|" \
          -e "s|^SOLUTION_FILENAME=.*|SOLUTION_FILENAME= $8|" \
          -e "s|^MESH_FILENAME=.*|MESH_FILENAME= $MESH|" > $out
  for k in "CONV_NUM_METHOD_FLOW= SLAU2" "SLOPE_LIMITER_FLOW= VENKATAKRISHNAN_WANG" "VENKAT_LIMITER_COEFF= 0.05" \
           "ITER= $5" "CFL_NUMBER= $6" "CFL_ADAPT= YES" "CONV_RESIDUAL_MINVAL= -20" "WRT_VOLUME_OVERWRITE= NO" \
           "OUTPUT_WRT_FREQ= 1000, 1000, 1000" "MUSCL_FLOW= YES" "RESTART_SOL= YES"; do
    grep -q "^$k" $out || { echo "CHECK FAILED: '$k' missing in $out"; exit 1; }
  done
  grep -q "^HISTORY_OUTPUT=.*CUSTOM" $out || { echo "CHECK FAILED: CUSTOM not in HISTORY_OUTPUT ($out)"; exit 1; }
}

build() {  # build <P> <run folder> <stage list "name:ITER:CFL0:cap ...">
  local P=$1 RUN=$2 S=$R/SST_p0j$1_p0inf3.247_meshFine/2_jet_staged/t3_cap3 prev=../../0_start/restart st n it c0 cap
  [ -s $S/restart.dat ] && [ -s $S/flow.vtu ] && [ -s $S/t3_cap3.cfg ] || { echo "missing source files in $S"; exit 1; }
  [ -e $RUN ] && { echo "$RUN exists: not touched"; exit 1; }
  mkdir -p $RUN/0_start
  cp $S/restart.dat $S/flow.vtu $RUN/0_start/
  echo "start: copies of $S/{restart.dat,flow.vtu} (AUSM+UP2 t3, iteration 19999)" > $RUN/0_start/SOURCE.txt
  for st in $3; do
    IFS=: read n it c0 cap <<< "$st"
    mkdir -p $RUN/2_jet_staged/$n
    mkcfg $S/t3_cap3.cfg $RUN/2_jet_staged/$n/$n.cfg $P $n $it $c0 $cap $prev
    prev=../$n/restart
  done
  echo "== $RUN"; ls $RUN/2_jet_staged
  echo "-- diff vs the AUSM+UP2 t3 config (last stage):"
  diff $S/t3_cap3.cfg $RUN/2_jet_staged/$n/$n.cfg | grep -E "^[<>] [A-Z]" || true
}

if [ "${TEST:-0}" = 1 ]; then
  build 7.4945 $R/SST_p0j7.4945_p0inf3.247_meshFine_SLAU2-VW_test "t1_cap0.5:2000:0.1:0.5"
else
  for P in 7.0 7.4945; do
    build $P $R/SST_p0j${P}_p0inf3.247_meshFine_SLAU2-VW "t1_cap0.5:1000:0.1:0.5 t2_cap1:1000:0.5:1.0 t3_cap3:8000:1.0:3.0"
  done
fi
