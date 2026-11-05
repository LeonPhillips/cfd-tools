#!/bin/bash
# Add a continuation stage to a SLAU2-VW run folder (user, 2026-11-05: +40k at cap 3 after t3).
# The new config is the previous stage's config with only ITER, CFL_NUMBER, the CFL cap and SOLUTION_FILENAME changed
# (restart from the previous stage). Probes are reset by run_chain_gadi.sh, which skips completed stages.
# Usage (login node): bash add_stage.sh <run folder> <previous stage> <new stage> <ITER> <CFL start> <CFL cap>
#   e.g. bash add_stage.sh $RUN t3_cap3 t4_cap3 40000 1.0 3.0
set -e
RUN=$(readlink -f ${1:?run folder}); P=${2:?previous stage}; N=${3:?new stage}; IT=${4:?ITER}; C0=${5:?CFL start}; CAP=${6:?CFL cap}
S=$RUN/2_jet_staged/$P; D=$RUN/2_jet_staged/$N
[ -s $S/restart.dat ] || { echo "no $S/restart.dat"; exit 1; }
[ -e $D ] && { echo "$D exists: not touched"; exit 1; }
mkdir -p $D
{
  echo "% Stage $N (continuation of $P, user 2026-11-05): $IT iterations, CFL starts at $C0, cap $CAP. Made by add_stage.sh"
  echo "% (cfd-tools $(git -C $(dirname $(readlink -f $0)) rev-parse --short HEAD)) from ../$P/$P.cfg; only ITER, CFL and the solution file changed."
  cat $S/$P.cfg
} | sed -e "s|^ITER=.*|ITER= $IT|" -e "s|^CFL_NUMBER=.*|CFL_NUMBER= $C0|" \
        -e "s|^CFL_ADAPT_PARAM=.*|CFL_ADAPT_PARAM= ( 0.5, 1.01, 0.1, $CAP, 0.001, 0 )|" \
        -e "s|^SOLUTION_FILENAME=.*|SOLUTION_FILENAME= ../$P/restart|" > $D/$N.cfg
echo "== $D"; diff $S/$P.cfg $D/$N.cfg | grep -E "^[<>] [A-Z]" || true
