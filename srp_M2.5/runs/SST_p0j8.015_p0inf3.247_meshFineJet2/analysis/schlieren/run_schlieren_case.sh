#!/bin/bash
# Synthetic schlieren for any run vs the NPR 501 mean experimental image, using the scripts in this folder (cfd-tools clone).
# Usage (on a compute/desktop node): run_schlieren_case.sh <run dir> <stage, e.g. 2_jet_staged/t3_cap3> <prefix> "<CFD label>" "<case label>"
# Output goes to <run dir>/analysis/schlieren/.
T=$(cd "$(dirname "$0")" && pwd)
RUN=$1; STAGE=$2; PRE=$3
export CFD_LABEL=$4 CASE_LABEL=$5
EXP=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs/SST_p0j8.015_p0inf3.247_meshFineJet2/analysis/schlieren/SCH1_ZOOM1_NPR501_220_800_1_mean.mat
A=$RUN/analysis/schlieren
mkdir -p $A && cd $A || exit 1
( unset LD_PRELOAD; module load paraview/5.13 >/dev/null 2>&1
  [ -f cfd_drhodx.npz ] || pvbatch --mesa $T/sample_schlieren.py $RUN/$STAGE/flow.vtu cfd_drhodx.npz 2>&1 | grep -vi warning )
module load anaconda/2024.02-1 >/dev/null 2>&1
python $T/schlieren_compare.py cfd_drhodx.npz $EXP $PRE
echo "cfd-tools $(git -C $T rev-parse --short HEAD), stage $STAGE" > source.txt
ls -l $A
