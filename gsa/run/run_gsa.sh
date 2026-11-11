#!/bin/bash
# Run run_gsa_srp.m for one or more grid spacings, each in its own MATLAB process (separate peak-memory log).
#   bash run_gsa.sh <case dir> <out prefix> "<h list>" "<m list>" "<shift list (MATLAB)>" [extra MATLAB args]
# e.g. bash run_gsa.sh $C gridcheck2 "0.1 0.05" "[0 1]" "2+0.1i"
# Writes <case dir>/<out prefix>_h<h>/{run.log,run.time,gsa_*.mat,*.png} and <case dir>/<out prefix>.done
set -u
C=$(readlink -f "$1"); PFX=$2; HS=$3; MS=$4; SH=$5; EXTRA=${6:-}
# cluster paths (defaults: M3); on Gadi set GSA_ROOT, GSA_TOOLS, MUMPS_DIR, MATLAB_MODULES (see run_gsa_gadi.pbs)
G=${GSA_ROOT:-/scratch/qk75/LeonPhillips/PhD/GSA}
T=${GSA_TOOLS:-/scratch/qk75/LeonPhillips/PhD/cfd/su2/cfd-tools/gsa}
NT=${NT:-16}
module load ${MATLAB_MODULES:-matlab/r2023b} >/dev/null 2>&1
export OMP_NUM_THREADS=$NT FLS_MUMPS_DIR=${MUMPS_DIR:-$G/mumps_m3}
rm -f "$C/$PFX.done"
for h in $HS; do
    d="$C/${PFX}_h${h/./p}"; mkdir -p "$d"; cd "$d"
    echo "h = $h  m = $MS  shifts = $SH  cfd-tools $(git -C $T rev-parse --short HEAD)  $(date)" > run.log
    /usr/bin/time -v -o run.time matlab -batch "addpath(genpath('$G/LinStab2D/aux_matlab')); addpath('$T/baseflow','$T/operator','$T/fast_solver','$T/run'); maxNumCompThreads($NT); try; run_gsa_srp(pwd, $h, $MS, $SH $EXTRA); catch e; disp(getReport(e)); end" >> run.log 2>&1
done
echo ALLDONE > "$C/$PFX.done"
