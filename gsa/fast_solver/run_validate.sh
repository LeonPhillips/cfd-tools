#!/bin/bash
# Validate the MUMPS shift-invert wrapper on LinStab2D's cavity example, then time it at larger size.
# Run from the output folder (on a compute node you hold, or inside a job):
#   bash run_validate.sh <LinStab2D repo> <MUMPS dir with zmumpsmex.mexa64>
#   NT=16 (threads), CASES="1:builtin,umfpack,mumps 2.5:mumps" (factor:solvers ...)
set -u
REPO=$(readlink -f "$1"); MUMPS=$(readlink -f "$2")
NT=${NT:-${SLURM_CPUS_PER_TASK:-16}}
CASES=${CASES:-"1:builtin,umfpack,mumps 2.5:mumps"}
HERE=${FLS_DIR:-/scratch/qk75/LeonPhillips/PhD/cfd/su2/cfd-tools/gsa/fast_solver}
OUT=$(pwd)
module load ${MATLAB_MODULES:-matlab/r2023b}
export OMP_NUM_THREADS=$NT FLS_MUMPS_DIR=$MUMPS
echo "host $(hostname)  threads $NT  repo $(git -C "$REPO" rev-parse --short HEAD)  tools $(git -C "$HERE" rev-parse --short HEAD)"
echo "case,exit,wall_s,maxrss_MB" > "$OUT/validate_summary.csv"
for c in $CASES; do
    F=${c%%:*}; M=${c#*:}
    tag="F${F/./p}_${M//,/-}"
    mlist="{'$(echo "$M" | sed "s/,/','/g")'}"
    /usr/bin/time -v -o "$OUT/$tag.time" matlab -batch \
        "addpath('$HERE'); maxNumCompThreads($NT); validate_fast_solver_cavity('$REPO','$OUT',$F,$mlist)" \
        > "$OUT/$tag.log" 2>&1
    rc=$?
    wall=$(awk -F': ' '/Elapsed \(wall/ {split($2,a,":"); s=0; for(i=1;i<=length(a);i++) s=s*60+a[i]; print s}' "$OUT/$tag.time")
    rss=$(awk -F': ' '/Maximum resident/ {printf "%.0f", $2/1024}' "$OUT/$tag.time")
    echo "$tag,$rc,$wall,$rss" | tee -a "$OUT/validate_summary.csv"
    grep -v Fontconfig "$OUT/$tag.log" | tail -12
done
echo "done $(date +%T)"
