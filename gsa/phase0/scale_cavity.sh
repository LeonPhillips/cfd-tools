#!/bin/bash
# Memory / time scaling of a LinStab2D GSA solve (direct + adjoint eigs, 50 modes, p = 300),
# using the cavity example refined by factors F. Run from the output folder:
#   bash scale_cavity.sh <LinStab2D repo> [factors...]
# Stops refining once a run peaks above MAXGB (default 60) of resident memory.
set -u
REPO=$(readlink -f "$1"); shift
FACTORS=${*:-"1 1.5 2 2.5"}
MAXGB=${MAXGB:-60}
NT=${NT:-${SLURM_CPUS_PER_TASK:-16}}
OUT=$(pwd)
module load matlab/r2023b
echo "factor,exit,dof,nnz,eigs_s,wall_s,maxrss_MB" > "$OUT/scale_summary.csv"
for F in $FACTORS; do
    ex="cav_F${F/./p}"
    sed -e "s/^Nx_cavity *=.*/Nx_cavity = round(80*$F);/" \
        -e "s/^Ny_cavity *=.*/Ny_cavity = round(50*$F);/" \
        -e "s/^verbose *=.*/verbose = false;/" \
        -e "s/^time_eigsMatlab = toc;/time_eigsMatlab = toc; fprintf('CALIB dof %d nnz %d eigs_s %.1f\\\\n', size(L0,1), nnz(L0), time_eigsMatlab);/" \
        "$REPO/example_4_Cavity.m" > "$OUT/$ex.m"
    /usr/bin/time -v -o "$OUT/$ex.time" matlab -batch \
        "cd('$REPO'); maxNumCompThreads($NT); run('$OUT/$ex.m')" > "$OUT/$ex.log" 2>&1
    rc=$?
    read dof nnzL eigs_s < <(awk '/^CALIB/ {print $3, $5, $7}' "$OUT/$ex.log")
    wall=$(awk -F': ' '/Elapsed \(wall/ {split($2,a,":"); s=0; for(i=1;i<=length(a);i++) s=s*60+a[i]; print s}' "$OUT/$ex.time")
    rss=$(awk -F': ' '/Maximum resident/ {printf "%.0f", $2/1024}' "$OUT/$ex.time")
    echo "$F,$rc,${dof:-NA},${nnzL:-NA},${eigs_s:-NA},$wall,$rss" | tee -a "$OUT/scale_summary.csv"
    [ "$rc" -ne 0 ] && break
    [ "$rss" -gt $((MAXGB*1024)) ] && { echo "peak above ${MAXGB} GB, stopping"; break; }
done
echo "done $(date +%T)"
