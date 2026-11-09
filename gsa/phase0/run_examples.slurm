#!/bin/bash
# Phase 0: run LinStab2D's own examples on M3 to check the install and size later jobs.
# Usage (from the output folder):  sbatch <cfd-tools>/gsa/phase0/run_examples.slurm <LinStab2D repo> [examples...]
# Each example runs in its own MATLAB process, so its wall time and peak memory are logged separately.
#SBATCH --job-name=ls2d_phase0
#SBATCH --account=qk75
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=phase0_%j.out

set -u
REPO=$(readlink -f "$1"); shift
EXAMPLES=${*:-"example_0_TestCurvilinearDerivative example_2_Jet example_4_Cavity example_5_Paraboloid"}
OUT=$(pwd)
# sbatch runs a spooled copy of this script, so point at the M3 clone explicitly
TOOLS=${CFD_TOOLS_GSA:-/scratch/qk75/LeonPhillips/PhD/cfd/su2/cfd-tools/gsa/phase0}

module load matlab/r2023b
echo "host $(hostname)  cpus $SLURM_CPUS_PER_TASK  repo $REPO ($(git -C "$REPO" rev-parse --short HEAD))"
echo "example,exit,wall_s,maxrss_MB" > "$OUT/summary.csv"

for ex in $EXAMPLES; do
    echo "=================== $ex  $(date +%T)"
    /usr/bin/time -v -o "$OUT/${ex}.time" \
        matlab -batch "cd('$REPO'); addpath('$TOOLS'); maxNumCompThreads($SLURM_CPUS_PER_TASK); run('${ex}.m'); save_all_figs('$OUT','$ex')" \
        > "$OUT/${ex}.log" 2>&1
    rc=$?
    wall=$(awk -F': ' '/Elapsed \(wall/ {split($2,a,":"); s=0; for(i=1;i<=length(a);i++) s=s*60+a[i]; print s}' "$OUT/${ex}.time")
    rss=$(awk -F': ' '/Maximum resident/ {printf "%.0f", $2/1024}' "$OUT/${ex}.time")
    echo "$ex,$rc,$wall,$rss" | tee -a "$OUT/summary.csv"
    tail -5 "$OUT/${ex}.log"
done
echo "done $(date +%T)"
