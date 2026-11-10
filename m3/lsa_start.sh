#!/bin/bash
# Pull cfd-tools on M3, copy lsa_local.m into the LSA work folder stamped with the commit hash, and start one
# lsa_local MODE as an srun step inside an existing job (default: the desktop job), held in tmux on the login node.
# ssh straight to the compute node fails host-key verification (2026-11-10), so the step runs with srun --overlap.
# srun counts as submitting: ask the user first (M3 rules).
# Usage (from the PC):  ssh ... 'MODE=sensitivity JOBID=62391826 CPUS=8 MEM=32G bash -s' < scripts/m3/lsa_start.sh
set -euo pipefail
MODE=${MODE:?set MODE}; JOBID=${JOBID:?set JOBID}; CPUS=${CPUS:-8}; MEM=${MEM:-32G}; NWORKERS=${NWORKERS:-0}
TOOLS=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cfd-tools
WORK=/scratch/qk75/LeonPhillips/PhD/LSA

git -C "$TOOLS" pull -q --ff-only
HASH=$(git -C "$TOOLS" log -1 --format=%h -- lsa/lsa_local.m)     # last commit that changed lsa_local.m
{ echo "% Copied from cfd-tools lsa/lsa_local.m at commit $HASH on $(date +%F). Edit in cfd-tools, not here."
  cat "$TOOLS/lsa/lsa_local.m"; } > "$WORK/lsa_local.m"
cp "$TOOLS/lsa/README.md" "$WORK/README.md"
echo "copied lsa_local.m at $HASH to $WORK"

tmux new-session -d -s "lsa_$MODE" "srun --jobid=$JOBID --overlap --ntasks=1 --cpus-per-task=$CPUS --mem=$MEM \
  bash -lc 'cd $WORK && module load matlab/r2023b && matlab -batch \"MODE='\''$MODE'\''; NWORKERS=$NWORKERS; lsa_local\"' \
  > $WORK/run_$MODE.out 2>&1"
echo "started MODE=$MODE as a step of job $JOBID (tmux lsa_$MODE on $(hostname)); output $WORK/run_$MODE.out, $WORK/$MODE/"
