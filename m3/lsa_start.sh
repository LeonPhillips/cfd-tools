#!/bin/bash
# Pull cfd-tools on M3, copy lsa_local.m into the LSA work folder stamped with the commit hash, and start one
# lsa_local MODE in tmux on the compute node of an existing job (default: the desktop job on m3f102).
# Usage (from the PC):  ssh ... 'MODE=sensitivity NODE=m3f102 bash -s' < scripts/m3/lsa_start.sh
set -euo pipefail
MODE=${MODE:?set MODE}; NODE=${NODE:-m3f102}; NWORKERS=${NWORKERS:-0}
TOOLS=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cfd-tools
WORK=/scratch/qk75/LeonPhillips/PhD/LSA

git -C "$TOOLS" pull -q --ff-only
HASH=$(git -C "$TOOLS" rev-parse --short HEAD)
{ echo "% Copied from cfd-tools lsa/lsa_local.m at commit $HASH on $(date +%F). Edit in cfd-tools, not here."
  cat "$TOOLS/lsa/lsa_local.m"; } > "$WORK/lsa_local.m"
cp "$TOOLS/lsa/README.md" "$WORK/README.md"
echo "copied lsa_local.m at $HASH to $WORK"

ssh -o BatchMode=yes "$NODE" "tmux new-session -d -s lsa_$MODE \"cd $WORK && module load matlab/r2023b && \
  matlab -batch \\\"MODE='$MODE'; NWORKERS=$NWORKERS; lsa_local\\\" > $WORK/run_$MODE.out 2>&1\""
echo "started MODE=$MODE in tmux lsa_$MODE on $NODE; output $WORK/run_$MODE.out and $WORK/$MODE/"
