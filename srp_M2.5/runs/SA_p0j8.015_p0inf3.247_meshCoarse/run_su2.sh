#!/bin/bash
# Run Phase 1 (freestream only) of case SA_p0j8.015_p0inf3.247 on mesh v4, on a Strudel desktop node (no sbatch).
# Trimmed on 2026-10-28 from the old ../SA_p0j8.015_p0inf3.247/run_su2.sh: the other stages belonged to deleted runs.
# Usage, on the desktop node: tmux new -d -s srp_v4_phase1 "bash <this folder>/run_su2.sh v4_phase1"
# Ranks = CPUs in this shell's allocation (24 on a desktop). Phase 2 runs use run_staged.sh / run_s5.sh / run_s6.sh.
R=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs/SA_p0j8.015_p0inf3.247_meshCoarse
STAGE=${1:?stage: v4_phase1}
NP=${NP:-$(nproc)}
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh
set -u   # only after sourcing: the module scripts use unset variables

run() {  # run <folder> <cfg>
  cd "$1" || exit 1
  local log=${2%.cfg}.log
  { echo "host $(hostname)  ranks $NP  cpus $(taskset -cp $$ | cut -d: -f2)"; date; } > "$log"
  mpirun -n "$NP" SU2_CFD -t 1 "$2" >> "$log" 2>&1
  echo "SU2_CFD exit=$?  $(date)" >> "$log"
}

case $STAGE in
  v4_phase1)  # geometry/mesh v4 (outlet at x = -308.1 mm, supersonic outlet), freestream only
          run $R/1_freestream phase1.cfg ;;
  *) echo "unknown stage $STAGE"; exit 1 ;;
esac
