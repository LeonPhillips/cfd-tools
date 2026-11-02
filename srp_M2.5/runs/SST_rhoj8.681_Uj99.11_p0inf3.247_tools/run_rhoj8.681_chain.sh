#!/bin/bash
# k-omega SST with a MASS_FLOW jet inlet (rho 8.681 kg/m3, V 99.11 m/s), p0inf 3.247 bar (user, 2026-11-02).
# Configs built by setup_rhoj8.681.sh.
#   meshCoarse from scratch: Phase 1 (20000, nozzle plugged), then t1 cap 0.3 / t2 cap 0.5 / t3 cap 1 / t4 cap 3, 50000 each.
#   meshFineJet: coarse t4 interpolated, then t1 cap 0.5 (5000), t2 cap 1 (25000), t3 cap 3 (20000).
#   Every stage runs its full length; the watchers only log and plot.
#   After each mesh: its line check (coarse t3 vs t4, fineJet t2 vs t3) and total-conditions inlet (SST_p0j7.4945, t3)
#   vs MASS_FLOW lines on the same mesh.
# Run on the desktop node from the run-tree copy of this folder, after the p0j 7.4945 chain:
#   SU2_CORES=<list> bash <this folder>/run_rhoj8.681_chain.sh   (log: chain.log, appended)
# SU2_CORES=<comma list>: one rank per listed core (Intel MPI pin list); keep the count even (interpolation).
# Restartable: a step counts as done only if its log says "SU2_CFD exit=0" and its history reaches ITER-1.
TOOLS=$(dirname $(readlink -f $0))
R=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs
C=$R/SA_p0j8.015_p0inf3.247_meshCoarse   # line_profiles.py, plot_lines.py, check_run.py
NEW=SST_rhoj8.681_Uj99.11_p0inf3.247; TC=SST_p0j7.4945_p0inf3.247
run() { echo $R/${NEW}_mesh$1; }
exec >> $TOOLS/chain.log 2>&1
if [ -n "$SU2_CORES" ]; then NP=$(echo $SU2_CORES | tr "," " " | wc -w); export I_MPI_PIN_PROCESSOR_LIST=$SU2_CORES; else NP=$(nproc); fi
echo "chain start $(date) on $(hostname), $NP ranks${SU2_CORES:+ on cores $SU2_CORES}"
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh

su2() {  # su2 <folder> <cfg>: run SU2_CFD there, log to <cfg>.log, return its exit code
  cd $1 || return 1
  { echo "host $(hostname)  ranks $NP${SU2_CORES:+  cores $SU2_CORES}"; date; } > ${2%.cfg}.log
  mpirun -n $NP SU2_CFD -t 1 $2 >> ${2%.cfg}.log 2>&1
  local rc=$?; echo "SU2_CFD exit=$rc  $(date)" >> ${2%.cfg}.log; return $rc
}
last_it() { tail -1 $1/history.csv 2>/dev/null | cut -d, -f3 | tr -d ' '; }
done_step() {  # done_step <folder> <cfg>: true if that SU2 run finished all its iterations
  local it n=$(sed -n 's/^ITER= *\([0-9]*\).*/\1/p' $1/$2)
  grep -q "SU2_CFD exit=0" $1/${2%.cfg}.log 2>/dev/null || return 1
  [ -s $1/restart.dat ] || return 1
  [ "$n" = 1 ] && return 0          # interpolation: 1 iteration
  it=$(last_it $1); [ -n "$it" ] && [ "$it" -eq $((n - 1)) ]
}
stage() {  # stage <run folder> <stage> <previous stage folder or ->
  local d=$1/2_jet_staged/$2 rc
  done_step $d $2.cfg && { echo "$(basename $1) $2 already complete: skipped"; return 0; }
  echo "$(basename $1) $2 start $(date), ITER $(sed -n 's/^ITER= *\([0-9]*\).*/\1/p' $d/$2.cfg)"
  (sleep 30; HOLD=0 PREV=$3 bash $TOOLS/watch_fj.sh $d) &
  su2 $d $2.cfg; rc=$?
  wait
  echo "$(basename $1) $2 end exit=$rc $(date), last iteration $(last_it $d)"
  done_step $d $2.cfg || { echo "STAGE FAILED: $(basename $1) $2. Chain stopped."; exit 1; }
}
lines() {  # lines <run folder> <mesh> <stage before last> <last stage>: into <run>/analysis
  local A=$1/analysis S=$1/2_jet_staged; mkdir -p $A; cd $A
  ( unset LD_PRELOAD; module load paraview/5.13 >/dev/null 2>&1
    pvbatch --mesa $C/line_profiles.py conv_lines.npz ${3%%_*}_end=$S/$3/flow.vtu ${4%%_*}_end=$S/$4/flow.vtu 2>&1 | grep -vi warning
    pvbatch --mesa $C/line_profiles.py inlet_lines.npz total_p0j7.4945=$R/${TC}_mesh$2/2_jet_staged/t3_cap3/flow.vtu \
      massflow=$S/$4/flow.vtu 2>&1 | grep -vi warning )
  echo "-- $(basename $1): end of $3 vs end of $4"
  python $C/plot_lines.py conv_lines.png conv_lines.npz "SST MASS_FLOW inlet, mesh$2: end of $3 vs end of $4"
  echo "-- $(basename $1): total-conditions inlet (p0j 7.4945 bar, t3) vs MASS_FLOW ($4) (differences against MASS_FLOW)"
  python $C/plot_lines.py inlet_lines.png inlet_lines.npz "SST mesh$2: total conditions p0j 7.4945 bar vs MASS_FLOW rho 8.681, V 99.11"
}

# ---- coarse from scratch ----
CO=$(run Coarse)
if ! done_step $CO/1_freestream phase1.cfg; then
  echo "Coarse phase1 start $(date)"
  (sleep 30; HOLD=0 bash $TOOLS/watch_fj.sh $CO/1_freestream) &
  su2 $CO/1_freestream phase1.cfg; rc=$?
  wait
  echo "Coarse phase1 end exit=$rc $(date), last iteration $(last_it $CO/1_freestream)"
  done_step $CO/1_freestream phase1.cfg || { echo "PHASE 1 FAILED. Chain stopped."; exit 1; }
fi
prev=-
for s in t1_cap0.3 t2_cap0.5 t3_cap1 t4_cap3; do stage $CO $s $prev; prev=$CO/2_jet_staged/$s; done
[ -s $CO/analysis/inlet_lines.png ] || lines $CO Coarse t3_cap1 t4_cap3

# ---- meshFineJet from coarse t4 ----
d=$(run FineJet)
if ! done_step $d/0_interp interp.cfg; then
  echo "FineJet interpolation start $(date)"
  su2 $d/0_interp interp.cfg; rc=$?
  grep -q "does not match the mesh" $d/0_interp/interp.log && echo "  (SU2 interpolated: point counts differ)"
  done_step $d/0_interp interp.cfg || { echo "INTERPOLATION FAILED: FineJet. Chain stopped."; exit 1; }
fi
prev=-
for s in t1_cap0.5 t2_cap1 t3_cap3; do stage $d $s $prev; prev=$d/2_jet_staged/$s; done
[ -s $d/analysis/inlet_lines.png ] || lines $d FineJet t2_cap1 t3_cap3
echo "ALL_DONE $(date)"
