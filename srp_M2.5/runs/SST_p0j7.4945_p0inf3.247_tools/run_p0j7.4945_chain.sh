#!/bin/bash
# k-omega SST chain for srp_M2.5 with the corrected jet inlet, p0j 7.4945 bar / p0inf 3.247 bar (user, 2026-11-02).
# Configs built by setup_p0j7.4945.sh. Coarse restarts from its own SST 8.015 bar final solution; every other mesh starts
# from the previous, coarser mesh's new t3 solution, interpolated (0_interp; user, 2026-11-03): medium <- coarse,
# fine <- medium, fineJet <- fine, fineJet2 <- fineJet, veryFine <- fineJet2. Then t1 cap 0.5 (5000), t2 cap 1 (25000),
# t3 cap 3 (20000). Order: coarse, medium, fine, fineJet, fineJet2, veryFine. Every stage runs
# its full length; the watchers only log and plot.
#   After each mesh: its line check (t2 vs t3) and old inlet (8.015 bar, final) vs new (t3) lines.
#   At the end: grid lines across the six meshes and the 0.025 mm box resample of the final solutions.
# Run on the desktop node from the run-tree copy of this folder:
#   SU2_CORES=<list> tmux new -d -s srp_p0j7 "bash <this folder>/run_p0j7.4945_chain.sh"   (log: chain.log, appended)
# SU2_CORES=<comma list>: one rank per listed core (Intel MPI pin list); keep the count even. See the SST-Sarkar speed note.
# Restartable: a stage counts as done only if its log says "SU2_CFD exit=0" and its history reaches ITER-1.
TOOLS=$(dirname $(readlink -f $0))
R=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs
C=$R/SA_p0j8.015_p0inf3.247_meshCoarse   # line_profiles.py, plot_lines.py, check_run.py
NEW=SST_p0j7.4945_p0inf3.247; OLD=SST_p0j8.015_p0inf3.247
MESHES="Coarse Medium Fine FineJet FineJet2 VeryFine"
run() { echo $R/${NEW}_mesh$1; }
old_final() { [ $1 = Coarse ] && echo $R/${OLD}_meshCoarse/2_jet_staged/t4_cap3 || echo $R/${OLD}_mesh$1/2_jet_staged/t3_cap3; }
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
lines() {  # lines <run folder> <mesh>: t2 vs t3, and old vs new inlet, into <run>/analysis
  local A=$1/analysis S=$1/2_jet_staged; mkdir -p $A; cd $A
  ( unset LD_PRELOAD; module load paraview/5.13 >/dev/null 2>&1
    pvbatch --mesa $C/line_profiles.py conv_lines.npz t2_end=$S/t2_cap1/flow.vtu t3_end=$S/t3_cap3/flow.vtu 2>&1 | grep -vi warning
    pvbatch --mesa $C/line_profiles.py inlet_lines.npz p0j8.015=$(old_final $2)/flow.vtu p0j7.4945=$S/t3_cap3/flow.vtu 2>&1 | grep -vi warning )
  echo "-- $(basename $1): end of t2 vs end of t3"
  python $C/plot_lines.py conv_lines.png conv_lines.npz "SST p0j 7.4945 bar, mesh$2: end of t2 vs end of t3"
  echo "-- $(basename $1): old inlet (p0j 8.015 bar) vs new (7.4945 bar) (differences against the new one)"
  python $C/plot_lines.py inlet_lines.png inlet_lines.npz "SST mesh$2: p0j 8.015 bar (final) vs 7.4945 bar (t3)"
}

for m in $MESHES; do
  d=$(run $m); prev=-
  if [ $m != Coarse ] && ! done_step $d/0_interp interp.cfg; then
    echo "$m interpolation start $(date): $(sed -n 's/^SOLUTION_FILENAME= *//p' $d/0_interp/interp.cfg)"
    su2 $d/0_interp interp.cfg
    grep -q "does not match the mesh" $d/0_interp/interp.log && echo "  (SU2 interpolated: point counts differ)"
    done_step $d/0_interp interp.cfg || { echo "INTERPOLATION FAILED: $m. Chain stopped."; exit 1; }
  fi
  for s in t1_cap0.5 t2_cap1 t3_cap3; do stage $d $s $prev; prev=$d/2_jet_staged/$s; done
  [ -s $d/analysis/inlet_lines.png ] || lines $d $m
done
echo "all meshes done $(date)"

# ---- grid study and the box data for the contour pages ----
G=$R/grid_compare; cd $G
( unset LD_PRELOAD; module load paraview/5.13 >/dev/null 2>&1
  pvbatch --mesa $C/line_profiles.py sst_p0j7.4945_grid_lines.npz \
    $(for m in $MESHES; do echo "${m,}=$(run $m)/2_jet_staged/t3_cap3/flow.vtu"; done) 2>&1 | grep -vi warning
  H_MM=0.025 pvbatch --mesa $G/resample_box.py box_SST_p0j7.4945_0.025mm.npz \
    $(for m in $MESHES; do echo "SSTp0j7.4945_${m,}=$(run $m)/2_jet_staged/t3_cap3/flow.vtu"; done) 2>&1 | grep -vi warning )
echo "-- SST p0j 7.4945 bar grid study (differences against veryFine)"
python $C/plot_lines.py sst_p0j7.4945_grid_lines.png sst_p0j7.4945_grid_lines.npz "SST p0j 7.4945 bar final (t3) solutions: coarse to veryFine"
ls -l $G/sst_p0j7.4945_grid_lines.* $G/box_SST_p0j7.4945_0.025mm.npz
echo "ALL_DONE $(date)"
