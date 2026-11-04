#!/bin/bash
# k-omega SST chain for srp_M2.5 at jet total pressure p0j 7.0 bar / p0inf 3.247 bar (user, 2026-11-04).
# Configs built by setup_p0j7.0.sh. Coarse restarts from the 7.4945 bar coarse final solution; medium <- coarse,
# fine <- medium, veryFine <- fine, interpolated (0_interp). Then t1 cap 0.5 (5000), t2 cap 1 (25000), t3 cap 3 (20000).
# Every stage runs its full length. Convergence: the live monitor (monitor/watch_conv.sh, every 1000 iterations); the
# Mach-disk probes are re-placed before each stage from its starting solution (monitor/probes_md.py). Nothing stops.
#   After each mesh: Mach/density on the monitor lines, 7.4945 bar t3 vs 7.0 bar t3 on that mesh (veryFine: the SST
#   8.015 bar veryFine t3, as there is no 7.4945 bar veryFine run). No t2-vs-t3 comparison (user, 2026-11-04).
#   At the end: grid_lines.png and grid_features.png, the monitor's line and feature plots with one curve / point per
#   mesh (final t3 snapshot of each), into ../grid_compare/sst_p0j7.0_grid_*.
# Run on the desktop node from the run-tree copy of this folder, e.g. from the login node:
#   nohup srun --jobid=<desktop job> --overlap -n1 -c24 bash <this folder>/run_p0j7.0_chain.sh > /dev/null 2>&1 &
# (log: chain.log, appended). NP ranks (default 22). SU2_CORES=<comma list> pins one rank per listed core instead.
# Restartable: a stage counts as done only if its log says "SU2_CFD exit=0" and its history reaches ITER-1.
TOOLS=$(dirname $(readlink -f $0)); MON=$TOOLS/monitor
R=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs
C=$R/SA_p0j8.015_p0inf3.247_meshCoarse   # line_profiles.py, plot_lines.py, check_run.py
NEW=SST_p0j7.0_p0inf3.247; OLD=SST_p0j7.4945_p0inf3.247
MESHES=${MESHES:-"Coarse Medium Fine VeryFine"}
LINES="1,0.01,20,0.01;1,1.3,20,1.3;1,2.3,20,2.3;1,3,20,3;1,5,20,5;2.0,2.34,7.84,4.70"   # the monitor lines, mm
run() { echo $R/${NEW}_mesh$1; }
ref_final() { [ $1 = VeryFine ] && echo $R/SST_p0j8.015_p0inf3.247_meshVeryFine/2_jet_staged/t3_cap3 || echo $R/${OLD}_mesh$1/2_jet_staged/t3_cap3; }
ref_label() { [ $1 = VeryFine ] && echo p0j8.015 || echo p0j7.4945; }
exec >> $TOOLS/chain.log 2>&1
if [ -n "$SU2_CORES" ]; then NP=$(echo $SU2_CORES | tr "," " " | wc -w); export I_MPI_PIN_PROCESSOR_LIST=$SU2_CORES; else NP=${NP:-22}; fi
[ -n "$SLURM_STEP_ID" ] && export I_MPI_HYDRA_BOOTSTRAP=fork   # inside an srun step: start the ranks locally
echo "chain start $(date) on $(hostname), $NP ranks${SU2_CORES:+ on cores $SU2_CORES}, cpus $(grep Cpus_allowed_list /proc/self/status | cut -f2)"
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
probes() {  # probes <stage folder> <cfg> <starting flow.vtu>: CUSTOM_OUTPUTS 1 mm either side of the Mach disk
  local p
  p=$(cd $1; unset LD_PRELOAD; module load paraview/5.13 >/dev/null 2>&1
      pvbatch --mesa $MON/probes_md.py $3 2> probes.log | grep -v -i warning | tail -1)
  echo "  probes from $3:"; sed 's/^/    /' $1/probes.log
  [[ "$p" == up0* ]] || { echo "  WARNING: probes_md.py failed; the config keeps its previous probes"; return 1; }
  sed -i "s|^CUSTOM_OUTPUTS=.*|CUSTOM_OUTPUTS= '$p'|" $1/$2
}
stage() {  # stage <run folder> <stage> <previous stage folder or -> <starting flow.vtu>
  local d=$1/2_jet_staged/$2 rc
  done_step $d $2.cfg && { echo "$(basename $1) $2 already complete: skipped"; return 0; }
  echo "$(basename $1) $2 start $(date), ITER $(sed -n 's/^ITER= *\([0-9]*\).*/\1/p' $d/$2.cfg)"
  probes $d $2.cfg $4
  rm -f $d/flow_*.vtu   # leftovers of an interrupted attempt
  (sleep 30; HOLD=0 PREV=$3 bash $MON/watch_conv.sh $d) &
  su2 $d $2.cfg; rc=$?
  wait
  echo "$(basename $1) $2 end exit=$rc $(date), last iteration $(last_it $d), jet mass flow $(python - $d/history.csv <<'PY'
import csv, sys
r = list(csv.reader(open(sys.argv[1]))); h = [c.strip().strip('"') for c in r[0]]
print(r[-1][h.index("Avg_Massflow(inlet_jet)")].strip() if "Avg_Massflow(inlet_jet)" in h else "n/a")
PY
) kg/s"
  done_step $d $2.cfg || { echo "STAGE FAILED: $(basename $1) $2. Chain stopped."; exit 1; }
}
inlet_lines() {  # inlet_lines <run folder> <mesh>: reference-inlet t3 vs 7.0 bar t3 on the monitor lines, into <run>/analysis
  local A=$1/analysis S=$1/2_jet_staged; mkdir -p $A; cd $A
  ( unset LD_PRELOAD; module load paraview/5.13 >/dev/null 2>&1
    LINES=$LINES pvbatch --mesa $C/line_profiles.py inlet_lines.npz $(ref_label $2)=$(ref_final $2)/flow.vtu \
      p0j7.0=$S/t3_cap3/flow.vtu 2>&1 | grep -vi warning )
  echo "-- $(basename $1): $(ref_label $2) (final) vs p0j 7.0 bar (t3) (differences against 7.0 bar)"
  python $C/plot_lines.py inlet_lines.png inlet_lines.npz "SST mesh$2: $(ref_label $2) final vs p0j 7.0 bar t3"
}

for m in $MESHES; do
  d=$(run $m)
  if [ $m = Coarse ]; then start=$R/${OLD}_meshCoarse/2_jet_staged/t3_cap3/flow.vtu
  else
    if ! done_step $d/0_interp interp.cfg; then
      echo "$m interpolation start $(date): $(sed -n 's/^SOLUTION_FILENAME= *//p' $d/0_interp/interp.cfg)"
      su2 $d/0_interp interp.cfg
      grep -q "does not match the mesh" $d/0_interp/interp.log && echo "  (SU2 interpolated: point counts differ)"
      done_step $d/0_interp interp.cfg || { echo "INTERPOLATION FAILED: $m. Chain stopped."; exit 1; }
    fi
    start=$d/0_interp/flow.vtu
  fi
  prev=-
  for s in t1_cap0.5 t2_cap1 t3_cap3; do
    stage $d $s $prev $start; prev=$d/2_jet_staged/$s; start=$prev/flow.vtu
  done
  [ -s $d/analysis/inlet_lines.png ] || inlet_lines $d $m
done
echo "all meshes done $(date)"

# ---- grid study: the monitor's line and feature plots, one curve / point per mesh (user, 2026-11-04) ----
G=$R/grid_compare; cd $G
ARGS=$(for m in $MESHES; do echo "${m,}=$(run $m)/2_jet_staged/t3_cap3"; done)
python $MON/plot_grid_lines.py sst_p0j7.0_grid_lines.png "SST p0j 7.0 bar, final (t3) solutions" $ARGS
python $MON/plot_grid_features.py sst_p0j7.0_grid_features.png "SST p0j 7.0 bar, final (t3) solutions" $ARGS
ls -l $G/sst_p0j7.0_grid_*
echo "ALL_DONE $(date)"
