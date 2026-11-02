#!/bin/bash
# k-omega SST + Sarkar compressibility correction chain for srp_M2.5, case p0j 8.015 bar / p0inf 3.247 bar (user, 2026-10-31).
# Same stages as SST_p0j8.015_p0inf3.247_tools/run_sst_chain.sh; configs built by setup_sstsarkar.sh.
#   coarse: Phase 1 (20000, CFL 0.1 -> cap 0.5, floor 0.1), then t1 cap 0.3 / t2 cap 0.5 / t3 cap 1 / t4 cap 3, 50000 each
#   then medium, fine, fineJet, fineJet2, veryFine (veryFine from fine): interpolation, t1 cap 0.5 (5000), t2 cap 1 (25000),
#   t3 cap 3 (20000). Every stage runs its full length; the watchers only log and plot (nothing stops on convergence).
#   After each mesh: its line check (t2 vs t3; coarse t3 vs t4) and SA vs SST vs SSTsarkar lines. At the end: SSTsarkar grid
#   lines and the 0.025 mm box resample of the SSTsarkar final solutions.
# Run on the desktop node: tmux new -d -s srp_sstc "bash <this folder>/run_sstsarkar_chain.sh"   (log: chain.log, appended)
# Optional SU2_CORES=<comma list>: run SU2 on exactly these cores, one rank each (Intel MPI pin list). On m3f100
# (2026-11-01) 24 ranks on all 24 allocated cores ran at 2.3 s/iteration; 22 ranks without cores 61 and 84 (busy with
# outside load) ran at 0.042 s/iteration (speedtest/results*.txt). Keep the rank count even (interpolation).
# Restartable: a step counts as done only if its log says "SU2_CFD exit=0" and its history reaches ITER-1
# (restart.dat alone is not enough: SU2 writes it every 1000 iterations).
TOOLS=$(dirname $(readlink -f $0))
R=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs
C=$R/SA_p0j8.015_p0inf3.247_meshCoarse   # line_profiles.py, plot_lines.py, check_run.py
run() { echo $R/SSTsarkar_p0j8.015_p0inf3.247_mesh$1; }
sst() { echo $R/SST_p0j8.015_p0inf3.247_mesh$1; }
declare -A SA=([Coarse]=SA_p0j8.015_p0inf3.247_meshCoarse/2_jet_staged/s6_cap3 [Medium]=SA_p0j8.015_p0inf3.247_meshMedium/2_jet_staged/m3_cap3
  [Fine]=SA_p0j8.015_p0inf3.247_meshFine/2_jet_staged/f3_cap3 [FineJet]=SA_p0j8.015_p0inf3.247_meshFineJet/2_jet_staged/j3_cap3
  [FineJet2]=SA_p0j8.015_p0inf3.247_meshFineJet2/2_jet_staged/j3_cap3 [VeryFine]=SA_p0j8.015_p0inf3.247_meshVeryFine/2_jet_staged/v3_cap3)
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
lines() {  # lines <run folder> <mesh> <stage before last> <last stage> <SST last stage>: into <run>/analysis
  local A=$1/analysis S=$1/2_jet_staged; mkdir -p $A; cd $A
  ( unset LD_PRELOAD; module load paraview/5.13 >/dev/null 2>&1
    pvbatch --mesa $C/line_profiles.py conv_lines.npz ${3%%_*}_end=$S/$3/flow.vtu ${4%%_*}_end=$S/$4/flow.vtu 2>&1 | grep -vi warning
    pvbatch --mesa $C/line_profiles.py models_lines.npz SA=$R/${SA[$2]}/flow.vtu SST=$(sst $2)/2_jet_staged/$5/flow.vtu \
      SSTsarkar=$S/$4/flow.vtu 2>&1 | grep -vi warning )
  echo "-- $(basename $1): end of $3 vs end of $4"
  python $C/plot_lines.py conv_lines.png conv_lines.npz "SSTsarkar mesh$2: end of $3 vs end of $4"
  echo "-- $(basename $1): SA and SST vs SSTsarkar (differences against SSTsarkar)"
  python $C/plot_lines.py models_lines.png models_lines.npz "mesh$2 final solutions: SA, SST, SST + Sarkar"
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
[ -s $CO/analysis/models_lines.png ] || lines $CO Coarse t3_cap1 t4_cap3 t4_cap3

# ---- the other meshes ----
for m in Medium Fine FineJet FineJet2 VeryFine; do
  d=$(run $m)
  if ! done_step $d/0_interp interp.cfg; then
    echo "$m interpolation start $(date)"
    su2 $d/0_interp interp.cfg; rc=$?
    grep -q "does not match the mesh" $d/0_interp/interp.log && echo "  (SU2 interpolated: point counts differ)"
    done_step $d/0_interp interp.cfg || { echo "INTERPOLATION FAILED: $m. Chain stopped."; exit 1; }
  fi
  prev=-
  for s in t1_cap0.5 t2_cap1 t3_cap3; do stage $d $s $prev; prev=$d/2_jet_staged/$s; done
  [ -s $d/analysis/models_lines.png ] || lines $d $m t2_cap1 t3_cap3 t3_cap3
done
echo "all meshes done $(date)"

# ---- grid study and the box data for the contour pages ----
G=$R/grid_compare; cd $G
( unset LD_PRELOAD; module load paraview/5.13 >/dev/null 2>&1
  pvbatch --mesa $C/line_profiles.py sstsarkar_grid_lines.npz coarse=$CO/2_jet_staged/t4_cap3/flow.vtu \
    $(for m in Medium Fine FineJet FineJet2 VeryFine; do echo "${m,}=$(run $m)/2_jet_staged/t3_cap3/flow.vtu"; done) 2>&1 | grep -vi warning
  H_MM=0.025 pvbatch --mesa $G/resample_box.py box_SSTsarkar_0.025mm.npz SSTsarkar_coarse=$CO/2_jet_staged/t4_cap3/flow.vtu \
    $(for m in Medium Fine FineJet FineJet2 VeryFine; do echo "SSTsarkar_${m,}=$(run $m)/2_jet_staged/t3_cap3/flow.vtu"; done) 2>&1 | grep -vi warning )
echo "-- SSTsarkar grid study (differences against veryFine)"
python $C/plot_lines.py sstsarkar_grid_lines.png sstsarkar_grid_lines.npz "SST + Sarkar final solutions: coarse t4, medium/fine/fineJet/fineJet2/veryFine t3"
ls -l $G/sstsarkar_grid_lines.* $G/box_SSTsarkar_0.025mm.npz
echo "ALL_DONE $(date)"
