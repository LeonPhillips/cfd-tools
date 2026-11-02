#!/bin/bash
# k-omega SST chain for srp_M2.5, case p0j 8.015 bar / p0inf 3.247 bar (user, 2026-10-30). Built by setup_sst.sh.
#   coarse: Phase 1 (20000), then t1 cap 0.3 / t2 cap 0.5 / t3 cap 1 / t4 cap 3, 50000 iterations each
#   then medium, fine, fineJet, fineJet2, veryFine: interpolation from the previous mesh, t1 cap 0.5 (5000), t2 cap 1 (25000),
#   t3 cap 3 (20000). Every stage runs its full length; the watchers only log and plot (nothing stops on convergence).
#   After each mesh: its line check (last two stages) and SST vs SA lines. At the end: SST grid lines and the
#   0.025 mm box resample of all SST and SA final solutions (for the difference-contour pages).
# Run on the desktop node: tmux new -d -s srp_sst "bash <this folder>/run_sst_chain.sh"   (log: chain.log, appended)
# A mesh whose last stage already has restart.dat is skipped, so the chain can be restarted after a stop.
TOOLS=$(dirname $(readlink -f $0))
R=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs
C=$R/SA_p0j8.015_p0inf3.247_meshCoarse   # line_profiles.py, plot_lines.py, check_run.py
run() { echo $R/SST_p0j8.015_p0inf3.247_mesh$1; }
declare -A SA=([Coarse]=SA_p0j8.015_p0inf3.247_meshCoarse/2_jet_staged/s6_cap3 [Medium]=SA_p0j8.015_p0inf3.247_meshMedium/2_jet_staged/m3_cap3
  [Fine]=SA_p0j8.015_p0inf3.247_meshFine/2_jet_staged/f3_cap3 [FineJet]=SA_p0j8.015_p0inf3.247_meshFineJet/2_jet_staged/j3_cap3
  [FineJet2]=SA_p0j8.015_p0inf3.247_meshFineJet2/2_jet_staged/j3_cap3 [VeryFine]=SA_p0j8.015_p0inf3.247_meshVeryFine/2_jet_staged/v3_cap3)
exec >> $TOOLS/chain.log 2>&1
echo "chain start $(date) on $(hostname), $(nproc) ranks"
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh

su2() {  # su2 <folder> <cfg>: run SU2_CFD there, log to <cfg>.log, return its exit code
  cd $1 || return 1
  { echo "host $(hostname)  ranks $(nproc)"; date; } > ${2%.cfg}.log
  mpirun -n $(nproc) SU2_CFD -t 1 $2 >> ${2%.cfg}.log 2>&1
  local rc=$?; echo "SU2_CFD exit=$rc  $(date)" >> ${2%.cfg}.log; return $rc
}
stage() {  # stage <run folder> <stage> <previous stage folder or ->
  local d=$1/2_jet_staged/$2 rc
  echo "$(basename $1) $2 start $(date), ITER $(sed -n 's/^ITER= *\([0-9]*\).*/\1/p' $d/$2.cfg)"
  (sleep 30; HOLD=0 PREV=$3 bash $TOOLS/watch_fj.sh $d) &
  su2 $d $2.cfg; rc=$?
  wait
  echo "$(basename $1) $2 end exit=$rc $(date), last iteration $(tail -1 $d/history.csv | cut -d, -f3 | tr -d ' ')"
  [ $rc -eq 0 ] && [ -s $d/restart.dat ] || { echo "STAGE FAILED: $(basename $1) $2. Chain stopped."; exit 1; }
}
lines() {  # lines <run folder> <mesh> <stage before last> <last stage>: line check and SST vs SA, into <run>/analysis
  local A=$1/analysis S=$1/2_jet_staged; mkdir -p $A; cd $A
  ( unset LD_PRELOAD; module load paraview/5.13 >/dev/null 2>&1
    pvbatch --mesa $C/line_profiles.py conv_lines.npz ${3%%_*}_end=$S/$3/flow.vtu ${4%%_*}_end=$S/$4/flow.vtu 2>&1 | grep -vi warning
    pvbatch --mesa $C/line_profiles.py sst_vs_sa_lines.npz SA=$R/${SA[$2]}/flow.vtu SST=$S/$4/flow.vtu 2>&1 | grep -vi warning )
  echo "-- $(basename $1): end of $3 vs end of $4"
  python $C/plot_lines.py conv_lines.png conv_lines.npz "SST mesh$2: end of $3 vs end of $4"
  echo "-- $(basename $1): SA vs SST (differences against SST)"
  python $C/plot_lines.py sst_vs_sa_lines.png sst_vs_sa_lines.npz "mesh$2: SA (final) vs SST $4"
}

# ---- coarse from scratch ----
CO=$(run Coarse)
if [ ! -s $CO/2_jet_staged/t4_cap3/restart.dat ]; then
  if [ ! -s $CO/1_freestream/restart.dat ]; then
    echo "Coarse phase1 start $(date)"
    su2 $CO/1_freestream phase1.cfg; rc=$?
    echo "Coarse phase1 end exit=$rc $(date), last iteration $(tail -1 $CO/1_freestream/history.csv | cut -d, -f3 | tr -d ' ')"
    python $TOOLS/plot_conv_roi.py $CO/1_freestream > /dev/null 2>&1
    [ $rc -eq 0 ] && [ -s $CO/1_freestream/restart.dat ] || { echo "PHASE 1 FAILED. Chain stopped."; exit 1; }
  fi
  prev=-
  for s in t1_cap0.3 t2_cap0.5 t3_cap1 t4_cap3; do
    [ -s $CO/2_jet_staged/$s/restart.dat ] || stage $CO $s $prev
    prev=$CO/2_jet_staged/$s
  done
  lines $CO Coarse t3_cap1 t4_cap3
fi

# ---- the other meshes ----
for m in Medium Fine FineJet FineJet2 VeryFine; do
  d=$(run $m)
  [ -s $d/2_jet_staged/t3_cap3/restart.dat ] && { echo "$m already done: skipped"; continue; }
  if [ ! -s $d/0_interp/restart.dat ]; then
    echo "$m interpolation start $(date)"
    su2 $d/0_interp interp.cfg; rc=$?
    grep -q "does not match the mesh" $d/0_interp/interp.log && echo "  (SU2 interpolated: point counts differ)"
    [ $rc -eq 0 ] && [ -s $d/0_interp/restart.dat ] || { echo "INTERPOLATION FAILED: $m. Chain stopped."; exit 1; }
  fi
  prev=-
  for s in t1_cap0.5 t2_cap1 t3_cap3; do
    [ -s $d/2_jet_staged/$s/restart.dat ] || stage $d $s $prev
    prev=$d/2_jet_staged/$s
  done
  lines $d $m t2_cap1 t3_cap3
done
echo "all meshes done $(date)"

# ---- grid study and the box data for the contour pages ----
G=$R/grid_compare; cd $G
( unset LD_PRELOAD; module load paraview/5.13 >/dev/null 2>&1
  pvbatch --mesa $C/line_profiles.py sst_grid_lines.npz coarse=$CO/2_jet_staged/t4_cap3/flow.vtu \
    $(for m in Medium Fine FineJet FineJet2 VeryFine; do echo "${m,}=$(run $m)/2_jet_staged/t3_cap3/flow.vtu"; done | sed 's/^fine\b/fine/') 2>&1 | grep -vi warning
  H_MM=0.025 pvbatch --mesa $G/resample_box.py box_SST_SA_0.025mm.npz SST_coarse=$CO/2_jet_staged/t4_cap3/flow.vtu \
    $(for m in Medium Fine FineJet FineJet2 VeryFine; do echo "SST_${m,}=$(run $m)/2_jet_staged/t3_cap3/flow.vtu"; done) \
    $(for m in Coarse Medium Fine FineJet FineJet2 VeryFine; do echo "SA_${m,}=$R/${SA[$m]}/flow.vtu"; done) 2>&1 | grep -vi warning )
echo "-- SST grid study (differences against veryFine)"
python $C/plot_lines.py sst_grid_lines.png sst_grid_lines.npz "SST final solutions: coarse t4, medium/fine/fineJet/fineJet2/veryFine t3"
ls -l $G/sst_grid_lines.* $G/box_SST_SA_0.025mm.npz
echo "ALL_DONE $(date)"
