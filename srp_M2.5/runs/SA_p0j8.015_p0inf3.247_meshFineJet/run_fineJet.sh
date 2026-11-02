#!/bin/bash
# meshFineJet Phase 2. Revised by the user on 2026-10-29: j1 (cap 0.5, 5000 iterations) is done and kept. The early exit
# on a held CFL is removed, as are the cap 5 / cap 10 stages. j2 (cap 1) runs its full 25000 iterations from the end of
# j1, then j3 (cap 3) its full 20000. Nothing stops on convergence. Then the line-profile analysis.
# Run on the desktop node: tmux new -d -s srp_finejet "bash <this folder>/run_fineJet.sh"   (log: run_fineJet.log, appended)
FJ=$(dirname $(readlink -f $0))
B=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs
C=$B/SA_p0j8.015_p0inf3.247_meshCoarse; F=$B/SA_p0j8.015_p0inf3.247_meshFine; V=$B/SA_p0j8.015_p0inf3.247_meshVeryFine
exec >> $FJ/run_fineJet.log 2>&1
echo "restart $(date): j2 and j3 at full length, no early exit (user, 2026-10-29)"
if [ ! -s $FJ/0_interp/restart.dat ]; then
  bash $FJ/0_interp/interp.sh; tail -2 $FJ/0_interp/interp.log
  grep -q "SU2_CFD exit=0" $FJ/0_interp/interp.log && [ -s $FJ/0_interp/restart.dat ] || { echo "INTERPOLATION FAILED"; exit 1; }
fi
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh

# run_stage <stage> <previous stage folder>
run_stage() {
  local s=$1 d=$FJ/2_jet_staged/$1 rc
  cd $d || exit 1
  echo "$s start $(date), ITER $(sed -n 's/^ITER= *\([0-9]*\).*/\1/p' $s.cfg)"
  { echo "host $(hostname)  ranks $(nproc)"; date; } > $s.log
  (sleep 30; HOLD=0 PREV=$2 bash $FJ/watch_fj.sh $d) &
  mpirun -n $(nproc) SU2_CFD -t 1 $s.cfg >> $s.log 2>&1
  rc=$?
  echo "SU2_CFD exit=$rc  $(date)" >> $s.log
  wait
  echo "$s end exit=$rc $(date), last iteration $(tail -1 history.csv | cut -d, -f3 | tr -d ' ')"
  [ $rc -eq 0 ] && [ -f restart.dat ] || { echo "stage $s failed: chain stopped"; exit 1; }
}
S=$FJ/2_jet_staged
run_stage j2_cap1 $S/j1_cap0.5
run_stage j3_cap3 $S/j2_cap1
ran=(j1_cap0.5 j2_cap1 j3_cap3)
last=${ran[-1]}; before=${ran[-2]}
echo "all stages done $(date): last stage $last"

echo "== analysis $(date)"
A=$FJ/analysis; mkdir -p $A; cd $A
unset LD_PRELOAD; module load paraview/5.13 >/dev/null 2>&1
pvbatch --mesa $C/line_profiles.py conv_lines.npz ${before%%_*}_end=$S/$before/flow.vtu ${last%%_*}_end=$S/$last/flow.vtu 2>&1 | grep -vi warning
pvbatch --mesa $C/line_profiles.py grid_lines.npz fine=$F/2_jet_staged/f3_cap3/flow.vtu fineJet=$S/$last/flow.vtu \
  veryFine=$V/2_jet_staged/v3_cap3/flow.vtu 2>&1 | grep -vi warning
module unload paraview/5.13 >/dev/null 2>&1
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh
echo "-- convergence: end of $before vs end of $last"
python $C/plot_lines.py conv_lines.png conv_lines.npz "meshFineJet: end of $before vs end of $last"
echo "-- grid: fine f3, fineJet $last, veryFine v3 (differences against veryFine)"
python $C/plot_lines.py grid_lines.png grid_lines.npz "Final solutions: fine f3, fineJet $last, very fine v3"
echo "ALL_DONE $(date)"
