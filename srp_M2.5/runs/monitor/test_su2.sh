#!/bin/bash
# Live test of the monitor with SU2 (user, 2026-11-04): SST p0j 7.4945 bar coarse, restart from its t3 solution, the t3
# config with ITER=${ITER:-3000} and the Mach-disk probes from probes_md.py; watch_conv.sh runs alongside.
# Nothing physical changes from t3 (same CFL ramp to cap 3, same inlet); it only tests the monitor.
# Run inside the desktop job: srun --jobid=<job> --overlap -n1 -c24 bash test_su2.sh   (NP ranks, default 22)
# Output: <coarse run>/test_monitor_su2/2_jet_staged/m_su2/ (remove when done).
MON=$(dirname $(readlink -f $0))
R=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs/SST_p0j7.4945_p0inf3.247_meshCoarse
S=$R/2_jet_staged/t3_cap3; T=$R/test_monitor_su2; D=$T/2_jet_staged/m_su2; NP=${NP:-22}
[ -e $T ] && { echo "$T exists: remove it first"; exit 1; }
mkdir -p $D; cd $D
echo "host $(hostname), cpus allowed $(grep Cpus_allowed_list /proc/self/status | cut -f2), $NP ranks, $(date)"
( unset LD_PRELOAD; module load paraview/5.13 >/dev/null 2>&1
  pvbatch --mesa $MON/probes_md.py $S/flow.vtu 2> probes.log | grep -v -i warning | tail -1 > probes.txt )
cat probes.log
[ -s probes.txt ] || { echo "probes_md.py failed"; exit 1; }
{ echo "% Monitor test (2026-11-04): copy of ../../../2_jet_staged/t3_cap3/t3_cap3.cfg, restart from its solution,"
  echo "% ITER ${ITER:-3000}, probes 1 mm either side of the Mach disk (probes_md.py). Made by cfd-tools monitor/test_su2.sh."
  sed -e "s|^SOLUTION_FILENAME=.*|SOLUTION_FILENAME= $S/restart|" \
      -e "s|^ITER=.*|ITER= ${ITER:-3000}|" \
      -e "s|^CUSTOM_OUTPUTS=.*|CUSTOM_OUTPUTS= '$(cat probes.txt)'|" $S/t3_cap3.cfg
} > m_su2.cfg
grep -E "^(SOLUTION_FILENAME|ITER|OUTPUT_WRT_FREQ|WRT_VOLUME_OVERWRITE|CONV_RESIDUAL_MINVAL)=" m_su2.cfg
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh >/dev/null
export I_MPI_HYDRA_BOOTSTRAP=fork   # inside an srun step: start the ranks locally, not as a new Slurm step
(sleep 30; PREV=$S bash $MON/watch_conv.sh $D) &
{ echo "host $(hostname)  ranks $NP"; date; } > m_su2.log
mpirun -n $NP SU2_CFD -t 1 m_su2.cfg >> m_su2.log 2>&1
rc=$?; echo "SU2_CFD exit=$rc  $(date)" | tee -a m_su2.log
wait
echo "== watch.log"; cat watch.log
for f in box_change line_change features; do echo "== $f.csv"; cat $f.csv; done
echo "== probe columns in history"; head -1 history.csv | tr ',' '\n' | grep -n -E '"?(up|dn)[0-9]' | tr '\n' ' '; echo
ls -l $D $D/snaps $T/*_latest.png; du -sh $D
