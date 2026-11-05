#!/bin/bash
# Run the staged chain of one run folder on Gadi, inside a PBS job (user, 2026-11-05), with the live monitor:
#   for each stage folder in <run>/2_jet_staged (sorted: t1, t2, t3):
#     probes_md.py (ParaView 6.1) on the starting solution -> CUSTOM_OUTPUTS 1 mm either side of the Mach disk
#     SU2_CFD on NP ranks (default 46 of the 48 cores; the other two run the watcher), watch_conv_gadi.sh beside it
#   then plot_chain.py: convergence_chain.png, conv_lines_chain.png, conv_features_chain.png in <run>.
# Start of the first stage: <run>/0_start/flow.vtu; later stages: the previous stage's flow.vtu.
# Usage: bash run_chain_gadi.sh <run folder> "<chain title>"      Log: <run>/chain.log (appended)
# Restartable: a stage counts as done only if its log says "SU2_CFD exit=0" and its history reaches ITER-1.
# Nothing here stops on convergence.
TOOLS=$(dirname $(readlink -f $0)); MON=$(readlink -f $TOOLS/../monitor)
RUN=$(readlink -f ${1:?run folder}); TITLE=${2:-$(basename $RUN)}
NP=${NP:-46}; PV_MODULE=${PV_MODULE:-paraview/6.1.1-mesa}
exec >> $RUN/chain.log 2>&1
echo "chain start $(date) on $(hostname), job ${PBS_JOBID:-none}, $NP ranks, cfd-tools $(git -C $TOOLS rev-parse --short HEAD)"
source /scratch/os85/lp3540/su2/SU2-v8.5.0/su2_env.sh
type module >/dev/null 2>&1 || source /opt/Modules/v4.3.0/init/bash

su2() {  # su2 <folder> <cfg>: run SU2_CFD there, log to <cfg>.log, return its exit code
  cd $1 || return 1
  { echo "host $(hostname)  ranks $NP  job ${PBS_JOBID:-none}"; date; } > ${2%.cfg}.log
  mpirun -np $NP SU2_CFD -t 1 $2 >> ${2%.cfg}.log 2>&1
  local rc=$?; echo "SU2_CFD exit=$rc  $(date)" >> ${2%.cfg}.log; return $rc
}
last_it() { tail -1 $1/history.csv 2>/dev/null | cut -d, -f3 | tr -d ' '; }
done_step() {  # done_step <folder> <cfg>: true if that SU2 run finished all its iterations
  local it n=$(sed -n 's/^ITER= *\([0-9]*\).*/\1/p' $1/$2)
  grep -q "SU2_CFD exit=0" $1/${2%.cfg}.log 2>/dev/null || return 1
  [ -s $1/restart.dat ] || return 1
  it=$(last_it $1); [ -n "$it" ] && [ "$it" -eq $((n - 1)) ]
}
probes() {  # probes <stage folder> <cfg> <starting flow.vtu>
  local p
  p=$(cd $1; unset LD_PRELOAD PYTHONPATH PYTHONHOME VIRTUAL_ENV; module purge >/dev/null 2>&1; module load $PV_MODULE >/dev/null 2>&1
      pvbatch $MON/probes_md.py $3 2> probes.log | grep -v -i warning | tail -1)
  echo "  probes from $3:"; sed 's/^/    /' $1/probes.log
  [[ "$p" == up0* ]] || { echo "  WARNING: probes_md.py failed; the config keeps its previous probes"; return 1; }
  sed -i "s|^CUSTOM_OUTPUTS=.*|CUSTOM_OUTPUTS= '$p'|" $1/$2
}
massflow() {
  python - $1/history.csv <<'PY'
import csv, sys
r = list(csv.reader(open(sys.argv[1]))); h = [c.strip().strip('"') for c in r[0]]
print(r[-1][h.index("Avg_Massflow(inlet_jet)")].strip() if "Avg_Massflow(inlet_jet)" in h else "n/a")
PY
}
stage() {  # stage <stage folder> <previous stage folder or -> <starting flow.vtu>
  local d=$1 s=$(basename $1) rc t0=$(date +%s)
  done_step $d $s.cfg && { echo "$s already complete: skipped"; return 0; }
  echo "$s start $(date), ITER $(sed -n 's/^ITER= *\([0-9]*\).*/\1/p' $d/$s.cfg)"
  probes $d $s.cfg $3
  rm -f $d/flow_*.vtu   # leftovers of an interrupted attempt
  (sleep 30; PREV=$2 bash $MON/watch_conv_gadi.sh $d) &
  su2 $d $s.cfg; rc=$?
  wait
  echo "$s end exit=$rc $(date), $(( $(date +%s) - t0 )) s, last iteration $(last_it $d), jet mass flow $(massflow $d) kg/s"
  echo "  $(python $MON/check_run.py $d 2>&1)"
  done_step $d $s.cfg || { echo "STAGE FAILED: $s. Chain stopped."; exit 1; }
}

prev=-; start=$RUN/0_start/flow.vtu
for d in $(ls -d $RUN/2_jet_staged/t* | sort); do
  stage $d $prev $start; prev=$d; start=$d/flow.vtu
done
echo "all stages done $(date)"
python $MON/plot_chain.py $RUN "$TITLE" && ls -l $RUN/*_chain.png
echo "ALL_DONE $(date)"
