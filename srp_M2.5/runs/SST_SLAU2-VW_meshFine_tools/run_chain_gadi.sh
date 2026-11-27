#!/bin/bash
# Run the staged chain of one run folder on Gadi, inside a PBS job (user, 2026-11-05), with the live monitor:
#   for each stage folder in <run>/2_jet_staged (sorted: t1, t2, t3):
#     probes_md.py (ParaView 6.1) on the starting solution -> CUSTOM_OUTPUTS 1 mm either side of the Mach disk
#     SU2_CFD on NP ranks (default 46 of the 48 cores; the other two run the watcher), watch_conv_gadi.sh beside it
#   then plot_chain.py: convergence_chain.png, conv_lines_chain.png, conv_features_chain.png in <run>.
# Start of the first stage: <run>/0_start/flow.vtu; later stages: the previous stage's flow.vtu.
# Usage: bash run_chain_gadi.sh <run folder> "<chain title>"      Log: <run>/chain.log (appended)
# Restartable: a stage counts as done only if its log says "SU2_CFD exit=0" and its history reaches ITER-1.
# Nothing here stops on convergence, except with AUTOSTOP=1 (CT sweep only, see the end of this file).
TOOLS=$(dirname $(readlink -f $0)); MON=$(readlink -f $TOOLS/../monitor)
RUN=$(readlink -f ${1:?run folder}); TITLE=${2:-$(basename $RUN)}
NP=${NP:-46}   # chain.pbs: 46 of 48 (normal); chain_sr.pbs: 102 of 104 (normalsr); PV_MODULE=${PV_MODULE:-paraview/6.1.1-mesa}
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
  local p try
  # 2026-11-06: on a normalsr node "pvbatch: command not found" here although the watcher's identical call worked, so
  # the module system is (re)initialised inside the subshell, module errors go to probes.log, and it is tried twice.
  for try in 1 2; do
    p=$(cd $1; unset LD_PRELOAD PYTHONPATH PYTHONHOME VIRTUAL_ENV
        source /opt/Modules/v4.3.0/init/bash >/dev/null 2>&1
        { module purge; module load $PV_MODULE; } > probes_module.log 2>&1
        command -v pvbatch >> probes_module.log 2>&1
        pvbatch $MON/probes_md.py $3 2> probes.log | grep -v -i warning | tail -1)
    [[ "$p" == up0* ]] && break
    echo "  probes try $try failed: $(tr '\n' ' ' < $1/probes_module.log | cut -c1-300)"; sleep 5
  done
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

# Deferred start (CT sweep waves, 2026-11-27): 0_start/PARENT.txt names another sweep case; copy the restart.dat and
# flow.vtu of its last stage, but only if that case finished (AUTOSTOP_PASS.txt or AUTOSTOP_CAP.txt).
if [ ! -s $RUN/0_start/restart.dat ] && [ -s $RUN/0_start/PARENT.txt ]; then
  par=$(cat $RUN/0_start/PARENT.txt)
  [ -e $par/AUTOSTOP_PASS.txt ] || [ -e $par/AUTOSTOP_CAP.txt ] || { echo "parent $par has not finished (no AUTOSTOP_PASS/CAP): chain stopped"; exit 1; }
  ps=$(ls -d $par/2_jet_staged/t* | sort | tail -1)
  done_step $ps $(basename $ps).cfg && [ -s $ps/flow.vtu ] || { echo "parent's last stage $ps is not complete: chain stopped"; exit 1; }
  cp $ps/restart.dat $ps/flow.vtu $RUN/0_start/ || exit 1
  echo "start: copies of $ps/{restart.dat,flow.vtu} ($(cat $par/AUTOSTOP_*.txt | cut -c1-120)) $(date)" > $RUN/0_start/SOURCE.txt
  echo "deferred start: $(cat $RUN/0_start/SOURCE.txt)"
fi
[ -s $RUN/0_start/restart.dat ] || { echo "no 0_start/restart.dat: chain stopped"; exit 1; }
prev=-; start=$RUN/0_start/flow.vtu; cum=0
for d in $(ls -d $RUN/2_jet_staged/t* | sort); do
  stage $d $prev $start; prev=$d; start=$d/flow.vtu
  cum=$((cum + $(sed -n 's/^ITER= *\([0-9]*\).*/\1/p' $d/$(basename $d).cfg)))
done
echo "all stages done $(date), cumulative $cum iterations"

# AUTOSTOP=1 (CT sweep, user 2026-11-27: "Auto-stop, cap 250k", limits 0.02 mm / 0.01): after the stages above, run
# CHUNK-iteration continuation stages t4_c01, t4_c02, ... (CFL 1 -> cap 3, from the previous restart) until
# monitor/autostop_check.py PASSes on the latest chunk (checked from MIN_CUM cumulative iterations on) or CAP is reached.
# PASS -> <run>/AUTOSTOP_PASS.txt; cap -> <run>/AUTOSTOP_CAP.txt, the chain is flagged and exits 3. Restartable.
rc=0
if [ "${AUTOSTOP:-0}" = 1 ]; then
  CHUNK=${CHUNK:-10000}; MIN_CUM=${MIN_CUM:-50000}; CAP=${CAP:-250000}
  echo "autostop: chunk $CHUNK, min $MIN_CUM, cap $CAP, limits TOL_MM ${TOL_MM:-0.02} TOL_DMACH ${TOL_DMACH:-0.01} WINDOW ${WINDOW:-10000}"
  pp=$(ls -d $RUN/2_jet_staged/t* | sort | tail -2 | head -1)   # the stage before the last one
  while :; do
    if [ $cum -ge $MIN_CUM ]; then
      res=$(python $MON/autostop_check.py $prev $pp); arc=$?
      echo "  autostop check after $(basename $prev) (cumulative $cum): $res"
      if [ $arc -eq 0 ]; then
        echo "PASS at cumulative iteration $cum ($(basename $prev)) $(date): $res" > $RUN/AUTOSTOP_PASS.txt; break
      fi
    fi
    if [ $cum -ge $CAP ]; then
      echo "CAP reached at cumulative iteration $cum without a PASS $(date): $res" > $RUN/AUTOSTOP_CAP.txt
      echo "AUTOSTOP CAP REACHED: flagged, not converged by the auto-stop test"; rc=3; break
    fi
    n=$(ls -d $RUN/2_jet_staged/t4_c* 2>/dev/null | wc -l); nm=$(printf "t4_c%02d" $((n + 1)))
    bash $TOOLS/add_stage.sh $RUN $(basename $prev) $nm $CHUNK 1.0 3.0 || { echo "add_stage failed"; exit 1; }
    stage $RUN/2_jet_staged/$nm $prev $start
    pp=$prev; prev=$RUN/2_jet_staged/$nm; start=$prev/flow.vtu; cum=$((cum + CHUNK))
  done
fi
python $MON/plot_chain.py $RUN "$TITLE" && ls -l $RUN/*_chain.png
echo "ALL_DONE $(date)"
exit $rc
