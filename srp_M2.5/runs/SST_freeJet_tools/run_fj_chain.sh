#!/bin/bash
# Run one free-jet mesh's chain on Gadi inside a PBS job (user, 2026-11-23):
#   0_interp (medium and fine only): SU2_CFD reads the previous mesh's final restart, interpolates it itself because
#     the point count differs (nearest point + averaging; needs an even number of ranks), and writes it after
#     1 iteration at CFL 1e-8
#   then each stage in <run>/2_jet_staged (sorted: t1, t2, t3, t4) with watch_fj_gadi.sh beside SU2
#   then plot_fj.py: the chain figures in <run>.
# Usage: [NP=102] [STOP_AFTER=<stage name>] bash run_fj_chain.sh <run folder> "<title>"     Log: <run>/chain.log
# Restartable: a step counts as done only if its log says "SU2_CFD exit=0" and its history reaches ITER-1.
# Nothing here stops on convergence.
T=$(dirname $(readlink -f $0))
RUN=$(readlink -f ${1:?run folder}); TITLE=${2:-$(basename $RUN)}
NP=${NP:-102}; STOP_AFTER=${STOP_AFTER:-}
exec >> $RUN/chain.log 2>&1
echo "chain start $(date) on $(hostname), job ${PBS_JOBID:-none}, $NP ranks, cfd-tools $(git -C $T rev-parse --short HEAD)"
source /scratch/os85/lp3540/su2/SU2-v8.5.0/su2_env.sh
su2() {  # su2 <folder> <cfg>
  cd $1 || return 1
  { echo "host $(hostname)  ranks $NP  job ${PBS_JOBID:-none}"; date; } > ${2%.cfg}.log
  mpirun -np $NP SU2_CFD -t 1 $2 >> ${2%.cfg}.log 2>&1
  local rc=$?; echo "SU2_CFD exit=$rc  $(date)" >> ${2%.cfg}.log; return $rc
}
last_it() { tail -1 $1/history.csv 2>/dev/null | cut -d, -f3 | tr -d ' '; }
done_step() {
  local it n=$(sed -n 's/^ITER= *\([0-9]*\).*/\1/p' $1/$2)
  grep -q "SU2_CFD exit=0" $1/${2%.cfg}.log 2>/dev/null || return 1
  [ -s $1/restart.dat ] || return 1
  it=$(last_it $1); [ -n "$it" ] && [ "$it" -eq $((n - 1)) ]
}
if [ -d $RUN/0_interp ]; then
  if done_step $RUN/0_interp interp.cfg; then echo "0_interp already complete: skipped"
  else
    src=$(sed -n 's/^SOLUTION_FILENAME= *//p' $RUN/0_interp/interp.cfg)
    [ -s $src.dat ] || { echo "missing $src.dat: chain stopped"; exit 1; }
    echo "0_interp start $(date) from $src.dat"
    su2 $RUN/0_interp interp.cfg; rc=$?
    grep -i -m3 "interpolat" $RUN/0_interp/interp.log
    echo "0_interp end exit=$rc $(date)"
    done_step $RUN/0_interp interp.cfg || { echo "STEP FAILED: 0_interp. Chain stopped."; exit 1; }
  fi
fi
for d in $(ls -d $RUN/2_jet_staged/t* | sort); do
  s=$(basename $d); t0=$(date +%s)
  if done_step $d $s.cfg; then echo "$s already complete: skipped"
  else
    echo "$s start $(date), ITER $(sed -n 's/^ITER= *\([0-9]*\).*/\1/p' $d/$s.cfg)"
    rm -f $d/flow_*.vtu
    (sleep 30; bash $T/watch_fj_gadi.sh $d "$TITLE") &
    su2 $d $s.cfg; rc=$?
    wait
    echo "$s end exit=$rc $(date), $(( $(date +%s) - t0 )) s, last iteration $(last_it $d)"
    echo "  $(tail -1 $d/watch.log)"
    done_step $d $s.cfg || { echo "STAGE FAILED: $s. Chain stopped."; exit 1; }
  fi
  [ "$s" = "$STOP_AFTER" ] && { echo "STOP_AFTER=$s reached $(date)"; exit 0; }
done
echo "all stages done $(date)"
python $T/plot_fj.py $RUN "$TITLE" && ls -l $RUN/*_chain.png
echo "ALL_DONE $(date)"
