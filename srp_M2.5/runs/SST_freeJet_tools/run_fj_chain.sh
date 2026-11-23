#!/bin/bash
# Run one free-jet mesh's chain on Gadi inside a PBS job (user, 2026-11-23):
#   0_interp (medium and fine only): SU2_CFD reads the previous mesh's final restart, interpolates it itself because
#     the point count differs (nearest point + averaging; needs an even number of ranks), and writes it after
#     1 iteration at CFL 1e-8
#   t1 (CFL cap 0.5) in 1000-iteration chunks, 2_jet_staged/t1_cap0.5_c001, _c002, ... After each chunk,
#     jet_arrival.py checks the axial velocity at x = ARR_X (220 mm), r = 0.01 mm. Once it is above ARR_U (20 m/s),
#     one more chunk runs and t1 ends (user, 2026-11-23). For medium and fine the 0_interp solution is checked first,
#     so t1 is a single chunk if the interpolated jet already reaches the outlet. $RUN/t1_final -> the last chunk.
#     If the jet has not arrived after T1_CAP chunks (120 = 120,000 iterations), the chain stops with exit 3 and is
#     flagged (user, 2026-11-23); nothing after t1 runs.
#   then t2, t3, t4 in <run>/2_jet_staged (sorted) with watch_fj_gadi.sh beside SU2
#   then plot_fj.py: the chain figures in <run>.
# Usage: [NP=102] [STOP_AFTER=t1|<stage name>] [T1_CAP=120] bash run_fj_chain.sh <run folder> "<title>"
# Log: <run>/chain.log. Restartable: a step counts as done only if its log says "SU2_CFD exit=0" and its history
# reaches ITER-1. Nothing here stops on convergence; the t1 arrival check only decides when t1 ends.
T=$(dirname $(readlink -f $0))
RUN=$(readlink -f ${1:?run folder}); TITLE=${2:-$(basename $RUN)}
NP=${NP:-102}; STOP_AFTER=${STOP_AFTER:-}; T1_CAP=${T1_CAP:-120}; ARR_X=${ARR_X:-220}; ARR_U=${ARR_U:-20}
S=$RUN/2_jet_staged
exec >> $RUN/chain.log 2>&1
echo "chain start $(date) on $(hostname), job ${PBS_JOBID:-none}, $NP ranks, cfd-tools $(git -C $T rev-parse --short HEAD)"
source /scratch/os85/lp3540/su2/SU2-v8.5.0/su2_env.sh
type module >/dev/null 2>&1 || source /opt/Modules/v4.3.0/init/bash
su2() {  # su2 <folder> <cfg>
  cd $1 || return 1
  { echo "host $(hostname)  ranks $NP  job ${PBS_JOBID:-none}"; date; } > ${2%.cfg}.log
  mpirun -np $NP SU2_CFD -t 1 $2 >> ${2%.cfg}.log 2>&1
  local rc=$?; echo "SU2_CFD exit=$rc  $(date)" >> ${2%.cfg}.log; return $rc
}
pv() { ( unset LD_PRELOAD PYTHONPATH PYTHONHOME VIRTUAL_ENV; module purge >/dev/null 2>&1
         module load paraview/6.1.1-mesa >/dev/null 2>&1; pvbatch "$@" ) 2>&1 | grep -vi warning; }
last_it() { tail -1 $1/history.csv 2>/dev/null | cut -d, -f3 | tr -d ' '; }
done_step() {
  local it n=$(sed -n 's/^ITER= *\([0-9]*\).*/\1/p' $1/$2)
  grep -q "SU2_CFD exit=0" $1/${2%.cfg}.log 2>/dev/null || return 1
  [ -s $1/restart.dat ] || return 1
  it=$(last_it $1); [ -n "$it" ] && [ "$it" -eq $((n - 1)) ]
}
run_stage() {  # run_stage <stage folder> [nowatch]: nowatch = no live watcher (t1 chunks; one snapshot afterwards)
  local d=$1 s=$(basename $1) t0=$(date +%s) rc
  if done_step $d $s.cfg; then echo "$s already complete: skipped"; return 0; fi
  echo "$s start $(date), ITER $(sed -n 's/^ITER= *\([0-9]*\).*/\1/p' $d/$s.cfg)"
  rm -f $d/flow_*.vtu
  [ -z "${2:-}" ] && { (sleep 30; bash $T/watch_fj_gadi.sh $d "$TITLE") & }
  su2 $d $s.cfg; rc=$?
  [ -z "${2:-}" ] && wait
  echo "$s end exit=$rc $(date), $(( $(date +%s) - t0 )) s, last iteration $(last_it $d)"
  [ -z "${2:-}" ] && echo "  $(tail -1 $d/watch.log)"
  done_step $d $s.cfg || { echo "STAGE FAILED: $s. Chain stopped."; exit 1; }
}
arrival() {  # arrival <flow.vtu> <cache file>: prints the jet_arrival.py line (cached)
  [ -s $2 ] || pv $T/jet_arrival.py $1 $ARR_X $ARR_U | tail -1 > $2
  cat $2
}

# ---- 0_interp --------------------------------------------------------------------------------------------------------
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

# ---- 0_seed (coarse): the seeded jet column that t1 starts from ------------------------------------------------------
if [ -d $RUN/0_seed ] && [ ! -s $RUN/0_seed/restart.dat ]; then
  echo "0_seed start $(date) from $(cat $RUN/0_seed/base.txt)"
  python $T/seed_jet.py $(cat $RUN/0_seed/base.txt) $RUN/0_seed/restart.dat > $RUN/0_seed/seed.log 2>&1 \
    || { cat $RUN/0_seed/seed.log; echo "STEP FAILED: 0_seed. Chain stopped."; rm -f $RUN/0_seed/restart.dat; exit 1; }
  sed 's/^/  /' $RUN/0_seed/seed.log
fi

# ---- t1: chunks until the jet reaches the outlet, then one more -------------------------------------------------------
ka=""   # the chunk after which the jet had arrived (0 = already in the 0_interp solution)
if [ -d $RUN/0_interp ] && [ -s $RUN/0_interp/flow.vtu ]; then
  a=$(arrival $RUN/0_interp/flow.vtu $RUN/0_interp/arrival.txt); echo "  0_interp: $a"
  [[ $a == ARRIVED* ]] && ka=0
fi
k=1
while :; do
  c=$(printf 't1_cap0.5_c%03d' $k); d=$S/$c
  if [ ! -f $d/$c.cfg ]; then
    [ $k -gt 1 ] || { echo "missing $d/$c.cfg: chain stopped"; exit 1; }
    mkdir -p $d; sed "s|PREV_CHUNK|$(printf 't1_cap0.5_c%03d' $((k - 1)))|g" $RUN/t1_chunk.cfg > $d/$c.cfg
  fi
  run_stage $d nowatch
  ln -sfn 2_jet_staged/$c $RUN/t1_final
  if [ -s $d/flow.vtu ] && [ ! -f $d/snaps/flow_000999.npz ]; then   # one monitor snapshot per chunk, in the background
    ln -sf flow.vtu $d/flow_000999.vtu; (sleep 65; ONCE=1 bash $T/watch_fj_gadi.sh $d "$TITLE") &
  fi
  if [ -z "$ka" ]; then
    a=$(arrival $d/flow.vtu $d/arrival.txt); echo "  $c: $a"
    [[ $a == ARRIVED* ]] && ka=$k
  fi
  if [ -n "$ka" ] && [ $k -gt $ka ]; then echo "t1 done $(date): jet arrived after chunk $ka; $k chunks, $((k * 1000)) iterations"; break; fi
  if [ -z "$ka" ] && [ $k -ge $T1_CAP ]; then
    echo "T1 CAP REACHED $(date): no jet at x = $ARR_X mm after $k chunks ($((k * 1000)) iterations). FLAGGED: chain stopped."
    wait; exit 3
  fi
  k=$((k + 1))
done
wait
[ "$STOP_AFTER" = t1 ] && { echo "STOP_AFTER=t1 reached $(date)"; exit 0; }

# ---- t2, t3, t4 ------------------------------------------------------------------------------------------------------
for d in $(ls -d $S/t* | grep -v '/t1_cap0.5_c' | sort); do
  s=$(basename $d)
  run_stage $d
  [ "$s" = "$STOP_AFTER" ] && { echo "STOP_AFTER=$s reached $(date)"; exit 0; }
done
echo "all stages done $(date)"
python $T/plot_fj.py $RUN "$TITLE" && ls -l $RUN/*_chain.png
echo "ALL_DONE $(date)"
