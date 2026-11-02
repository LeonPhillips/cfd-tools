#!/bin/bash
# Short checks before the corrected-inlet runs (user, 2026-11-02). Needs setup_p0j7.4945.sh and setup_rhoj8.681.sh done.
#   A: total conditions p0j 7.4945 bar, meshMedium t1 for 100 iterations from the SST 8.015 bar t3 solution
#   B: MASS_FLOW, meshCoarse t1 (jet switched on) for 300 iterations from the SST 8.015 bar Phase 1 solution
#      (Phase 1 has the nozzle plugged, so it is the same field the MASS_FLOW Phase 1 will give)
#   C: MASS_FLOW, SST 8.015 bar coarse t4 interpolated onto meshFineJet, then t1 for 100 iterations
# Each test goes in a test/ folder of its run; prints exit code, last row and jet-inlet averages.
# Run on the desktop node: SU2_CORES=<list> bash test_inlets.sh   (log: test_inlets.log next to this script)
R=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs
OLD=SST_p0j8.015_p0inf3.247; TC=SST_p0j7.4945_p0inf3.247; MF=SST_rhoj8.681_Uj99.11_p0inf3.247
exec > $(dirname $(readlink -f $0))/test_inlets.log 2>&1
if [ -n "$SU2_CORES" ]; then NP=$(echo $SU2_CORES | tr "," " " | wc -w); export I_MPI_PIN_PROCESSOR_LIST=$SU2_CORES; else NP=$(nproc); fi
echo "tests start $(date) on $(hostname), $NP ranks"
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh

# mk <test folder> <source cfg> <iterations> <solution file or ->: copy the cfg with a short ITER and frequent output
mk() {
  mkdir -p $1; local g=$1/$(basename $2)
  sed -e "s|^ITER=.*|ITER= $3|" -e "s|^OUTPUT_WRT_FREQ=.*|OUTPUT_WRT_FREQ= 100000, 100000, 100000|" \
      -e "s|^WRT_VOLUME_OVERWRITE=.*|WRT_VOLUME_OVERWRITE= YES|" $2 > $g
  [ "$4" != - ] && sed -i "s|^SOLUTION_FILENAME=.*|SOLUTION_FILENAME= $4|" $g
  echo $g
}
run() {  # run <cfg>
  cd $(dirname $1); local c=$(basename $1)
  mpirun -n $NP SU2_CFD -t 1 $c > ${c%.cfg}.log 2>&1; local rc=$?
  echo "== $(dirname $1)/$c: exit=$rc"
  grep -i "error\|nan" ${c%.cfg}.log | grep -vi "linear solver\|LINEAR_SOLVER_ERROR" | head -5
  [ -f history.csv ] && awk -F, 'NR == 1 { for (i = 1; i <= NF; i++) { h = $i; gsub(/[" ]/, "", h); col[h] = i }; next }
    { last = $0; n++; split($0, v, ","); if (n == 1) for (k in col) first[k] = v[col[k]] }
    END { split(last, v, ",")
      printf "  iterations %d  rms[Rho] %s -> %s  AvgCFL %s  CD %s\n", n, first["rms[Rho]"], v[col["rms[Rho]"]], v[col["AvgCFL"]], v[col["CD"]]
      printf "  inlet_jet: p %s Pa  rho %s  T %s K  Vn %s m/s  M %s  p0 %s Pa  T0 %s K  mdot %s kg/s\n",
        v[col["Avg_Press(inlet_jet)"]], v[col["Avg_Density(inlet_jet)"]], v[col["Avg_Temp(inlet_jet)"]],
        v[col["Avg_NormalVel(inlet_jet)"]], v[col["Avg_Mach(inlet_jet)"]], v[col["Avg_TotalPress(inlet_jet)"]],
        v[col["Avg_TotalTemp(inlet_jet)"]], v[col["Avg_Massflow(inlet_jet)"]]
      printf "  outlet mdot %s, freestream inlet mdot %s\n", v[col["Avg_Massflow(outlet)"]], v[col["Avg_Massflow(inlet_freestream)"]] }' history.csv
  return $rc
}

run $(mk $R/${TC}_meshMedium/test $R/${TC}_meshMedium/2_jet_staged/t1_cap0.5/t1_cap0.5.cfg 100 -)
run $(mk $R/${MF}_meshCoarse/test $R/${MF}_meshCoarse/2_jet_staged/t1_cap0.3/t1_cap0.3.cfg 300 $R/${OLD}_meshCoarse/1_freestream/restart)
i=$(mk $R/${MF}_meshFineJet/test/0_interp $R/${MF}_meshFineJet/0_interp/interp.cfg 1 $R/${OLD}_meshCoarse/2_jet_staged/t4_cap3/restart)
sed -i "s|^OUTPUT_WRT_FREQ=.*|OUTPUT_WRT_FREQ= 1|" $i
run $i && run $(mk $R/${MF}_meshFineJet/test/t1 $R/${MF}_meshFineJet/2_jet_staged/t1_cap0.5/t1_cap0.5.cfg 100 ../0_interp/restart)
echo "tests done $(date)"
