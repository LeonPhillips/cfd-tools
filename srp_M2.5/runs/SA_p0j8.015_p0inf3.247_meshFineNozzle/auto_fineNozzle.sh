#!/bin/bash
# meshFineNozzle, unattended (user, 2026-10-28, so the laptop can be closed): wait for the mesh job, check it, interpolate
# fine f3 onto the new mesh, run n1 -> n2 -> n3 (run_fineNozzle.sh, never stops on convergence), then the line-profile
# analysis. Run on the desktop node: tmux new -d -s srp_auto "bash <this folder>/auto_fineNozzle.sh"  (log: auto.log)
R=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5
C=$R/runs/SA_p0j8.015_p0inf3.247_meshCoarse; F=$R/runs/SA_p0j8.015_p0inf3.247_meshFine
N=$R/runs/SA_p0j8.015_p0inf3.247_meshFineNozzle; ML=$R/mesh/mesh_job_fineNozzle.log
exec > $N/auto.log 2>&1
echo "start $(date)"
until grep -q MESH_JOB_DONE $ML; do sleep 30; done
echo "== mesh job done $(date)"; cat $ML
grep -q "mesh exit=0" $ML && [ -s $R/mesh/srp_M2.5_meshFineNozzle.su2 ] || { echo "MESH FAILED: chain stopped"; exit 1; }
grep -q "meshcheck exit=0" $ML || { echo "SU2 MESHCHECK FAILED: chain stopped"; exit 1; }

echo "== interpolation $(date)"
bash $N/0_interp/interp.sh
tail -3 $N/0_interp/interp.log
grep -q "SU2_CFD exit=0" $N/0_interp/interp.log && [ -s $N/0_interp/restart.dat ] || { echo "INTERPOLATION FAILED: chain stopped"; exit 1; }

echo "== runs $(date)"
bash $N/run_fineNozzle.sh
cat $N/run_fineNozzle.log
grep -q "all stages done" $N/run_fineNozzle.log || echo "RUN CHAIN DID NOT FINISH: analysing what exists"

echo "== analysis $(date)"
A=$N/analysis; cd $A || exit 1
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh >/dev/null 2>&1
cd $C
for s in n1_cap0.5 n2_cap1 n3_cap3; do
  [ -f $N/2_jet_staged/$s/history.csv ] || continue
  echo "-- check_run $s"; python check_run.py ../SA_p0j8.015_p0inf3.247_meshFineNozzle/2_jet_staged/$s
  python plot_history.py ../SA_p0j8.015_p0inf3.247_meshFineNozzle/2_jet_staged/$s 10000 | tail -6
  cp $N/2_jet_staged/$s/convergence.png $A/${s%%_*}_convergence.png
done
cd $A
F3=$F/2_jet_staged/f3_cap3/flow.vtu; N2=$N/2_jet_staged/n2_cap1/flow.vtu; N3=$N/2_jet_staged/n3_cap3/flow.vtu
NOZ="-8.2,0.01,2,0.01;1.9,0,1.9,3"   # axis through the nozzle; radial line across the jet exit
unset LD_PRELOAD; module load paraview/5.13 >/dev/null 2>&1
pvbatch --mesa $C/line_profiles.py grid_lines.npz fine=$F3 fineNozzle=$N3 2>&1 | grep -vi warning
pvbatch --mesa $C/line_profiles.py conv_lines_n2_n3.npz n2_end=$N2 n3_end=$N3 2>&1 | grep -vi warning
LINES="$NOZ" pvbatch --mesa $C/line_profiles.py nozzle_lines.npz fine=$F3 fineNozzle=$N3 2>&1 | grep -vi warning
LINES="$NOZ" pvbatch --mesa $C/line_profiles.py conv_nozzle_n2_n3.npz n2_end=$N2 n3_end=$N3 2>&1 | grep -vi warning
module unload paraview/5.13 >/dev/null 2>&1
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh >/dev/null 2>&1
echo "-- grid: fine f3 vs fineNozzle n3, standard lines"; python $C/plot_lines.py grid_lines.png grid_lines.npz "Standard lines: fine f3 vs fineNozzle n3"
echo "-- grid: fine f3 vs fineNozzle n3, nozzle lines"; python $C/plot_lines.py nozzle_lines.png nozzle_lines.npz "Nozzle lines: fine f3 vs fineNozzle n3"
echo "-- convergence: n2 end vs n3 end, standard lines"; python $C/plot_lines.py conv_lines_n2_n3.png conv_lines_n2_n3.npz "meshFineNozzle: end of n2 vs end of n3"
echo "-- convergence: n2 end vs n3 end, nozzle lines"; python $C/plot_lines.py conv_nozzle_n2_n3.png conv_nozzle_n2_n3.npz "meshFineNozzle nozzle lines: end of n2 vs end of n3"
echo "-- jet mass flow and imbalance at the end of each stage (fine f3 for reference)"
for h in $F/2_jet_staged/f3_cap3/history.csv $N/2_jet_staged/n*/history.csv; do
  python3 - $h <<'PY'
import sys
import numpy as np
p = sys.argv[1]
with open(p) as fh:
    h = [c.strip().strip('"') for c in fh.readline().split(",")]
d = np.genfromtxt(p, delimiter=",", skip_header=1)[-1]
v = dict(zip(h, d))
mf = [k for k in h if k.startswith("Avg_Massflow")]
print(p.split("/")[-3][-14:], p.split("/")[-2], " ".join(f"{k}={v[k]:.6g}" for k in mf), f"CD={v['CD']:.5f}")
PY
done
ls -l $A
echo "AUTO_DONE $(date)"
