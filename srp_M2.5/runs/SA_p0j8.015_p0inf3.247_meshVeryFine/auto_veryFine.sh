#!/bin/bash
# meshVeryFine, unattended (user, 2026-10-28): wait for the mesh job, check it, interpolate fineNozzle n3 onto the new mesh,
# run v1 -> v2 -> v3 (run_veryFine.sh, never stops on convergence), then the line-profile analysis.
# Run on the desktop node: tmux new -d -s srp_auto2 "bash <this folder>/auto_veryFine.sh"  (log: auto.log)
R=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5
C=$R/runs/SA_p0j8.015_p0inf3.247_meshCoarse; MD=$R/runs/SA_p0j8.015_p0inf3.247_meshMedium
F=$R/runs/SA_p0j8.015_p0inf3.247_meshFine; V=$R/runs/SA_p0j8.015_p0inf3.247_meshVeryFine
ML=$R/mesh/mesh_job_veryFine.log
exec > $V/auto.log 2>&1
echo "start $(date)"
until grep -q MESH_JOB_DONE $ML 2>/dev/null; do sleep 30; done
echo "== mesh job done $(date)"; cat $ML
grep -q "mesh exit=0" $ML && [ -s $R/mesh/srp_M2.5_meshVeryFine.su2 ] || { echo "MESH FAILED: chain stopped"; exit 1; }
grep -q "meshcheck exit=0" $ML || { echo "SU2 MESHCHECK FAILED: chain stopped"; exit 1; }

echo "== interpolation $(date)"
bash $V/0_interp/interp.sh
tail -3 $V/0_interp/interp.log
grep -q "SU2_CFD exit=0" $V/0_interp/interp.log && [ -s $V/0_interp/restart.dat ] || { echo "INTERPOLATION FAILED: chain stopped"; exit 1; }

echo "== runs $(date)"
bash $V/run_veryFine.sh
cat $V/run_veryFine.log
grep -q "all stages done" $V/run_veryFine.log || echo "RUN CHAIN DID NOT FINISH: analysing what exists"

echo "== analysis $(date)"
A=$V/analysis
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh >/dev/null 2>&1
cd $C
for s in v1_cap0.5 v2_cap1 v3_cap3; do
  [ -f $V/2_jet_staged/$s/history.csv ] || continue
  echo "-- check_run $s"; python check_run.py ../SA_p0j8.015_p0inf3.247_meshVeryFine/2_jet_staged/$s
  python plot_history.py ../SA_p0j8.015_p0inf3.247_meshVeryFine/2_jet_staged/$s 10000 | tail -6
  cp $V/2_jet_staged/$s/convergence.png $A/${s%%_*}_convergence.png
done
cd $A
V2=$V/2_jet_staged/v2_cap1/flow.vtu; V3=$V/2_jet_staged/v3_cap3/flow.vtu
unset LD_PRELOAD; module load paraview/5.13 >/dev/null 2>&1
pvbatch --mesa $C/line_profiles.py grid_lines.npz coarse=$C/2_jet_staged/s6_cap3/flow.vtu medium=$MD/2_jet_staged/m3_cap3/flow.vtu \
  fine=$F/2_jet_staged/f3_cap3/flow.vtu veryFine=$V3 2>&1 | grep -vi warning
pvbatch --mesa $C/line_profiles.py conv_lines_v2_v3.npz v2_end=$V2 v3_end=$V3 2>&1 | grep -vi warning
module unload paraview/5.13 >/dev/null 2>&1
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh >/dev/null 2>&1
echo "-- grid: coarse s6, medium m3, fine f3, veryFine v3 (differences against veryFine)"
python $C/plot_lines.py grid_lines.png grid_lines.npz "Final solutions: coarse s6, medium m3, fine f3, very fine v3"
echo "-- convergence: v2 end vs v3 end"
python $C/plot_lines.py conv_lines_v2_v3.png conv_lines_v2_v3.npz "meshVeryFine: end of v2 vs end of v3"
echo "-- jet mass flow at the end of each stage (fine f3 for reference)"
for h in $F/2_jet_staged/f3_cap3/history.csv $V/2_jet_staged/v*/history.csv; do
  python3 - $h <<'PY'
import sys
import numpy as np
p = sys.argv[1]
with open(p) as fh:
    h = [c.strip().strip('"') for c in fh.readline().split(",")]
d = np.genfromtxt(p, delimiter=",", skip_header=1)[-1]
v = dict(zip(h, d))
mf = [k for k in h if k.startswith("Avg_Massflow")]
print(p.split("/")[-2], " ".join(f"{k}={v[k]:.6g}" for k in mf), f"CD={v['CD']:.5f}")
PY
done
ls -l $A
echo "AUTO_DONE $(date)"
