#!/bin/bash
# Results of the meshMedium Phase 2 run (m1, m2) against coarse s6. Run on the desktop node:
#   bash <this folder>/analyse_medium.sh   (writes into analysis/)
C=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs/SA_p0j8.015_p0inf3.247_meshCoarse
M=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs/SA_p0j8.015_p0inf3.247_meshMedium
mkdir -p $M/analysis; cd $M/analysis || exit 1
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh >/dev/null 2>&1
cd $C
for s in m1_cap0.5 m2_cap1 m3_cap3; do
  echo "== check_run $s"; python check_run.py ../SA_p0j8.015_p0inf3.247_meshMedium/2_jet_staged/$s
  python plot_history.py ../SA_p0j8.015_p0inf3.247_meshMedium/2_jet_staged/$s 10000 | tail -8
done
echo "== check_run coarse s6"; python check_run.py 2_jet_staged/s6_cap3
cd $M/analysis
python3 - $M/2_jet_staged/m1_cap0.5/history.csv $M/2_jet_staged/m2_cap1/history.csv $M/2_jet_staged/m3_cap3/history.csv <<'PY'
import sys
import numpy as np
for p in sys.argv[1:]:
    with open(p) as fh:
        h = [c.strip().strip('"') for c in fh.readline().split(",")]
    d = np.genfromtxt(p, delimiter=",", skip_header=1)
    col = {k: d[:, i] for i, k in enumerate(h)}
    cfl = col.get("Avg CFL", col.get("Avg_CFL", None))
    drops = int(np.sum(np.diff(cfl) < -0.25 * cfl[:-1])) if cfl is not None else -1
    print(p.split("/")[-2], f"iterations {len(d)}  CFL drops {drops}  CD every 5000:",
          [round(float(col['CD'][i]), 5) for i in range(0, len(d), 5000)] + [round(float(col['CD'][-1]), 5)])
PY
unset LD_PRELOAD; module load paraview/5.13 >/dev/null 2>&1
F=$M/2_jet_staged/m3_cap3/flow.vtu
echo "== inspect m3"; pvbatch --mesa $C/inspect_field.py $F medium_m3.npz 2>&1 | grep -vi warning | grep -E "bow shock|outlet|inlet \(|backflow"
echo "== axis probe"; pvbatch --mesa $M/0_interp/axis_probe.py $C/2_jet_staged/s6_cap3/flow.vtu $F 2>&1 | grep -vi warning
echo "== bad interpolation points still there? (|dMach| > 0.5 vs coarse s6)"
pvbatch --mesa $M/0_interp/bad_points.py $C/2_jet_staged/s6_cap3/flow.vtu $F 2>&1 | grep -vi warning | head -14
echo "== field difference m3 vs coarse s6"; pvbatch --mesa $M/0_interp/diff_field.py $C/2_jet_staged/s6_cap3/flow.vtu $F 2>&1 | grep -vi warning
module unload paraview/5.13 >/dev/null 2>&1
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh >/dev/null 2>&1
python $C/plot_field.py coarse-s6_vs_medium-m3.png $M/0_interp/coarse_s6.npz medium_m3.npz
for s in m1_cap0.5 m2_cap1 m3_cap3; do cp $M/2_jet_staged/$s/convergence.png ${s%%_*}_convergence.png; done
ls -l
