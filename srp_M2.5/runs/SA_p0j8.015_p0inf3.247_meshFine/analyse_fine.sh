#!/bin/bash
# Results of the meshFine Phase 2 run (f1-f3) against meshMedium m3 and meshCoarse s6. Run on the desktop node:
#   bash <this folder>/analyse_fine.sh < /dev/null   (writes into analysis/)
R=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs
C=$R/SA_p0j8.015_p0inf3.247_meshCoarse; MD=$R/SA_p0j8.015_p0inf3.247_meshMedium; F=$R/SA_p0j8.015_p0inf3.247_meshFine
mkdir -p $F/analysis; cd $F/analysis || exit 1
cp $MD/0_interp/axis_probe.py $MD/0_interp/diff_field.py $MD/0_interp/bad_points.py .
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh >/dev/null 2>&1
cd $C
for s in f1_cap0.5 f2_cap1 f3_cap3; do
  echo "== check_run $s"; python check_run.py ../SA_p0j8.015_p0inf3.247_meshFine/2_jet_staged/$s
  python plot_history.py ../SA_p0j8.015_p0inf3.247_meshFine/2_jet_staged/$s 10000 | tail -6
done
cd $F/analysis
echo "== CD trend and CFL drops per stage"
python3 - $F/2_jet_staged/f1_cap0.5/history.csv $F/2_jet_staged/f2_cap1/history.csv $F/2_jet_staged/f3_cap3/history.csv <<'PY'
import sys
import numpy as np
for p in sys.argv[1:]:
    with open(p) as fh:
        h = [c.strip().strip('"') for c in fh.readline().split(",")]
    d = np.genfromtxt(p, delimiter=",", skip_header=1)
    col = {k: d[:, i] for i, k in enumerate(h)}
    cfl = col["Avg CFL"] if "Avg CFL" in col else col.get("Avg_CFL")
    drops = int(np.sum(np.diff(cfl) < -0.25 * cfl[:-1])) if cfl is not None else -1
    cd = col["CD"]
    print(p.split("/")[-2], f"iterations {len(d)}  CFL drops {drops}  CD every 5000:",
          [round(float(cd[i]), 5) for i in range(0, len(d), 5000)] + [round(float(cd[-1]), 5)],
          f" CD change over last 5000: {cd[-1]-cd[max(0,len(cd)-5001)]:+.5f}")
PY
cat > axis_wiggle.py <<'PY'
"""pvbatch --mesa axis_wiggle.py <flow.vtu> ...: smoothness of Mach on the axis nodes in the jet core (x 4.7-7.8 mm),
where the flow expands monotonically: number of local extrema and the mean |second difference|."""
import sys
import numpy as np
from paraview.simple import XMLUnstructuredGridReader
from paraview import servermanager as sm
from vtk.util.numpy_support import vtk_to_numpy as v2n
for f in sys.argv[1:]:
    g = sm.Fetch(XMLUnstructuredGridReader(FileName=[f]))
    p = v2n(g.GetPoints().GetData()) * 1e3; m = v2n(g.GetPointData().GetArray("Mach"))
    on = (np.abs(p[:, 1]) < 1e-6) & (p[:, 0] > 4.7) & (p[:, 0] < 7.8); o = np.argsort(p[on, 0])
    x, mm = p[on, 0][o], m[on][o]
    d = np.diff(mm); ext = int(np.sum(np.sign(d[1:]) != np.sign(d[:-1])))
    print(f"{f.split('/')[-4][-10:]:>10s} {f.split('/')[-2]:12s} axis nodes {len(x):4d}  local extrema {ext:4d}"
          f"  mean|d2M| {np.mean(np.abs(np.diff(mm, 2))):.4f}  Mach at x=4.7..7.8: {mm[0]:.3f} -> {mm[-1]:.3f}"
          f"  min after the 4.5 mm shock {mm[x < 5.5].min():.3f}")
PY
unset LD_PRELOAD; module load paraview/5.13 >/dev/null 2>&1
FF=$F/2_jet_staged/f3_cap3/flow.vtu
echo "== inspect fine"; pvbatch --mesa $C/inspect_field.py $FF fine_f3.npz 2>&1 | grep -vi warning | grep -E "bow shock|backflow|Mach on outlet"
echo "== bow shock on the axis"; pvbatch --mesa axis_probe.py $C/2_jet_staged/s6_cap3/flow.vtu $MD/2_jet_staged/m3_cap3/flow.vtu $FF 2>&1 | grep -vi warning | grep -E "flow.vtu|line:"
echo "== jet core on the axis (wiggles)"; pvbatch --mesa axis_wiggle.py $C/2_jet_staged/s6_cap3/flow.vtu $MD/2_jet_staged/m3_cap3/flow.vtu $FF 2>&1 | grep -vi warning
echo "== interpolation points still off vs medium m3 (|dMach| > 0.5)"; pvbatch --mesa bad_points.py $MD/2_jet_staged/m3_cap3/flow.vtu $FF 2>&1 | grep -vi warning | head -10
echo "== field difference fine f3 vs medium m3"; pvbatch --mesa diff_field.py $MD/2_jet_staged/m3_cap3/flow.vtu $FF 2>&1 | grep -vi warning
[ -f $MD/analysis/medium_m3.npz ] || pvbatch --mesa $C/inspect_field.py $MD/2_jet_staged/m3_cap3/flow.vtu $MD/analysis/medium_m3.npz > /dev/null 2>&1
module unload paraview/5.13 >/dev/null 2>&1
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh >/dev/null 2>&1
ln -sf $MD/0_interp/coarse_s6.npz coarse_s6.npz; ln -sf $MD/analysis/medium_m3.npz medium_m3.npz
python $C/plot_field.py coarse_medium_fine_Mach.png coarse_s6.npz medium_m3.npz fine_f3.npz
for s in f1_cap0.5 f2_cap1 f3_cap3; do cp $F/2_jet_staged/$s/convergence.png ${s%%_*}_convergence.png; done
ls -l
