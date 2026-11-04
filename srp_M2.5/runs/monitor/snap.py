"""pvbatch --mesa snap.py <out dir> <flow_NNNNN.vtu> ...: one <out dir>/flow_NNNNN.npz per snapshot (user, 2026-11-04):
- box: Mach, Density, Velocity_x on a uniform 0.025 mm grid, x 0-20 mm, r 0-9.5 mm (rows r, columns x; NaN outside
  the fluid). Used for the difference map (cropped to x 0-16, r 0-6), box_change.csv and the feature positions.
- lines: Mach along the convergence lines (LINES below), 2000 points each.
Used by watch_conv.sh. Keep LINES in step with conv_update.py and plot_conv_lines.py (they read lines_mm from the npz).
"""
import os
import sys
import numpy as np
from paraview.simple import XMLUnstructuredGridReader, ResampleToImage, PlotOverLine, Delete
from paraview import servermanager as sm
from vtk.util.numpy_support import vtk_to_numpy as v2n

X0, X1, R0, R1, H = 0.0, 20.0, 0.0, 9.5, 0.025
NX, NR = int(round((X1 - X0) / H)) + 1, int(round((R1 - R0) / H)) + 1
# (x0, r0, x1, r1) mm, user 2026-11-04: r = 0.01, 1.3, 2.3, 3, 5 mm from x = 1 to 20 mm, plus the jet shear layer
LINES = [(1.0, r, 20.0, r) for r in (0.01, 1.3, 2.3, 3.0, 5.0)] + [(2.0, 2.34, 7.84, 4.70)]
NAMES = ["r=0.01", "r=1.3", "r=2.3", "r=3", "r=5", "shear layer"]
N = 2000


def valid(pd, name):
    ok = v2n(pd.GetArray("vtkValidPointMask")).astype(bool)
    a = v2n(pd.GetArray(name)).astype(np.float32)
    if a.ndim > 1:
        a = a[:, 0]
    a[~ok] = np.nan
    return a, ok


for f in sys.argv[2:]:
    name = os.path.basename(f)[:-4]
    rd = XMLUnstructuredGridReader(FileName=[f])
    out = {"x_mm": np.linspace(X0, X1, NX), "r_mm": np.linspace(R0, R1, NR), "iter": int(name.split("_")[-1]),
           "lines_mm": np.array(LINES), "line_names": np.array(NAMES)}
    rs = ResampleToImage(Input=rd)
    rs.UseInputBounds = 0
    rs.SamplingBounds = [X0 * 1e-3, X1 * 1e-3, R0 * 1e-3, R1 * 1e-3, 0.0, 0.0]
    rs.SamplingDimensions = [NX, NR, 1]
    pd = sm.Fetch(rs).GetPointData()
    for v, key in (("Mach", "Mach"), ("Density", "Density"), ("Velocity", "Velocity_x")):
        a, ok = valid(pd, v)
        out[key] = a.reshape(NR, NX)
    Delete(rs)
    L = []
    for x0, r0, x1, r1 in LINES:
        pol = PlotOverLine(Input=rd); pol.Point1 = [x0 * 1e-3, r0 * 1e-3, 0]; pol.Point2 = [x1 * 1e-3, r1 * 1e-3, 0]
        pol.Resolution = N - 1
        a, _ = valid(sm.Fetch(pol).GetPointData(), "Mach")
        L.append(a); Delete(pol)
    out["line_Mach"] = np.array(L)
    np.savez_compressed(f"{sys.argv[1]}/{name}.npz", **out)
    print(f"{name}: box fluid points {int(ok.sum())} of {ok.size}, line points {[int(np.isfinite(a).sum()) for a in L]}")
    Delete(rd)
