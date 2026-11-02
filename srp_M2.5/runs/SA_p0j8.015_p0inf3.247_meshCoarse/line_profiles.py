"""pvbatch --mesa line_profiles.py <out.npz> <label=flow.vtu> [...]: Mach and Density along straight lines, 2000 points
each. Points outside the fluid are NaN. Default lines (standard convergence check, user 2026-10-28): r = 0.01, 1.3, 2.3,
3, 6, 9 mm from x = 1 to 20 mm. Other lines: LINES="x0,r0,x1,r1;..." in mm (env var)."""
import os
import sys
import numpy as np
from paraview.simple import XMLUnstructuredGridReader, PlotOverLine
from paraview import servermanager as sm
from vtk.util.numpy_support import vtk_to_numpy as v2n
DEFAULT = [(1.0, r, 20.0, r) for r in (0.01, 1.3, 2.3, 3.0, 6.0, 9.0)]
LINES = [tuple(float(v) for v in s.split(",")) for s in os.environ["LINES"].split(";")] if os.environ.get("LINES") else DEFAULT
N = 2000
out = {"lines_mm": np.array(LINES), "labels": []}
for arg in sys.argv[2:]:
    lab, f = arg.split("=", 1)
    rd = XMLUnstructuredGridReader(FileName=[f])
    M, D = [], []
    for x0, r0, x1, r1 in LINES:
        pol = PlotOverLine(Input=rd); pol.Point1 = [x0 * 1e-3, r0 * 1e-3, 0]; pol.Point2 = [x1 * 1e-3, r1 * 1e-3, 0]
        pol.Resolution = N
        d = sm.Fetch(pol); pd = d.GetPointData()
        ok = v2n(pd.GetArray("vtkValidPointMask")).astype(bool)
        m = v2n(pd.GetArray("Mach")).astype(float); rho = v2n(pd.GetArray("Density")).astype(float)
        m[~ok] = np.nan; rho[~ok] = np.nan
        M.append(m); D.append(rho)
    out[lab + "_Mach"] = np.array(M); out[lab + "_Density"] = np.array(D); out["labels"].append(lab)
    print(f"{lab}: {f}  valid points per line {[int(np.isfinite(m).sum()) for m in M]}")
out["labels"] = np.array(out["labels"])
np.savez(sys.argv[1], **out)
print("saved", sys.argv[1])
