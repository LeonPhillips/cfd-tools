"""pvbatch snap_fj.py <out dir> <flow_<iter>.vtu> ...   (ParaView 6.1 on Gadi; no --mesa flag)
Free-jet live monitor (user, 2026-11-23): Mach along the convergence lines, 2000 points each, one npz per snapshot:
    <out dir>/flow_<iter>.npz with x_mm, r_mm (per line), mach (lines x 2000), lines_mm, line_names, iteration.
Lines (mm): r = 0.01 (axis), 1.2 (half the exit radius), 2.4 (lip) and 4.8 (one exit diameter) from x = 1 to 100,
and the shear layer from the lip outwards, (2.0, 2.45) to (40.0, 5.0). Keep in step with conv_fj.py / plot_fj.py,
which read lines_mm and line_names from the npz.
"""
import os
import re
import sys
import numpy as np
from paraview.simple import XMLUnstructuredGridReader, PlotOverLine, Delete
from paraview import servermanager as sm
from vtk.util.numpy_support import vtk_to_numpy as v2n

LINES = [(1.0, r, 100.0, r) for r in (0.01, 1.2, 2.4, 4.8)] + [(2.0, 2.45, 40.0, 5.0)]
NAMES = ["r = 0.01 mm (axis)", "r = 1.2 mm", "r = 2.4 mm (lip)", "r = 4.8 mm", "shear layer"]
N = 2000
out_dir = sys.argv[1]
for f in sys.argv[2:]:
    it = int(re.search(r"flow_(\d+)\.vtu$", f).group(1))
    rd = XMLUnstructuredGridReader(FileName=[f])
    mach = np.full((len(LINES), N), np.nan)
    for i, (x0, r0, x1, r1) in enumerate(LINES):
        pl = PlotOverLine(Input=rd)
        pl.Point1, pl.Point2, pl.Resolution = [x0 * 1e-3, r0 * 1e-3, 0.0], [x1 * 1e-3, r1 * 1e-3, 0.0], N - 1
        pd = sm.Fetch(pl).GetPointData()
        ok = v2n(pd.GetArray("vtkValidPointMask")).astype(bool)
        m = v2n(pd.GetArray("Mach")).astype(float)
        m[~ok] = np.nan
        mach[i] = m
        Delete(pl)
    Delete(rd)
    np.savez_compressed(os.path.join(out_dir, f"flow_{it:06d}.npz"), mach=mach, lines_mm=np.array(LINES),
                        line_names=np.array(NAMES), iteration=it, npts=N)
    print(f"snap {it}: axis Mach max {np.nanmax(mach[0]):.3f}")
