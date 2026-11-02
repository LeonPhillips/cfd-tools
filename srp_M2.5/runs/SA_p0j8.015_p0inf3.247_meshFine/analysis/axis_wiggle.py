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
