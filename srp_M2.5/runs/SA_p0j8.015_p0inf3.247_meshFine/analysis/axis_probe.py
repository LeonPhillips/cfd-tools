"""pvbatch --mesa axis_probe.py <flow.vtu> [...]: Mach on the axis (r = 0.01 mm), x 14-19 mm, and the mesh points there."""
import sys
import numpy as np
from paraview.simple import XMLUnstructuredGridReader, PlotOverLine
from paraview import servermanager as sm
from vtk.util.numpy_support import vtk_to_numpy as v2n
for f in sys.argv[1:]:
    r = XMLUnstructuredGridReader(FileName=[f])
    pol = PlotOverLine(Input=r); pol.Point1 = [0.014, 1e-5, 0]; pol.Point2 = [0.019, 1e-5, 0]; pol.Resolution = 200
    d = sm.Fetch(pol); x = v2n(d.GetPoints().GetData())[:, 0] * 1e3; m = v2n(d.GetPointData().GetArray("Mach"))
    g = sm.Fetch(r); p = v2n(g.GetPoints().GetData()) * 1e3; mm = v2n(g.GetPointData().GetArray("Mach"))
    on = (np.abs(p[:, 1]) < 1e-6) & (p[:, 0] > 14) & (p[:, 0] < 19); o = np.argsort(p[on, 0])
    print(f)
    print("  line: x where Mach first < 2.425 from upstream:", x[::-1][np.argmax(m[::-1] < 2.425)].round(3))
    print("  axis mesh points (x mm, Mach) where 0.5 < Mach < 2.45:",
          [(round(a, 3), round(b, 3)) for a, b in zip(p[on, 0][o], mm[on][o]) if 0.5 < b < 2.45])
