"""pvbatch --mesa zero_probe.py <flow.vtu>: points with zero velocity that are not on a no-slip wall (wall nodes
legitimately have Mach 0), with their density/pressure/temperature."""
import sys
import numpy as np
from paraview.simple import XMLUnstructuredGridReader
from paraview import servermanager as sm
from vtk.util.numpy_support import vtk_to_numpy as v2n
for f in sys.argv[1:]:
    g = sm.Fetch(XMLUnstructuredGridReader(FileName=[f]))
    a = lambda k: v2n(g.GetPointData().GetArray(k))
    p = v2n(g.GetPoints().GetData()) * 1e3; m, rho, pr, t = a("Mach"), a("Density"), a("Pressure"), a("Temperature")
    cf = np.linalg.norm(a("Skin_Friction_Coefficient"), axis=1) if g.GetPointData().GetArray("Skin_Friction_Coefficient") else 0
    yp = a("Y_Plus")
    z = (m < 1e-9)
    wall = yp != 0   # Y_Plus is only non-zero on no-slip wall nodes
    bad = z & ~wall
    print(f, "Mach=0 points:", z.sum(), " on walls:", (z & wall).sum(), " off walls:", bad.sum())
    for i in np.flatnonzero(bad)[:40]:
        print(f"   x {p[i,0]:9.3f}  r {p[i,1]:7.4f}  rho {rho[i]:.4f}  p {pr[i]:9.1f}  T {t[i]:7.2f}")
