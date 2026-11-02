"""pvbatch --mesa blip_probe.py <flow.vtu> [<other flow.vtu>]: interior points whose Mach differs from every
neighbour by > 0.5 (isolated spikes), with all primitive values; the second file is sampled at the same spots."""
import sys
import numpy as np
from paraview.simple import XMLUnstructuredGridReader
from paraview import servermanager as sm
from vtk.util.numpy_support import vtk_to_numpy as v2n
g = sm.Fetch(XMLUnstructuredGridReader(FileName=[sys.argv[1]]))
a = lambda k: v2n(g.GetPointData().GetArray(k))
p = v2n(g.GetPoints().GetData()) * 1e3; m = a("Mach")
# neighbours from cell connectivity
nb = [set() for _ in range(len(m))]
for c in range(g.GetNumberOfCells()):
    ids = [g.GetCell(c).GetPointId(k) for k in range(g.GetCell(c).GetNumberOfPoints())]
    for i in ids: nb[i].update(ids)
yp = a("Y_Plus")
spikes = [i for i in range(len(m)) if yp[i] == 0 and len(nb[i]) > 1 and
          min(abs(m[i] - m[j]) for j in nb[i] if j != i) > 0.5]
print(sys.argv[1], "isolated Mach spikes (off walls):", len(spikes))
for i in spikes[:30]:
    v = a("Velocity")[i]
    print(f"   x {p[i,0]:8.3f} r {p[i,1]:7.4f}  Mach {m[i]:.4f} (neighbours {np.median([m[j] for j in nb[i] if j != i]):.3f})"
          f"  rho {a('Density')[i]:.4f} p {a('Pressure')[i]:.1f} T {a('Temperature')[i]:.2f} u {v[0]:.2f} nuT {a('Nu_Tilde')[i]:.3e}")
