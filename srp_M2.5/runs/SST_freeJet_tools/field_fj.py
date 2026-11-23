"""pvbatch field_fj.py <out.npz> <label>=<flow.vtu> ...   (ParaView 6.1 on Gadi)
Free-jet diagnostic (2026-11-23): resample each solution onto two uniform grids, the whole domain (x -308..240,
r 0..120 mm, 0.5 mm) and the near jet (x -20..100, r 0..40 mm, 0.1 mm), and save Ux, Ur, p and Mach for
plot_field_fj.py. Used to see where the large-scale circulation that made t1 diverge comes from.
"""
import sys
import numpy as np
from paraview.simple import XMLUnstructuredGridReader, ResampleToImage, Delete
from paraview import servermanager as sm
from vtk.util.numpy_support import vtk_to_numpy as v2n

GRIDS = {"full": (-308.0, 240.0, 0.0, 120.0, 0.5), "near": (-20.0, 100.0, 0.0, 40.0, 0.1)}
out = {}
labels = []
for arg in sys.argv[2:]:
    lab, f = arg.split("=", 1)
    labels.append(lab)
    rd = XMLUnstructuredGridReader(FileName=[f])
    for g, (x0, x1, r0, r1, h) in GRIDS.items():
        nx, nr = int(round((x1 - x0) / h)) + 1, int(round((r1 - r0) / h)) + 1
        rs = ResampleToImage(Input=rd)
        rs.UseInputBounds = 0
        rs.SamplingBounds = [x0 * 1e-3, x1 * 1e-3, r0 * 1e-3, r1 * 1e-3, 0.0, 0.0]
        rs.SamplingDimensions = [nx, nr, 1]
        pd = sm.Fetch(rs).GetPointData()
        ok = v2n(pd.GetArray("vtkValidPointMask")).astype(bool)
        sh = (nr, nx)
        v = v2n(pd.GetArray("Velocity")).astype(np.float32)
        for name, a in (("ux", v[:, 0]), ("ur", v[:, 1]), ("p", v2n(pd.GetArray("Pressure"))),
                        ("mach", v2n(pd.GetArray("Mach")))):
            a = np.asarray(a, dtype=np.float32).copy()
            a[~ok] = np.nan
            out[f"{g}_{name}_{lab}"] = a.reshape(sh)
        out[f"{g}_x"], out[f"{g}_r"] = np.linspace(x0, x1, nx), np.linspace(r0, r1, nr)
        Delete(rs)
    Delete(rd)
    print(f"{lab}: sampled {f}")
np.savez_compressed(sys.argv[1], labels=np.array(labels), **out)
