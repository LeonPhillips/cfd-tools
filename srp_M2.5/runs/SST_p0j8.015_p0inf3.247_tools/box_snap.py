"""pvbatch --mesa box_snap.py <out dir> <flow_NNNNN.vtu> ...: Mach and Density on a uniform 0.01 mm grid over the jet box
x 1.5-8.5 mm, r 0-3.6 mm (701 x 361 points), one <out dir>/<vtu name>.npz per snapshot. Points outside the fluid are NaN.
Used by watch_fj.sh for the box-change panel of plot_conv_roi.py."""
import os
import sys
import numpy as np
from paraview.simple import XMLUnstructuredGridReader, ResampleToImage, Delete
from paraview import servermanager as sm
from vtk.util.numpy_support import vtk_to_numpy as v2n
X0, X1, R0, R1, H = 1.5, 8.5, 0.0, 3.6, 0.01
NX, NR = int(round((X1 - X0) / H)) + 1, int(round((R1 - R0) / H)) + 1
for f in sys.argv[2:]:
    name = os.path.basename(f)[:-4]
    reader = XMLUnstructuredGridReader(FileName=[f])
    rs = ResampleToImage(Input=reader)
    rs.UseInputBounds = 0
    rs.SamplingBounds = [X0 * 1e-3, X1 * 1e-3, R0 * 1e-3, R1 * 1e-3, 0.0, 0.0]
    rs.SamplingDimensions = [NX, NR, 1]
    d = sm.Fetch(rs); pd = d.GetPointData()
    ok = v2n(pd.GetArray("vtkValidPointMask")).astype(bool)
    out = {"x_mm": np.linspace(X0, X1, NX), "r_mm": np.linspace(R0, R1, NR), "iter": int(name.split("_")[-1])}
    for v in ("Mach", "Density"):
        a = v2n(pd.GetArray(v)).astype(np.float32); a[~ok] = np.nan
        out[v] = a.reshape(NR, NX)   # rows r, columns x
    np.savez_compressed(f"{sys.argv[1]}/{name}.npz", **out)
    print(f"{name}: fluid points {int(ok.sum())} of {ok.size}")
    Delete(rs); Delete(reader)
