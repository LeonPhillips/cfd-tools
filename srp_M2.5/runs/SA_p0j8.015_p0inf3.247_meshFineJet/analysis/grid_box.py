"""pvbatch --mesa grid_box.py <out.npz> x0 x1 r0 r1 h label=<flow.vtu> ...: Mach and Density of each solution resampled on
one uniform grid (mm, step h), so solutions on different meshes can be subtracted. Points outside the fluid are NaN."""
import sys
import numpy as np
from paraview.simple import XMLUnstructuredGridReader, ResampleToImage, Delete
from paraview import servermanager as sm
from vtk.util.numpy_support import vtk_to_numpy as v2n
out, (X0, X1, R0, R1, H) = sys.argv[1], map(float, sys.argv[2:7])
NX, NR = int(round((X1 - X0) / H)) + 1, int(round((R1 - R0) / H)) + 1
res = {"x_mm": np.linspace(X0, X1, NX), "r_mm": np.linspace(R0, R1, NR), "labels": []}
for a in sys.argv[7:]:
    lab, f = a.split("=", 1)
    reader = XMLUnstructuredGridReader(FileName=[f])
    rs = ResampleToImage(Input=reader)
    rs.UseInputBounds = 0
    rs.SamplingBounds = [X0 * 1e-3, X1 * 1e-3, R0 * 1e-3, R1 * 1e-3, 0.0, 0.0]
    rs.SamplingDimensions = [NX, NR, 1]
    d = sm.Fetch(rs); pd = d.GetPointData()
    ok = v2n(pd.GetArray("vtkValidPointMask")).astype(bool)
    for v in ("Mach", "Density"):
        x = v2n(pd.GetArray(v)).astype(np.float32); x[~ok] = np.nan
        res[f"{lab}_{v}"] = x.reshape(NR, NX)   # rows r, columns x
    res["labels"].append(lab)
    print(f"{lab}: {f}  fluid points {int(ok.sum())} of {ok.size}")
    Delete(rs); Delete(reader)
np.savez_compressed(out, **res)
print("saved", out)
