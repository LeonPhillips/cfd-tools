"""pvbatch --mesa resample_box.py <out.npz> <label=flow.vtu> ...: Mach and Density on a uniform grid over the refinement
box x 0.5-20 mm, r 0-9.5 mm at 0.01 mm (1951 x 951 points). Points outside the fluid are NaN."""
import sys
import numpy as np
from paraview.simple import XMLUnstructuredGridReader, ResampleToImage
from paraview import servermanager as sm
from vtk.util.numpy_support import vtk_to_numpy as v2n
X0, X1, R0, R1, H = 0.5, 20.0, 0.0, 9.5, 0.01
NX, NR = int(round((X1 - X0) / H)) + 1, int(round((R1 - R0) / H)) + 1
out = {"x_mm": np.linspace(X0, X1, NX), "r_mm": np.linspace(R0, R1, NR), "labels": []}
for arg in sys.argv[2:]:
    lab, f = arg.split("=", 1)
    rs = ResampleToImage(Input=XMLUnstructuredGridReader(FileName=[f]))
    rs.UseInputBounds = 0
    rs.SamplingBounds = [X0 * 1e-3, X1 * 1e-3, R0 * 1e-3, R1 * 1e-3, 0.0, 0.0]
    rs.SamplingDimensions = [NX, NR, 1]
    d = sm.Fetch(rs); pd = d.GetPointData()
    ok = v2n(pd.GetArray("vtkValidPointMask")).astype(bool)
    for v in ("Mach", "Density"):
        a = v2n(pd.GetArray(v)).astype(np.float32); a[~ok] = np.nan
        out[f"{lab}_{v}"] = a.reshape(NR, NX)   # rows r, columns x
    out["labels"].append(lab)
    print(f"{lab}: {f}  grid {NX} x {NR}, fluid points {int(ok.sum())} of {ok.size}")
out["labels"] = np.array(out["labels"])
np.savez_compressed(sys.argv[1], **out)
print("saved", sys.argv[1])
