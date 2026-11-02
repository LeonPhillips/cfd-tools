"""pvbatch --mesa bad_points.py <coarse flow.vtu> <interpolated flow.vtu>: interpolated points whose Mach differs
from the coarse field (linearly resampled at the same spot) by > 0.5, grouped into clusters, with sample values."""
import sys
import numpy as np
from paraview.simple import XMLUnstructuredGridReader, ResampleWithDataset
from paraview import servermanager as sm
from vtk.util.numpy_support import vtk_to_numpy as v2n
c = XMLUnstructuredGridReader(FileName=[sys.argv[1]]); i = XMLUnstructuredGridReader(FileName=[sys.argv[2]])
a = sm.Fetch(ResampleWithDataset(SourceDataArrays=c, DestinationMesh=i)); b = sm.Fetch(i)
A = lambda d, k: v2n(d.GetPointData().GetArray(k))
p = v2n(b.GetPoints().GetData())[:, :2] * 1e3
mc, mi = A(a, "Mach"), A(b, "Mach"); ok = A(a, "vtkValidPointMask").astype(bool)
bad = np.flatnonzero(ok & (np.abs(mi - mc) > 0.5))
fs = np.flatnonzero(ok & (np.abs(A(b, "Pressure") - 19003.9) < 0.05) & (np.abs(A(b, "Velocity")[:, 0] + 567.16) < 0.005)
                    & (np.abs(A(a, "Pressure") - 19003.9) > 100))
print("points with |dMach| > 0.5:", len(bad), "of", len(mi))
print("points holding exact freestream values where the coarse field is not freestream:", len(fs))
# cluster bad points by 0.3 mm grid cells
keys = {}
for k in bad: keys.setdefault((round(p[k, 0] / 0.3), round(p[k, 1] / 0.3)), []).append(k)
groups = {}
for (gx, gy), ks in keys.items(): groups.setdefault((gx // 5, gy // 5), []).extend(ks)
print(f"{len(groups)} clusters (1.5 mm bins): x, r [mm], points, Mach interp vs coarse (median)")
for g, ks in sorted(groups.items(), key=lambda kv: -len(kv[1]))[:20]:
    ks = np.array(ks)
    print(f"   x {p[ks,0].min():7.2f}..{p[ks,0].max():7.2f}  r {p[ks,1].min():6.2f}..{p[ks,1].max():6.2f}  n {len(ks):4d}"
          f"  M {np.median(mi[ks]):.3f} vs {np.median(mc[ks]):.3f}")
