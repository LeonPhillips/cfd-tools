"""pvbatch --mesa diff_field.py <coarse flow.vtu> <interpolated flow.vtu>: Mach difference after resampling the
coarse solution onto the medium mesh points, by region."""
import sys
import numpy as np
from paraview.simple import XMLUnstructuredGridReader, ResampleWithDataset
from paraview import servermanager as sm
from vtk.util.numpy_support import vtk_to_numpy as v2n
coarse = XMLUnstructuredGridReader(FileName=[sys.argv[1]])
interp = XMLUnstructuredGridReader(FileName=[sys.argv[2]])
rs = ResampleWithDataset(SourceDataArrays=coarse, DestinationMesh=interp)
a = sm.Fetch(rs); b = sm.Fetch(interp)
xy = v2n(b.GetPoints().GetData())[:, :2] * 1e3
mc, mi = v2n(a.GetPointData().GetArray("Mach")), v2n(b.GetPointData().GetArray("Mach"))
ok = v2n(a.GetPointData().GetArray("vtkValidPointMask")).astype(bool)
d = np.abs(mi - mc)
regions = {"refined box x0.5-20 r<9.5": (xy[:, 0] > 0.5) & (xy[:, 0] < 20) & (xy[:, 1] < 9.5),
           "nozzle x<0.5 r<3": (xy[:, 0] < 0.5) & (xy[:, 0] > -8.2) & (xy[:, 1] < 3),
           "upstream x>20": xy[:, 0] > 20, "body side + wake x<-1": xy[:, 0] < -1, "all": np.ones(len(d), bool)}
print(f"{'region':28s} {'points':>8s} {'mean|dM|':>9s} {'p99|dM|':>9s} {'max|dM|':>8s}  where max (x, r mm)")
for k, m in regions.items():
    m = m & ok; i = np.argmax(np.where(m, d, -1))
    print(f"{k:28s} {m.sum():8d} {d[m].mean():9.4f} {np.percentile(d[m], 99):9.4f} {d[m].max():8.3f}  ({xy[i,0]:.2f}, {xy[i,1]:.2f})")
