"""[GSA_BOX="xmin,xmax,rmax" (mm)] pvbatch --mesa export_baseflow.py <run stage dir> <out prefix>

Exports an SU2 axisymmetric solution for LinStab2D (srp_M2.5 GSA):
  <out>_flow.bin/.json  : mesh nodes inside the box (x, r in m) with Density, Velocity_x, Velocity_r,
                          Temperature, Pressure, Mach, Laminar_Viscosity, Eddy_Viscosity, Turb_Kin_Energy
                          (SI units, as written by SU2 with REF_DIMENSIONALIZATION= DIMENSIONAL)
  <out>_wall.bin/.json  : wall_noslip surface segments inside the box: x0, r0, x1, r1 (m), from surface.vtu
Binary files are float64, column-major, shape (N, K); the .json gives N and the column names.
Default box: x -12..23 mm, r 0..30 mm. (The box is an env var: pvbatch mangles negative numbers in argv.)
"""
import json
import os
import sys

import numpy as np
from paraview.simple import XMLUnstructuredGridReader
from paraview import servermanager as sm
from vtk.util.numpy_support import vtk_to_numpy

stage, out = sys.argv[1], sys.argv[2]
xmin, xmax, rmax = [float(v) * 1e-3 for v in os.environ.get("GSA_BOX", "-12,23,30").split(",")]


def fetch(reader):
    reader.UpdatePipeline()
    return sm.Fetch(reader)


def write(prefix, cols, names, meta):
    a = np.column_stack(cols).astype(np.float64)
    a.T.tofile(prefix + ".bin")          # column-major for MATLAB fread(..., [N K])
    meta = dict(meta, N=int(a.shape[0]), names=names)
    with open(prefix + ".json", "w") as f:
        json.dump(meta, f, indent=1)
    print(f"wrote {prefix}.bin: {a.shape[0]} rows x {a.shape[1]} cols")


# --- volume nodes
g = fetch(XMLUnstructuredGridReader(FileName=[stage + "/flow.vtu"]))
pts = vtk_to_numpy(g.GetPoints().GetData())
pd = g.GetPointData()
avail = [pd.GetArrayName(i) for i in range(pd.GetNumberOfArrays())]
print("point arrays:", avail)
keep = (pts[:, 0] >= xmin) & (pts[:, 0] <= xmax) & (pts[:, 1] <= rmax)
cols, names = [pts[keep, 0], pts[keep, 1]], ["x", "r"]
vel = vtk_to_numpy(pd.GetArray("Velocity"))
cols += [vel[keep, 0], vel[keep, 1]]
names += ["Velocity_x", "Velocity_r"]
for n in ["Density", "Temperature", "Pressure", "Mach", "Laminar_Viscosity", "Eddy_Viscosity", "Turb_Kin_Energy"]:
    cols.append(vtk_to_numpy(pd.GetArray(n))[keep])
    names.append(n)
write(out + "_flow", cols, names, {"source": stage + "/flow.vtu", "box_m": [xmin, xmax, 0.0, rmax]})

# --- wall segments
s = fetch(XMLUnstructuredGridReader(FileName=[stage + "/surface.vtu"]))
sp = vtk_to_numpy(s.GetPoints().GetData())
segs = []
for c in range(s.GetNumberOfCells()):
    ids = s.GetCell(c).GetPointIds()
    p = [sp[ids.GetId(k)] for k in range(ids.GetNumberOfIds())]
    for a, b in zip(p[:-1], p[1:]):
        if min(a[0], b[0]) <= xmax and max(a[0], b[0]) >= xmin and min(a[1], b[1]) <= rmax:
            segs.append([a[0], a[1], b[0], b[1]])
segs = np.array(segs)
write(out + "_wall", [segs[:, k] for k in range(4)], ["x0", "r0", "x1", "r1"],
      {"source": stage + "/surface.vtu", "box_m": [xmin, xmax, 0.0, rmax]})
