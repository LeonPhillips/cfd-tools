"""Field check for one phase (pvbatch --mesa). Usage: pvbatch --mesa inspect_field.py <flow.vtu> <out.npz>
Prints the axis Mach/pressure profile ahead of the body, the bow-shock position on the axis, Mach along
the freestream inlet (x = 74.08 mm) and the outflow/backflow on the outlet; saves the model region to .npz for plot_field.py.
"""
import sys

import numpy as np
from paraview.simple import *
from paraview import servermanager
from vtk.numpy_interface import dataset_adapter as dsa

src, out = sys.argv[1], sys.argv[2]
r = XMLUnstructuredGridReader(FileName=[src])
r.UpdatePipeline()
names = r.PointData.keys()
print("fields:", ", ".join(names))


def line(p1, p2, n=2000):
    pol = PlotOverLine(Input=r)
    pol.Point1, pol.Point2, pol.Resolution = p1, p2, n
    d = dsa.WrapDataObject(servermanager.Fetch(pol))
    pts = d.Points
    ok = np.array(d.PointData["vtkValidPointMask"]).astype(bool)
    return pts[ok], {k: np.array(d.PointData[k])[ok] for k in ("Mach", "Pressure", "Velocity") if k in names}



# axis, just off r = 0, from the inlet to the lip/nozzle exit
pts, f = line([0.0741, 1e-5, 0], [0.0, 1e-5, 0])
x, m, p = pts[:, 0] * 1000, f["Mach"], f["Pressure"]
print("\naxis x[mm]  Mach  p[kPa]  (every ~4 mm)")
for i in range(0, len(x), 100):
    print(f"  {x[i]:7.2f}  {m[i]:6.3f}  {p[i]/1e3:8.2f}")
drop = np.where(m < 0.97 * m[0])[0]
print(f"bow shock on axis (first point from the inlet with Mach < 97% of inlet value): "
      f"{x[drop[0]]:.2f} mm" if len(drop) else "no shock found on the axis")

# freestream inlet edge 16 (x = 74.08 mm, r 0..45 mm), sampled just inside
pts, f = line([0.07405, 1e-5, 0], [0.07405, 0.04499, 0], 450)
m = f["Mach"]
print(f"\ninlet (x = 74.05 mm): Mach min {m.min():.3f} max {m.max():.3f}; "
      f"fraction with Mach < 2.4: {np.mean(m < 2.4):.2f}; first r with Mach < 2.4: "
      + (f"{pts[np.argmax(m < 2.4), 1]*1000:.2f} mm" if np.any(m < 2.4) else "none"))

# outlet: dump cells near the outlet edge from the whole field (outlet = curve 20, upper-left block edge)
d = dsa.WrapDataObject(servermanager.Fetch(r))
P = d.Points
xmin, rmax = P[:, 0].min(), P[:, 1].max()
sel = P[:, 1] > rmax - 1e-4
v = np.array(d.PointData["Velocity"])[sel] if "Velocity" in names else None
print(f"\ndomain x {xmin*1000:.1f}..{P[:,0].max()*1000:.1f} mm, r max {rmax*1000:.1f} mm")
if v is not None:
    print(f"outlet-row points (r = r max): {sel.sum()}, x {P[sel,0].min()*1000:.1f}..{P[sel,0].max()*1000:.1f} mm; "
          f"radial velocity min {v[:,1].min():.1f} max {v[:,1].max():.1f} m/s; backflow points {np.sum(v[:,1] < 0)}")
    print(f"  Mach on outlet row: min {np.array(d.PointData['Mach'])[sel].min():.3f} "
          f"max {np.array(d.PointData['Mach'])[sel].max():.3f}")

# export the model region for plotting (rendering needs a display; matplotlib does not)
box = (P[:, 0] > -0.07) & (P[:, 0] < 0.0745) & (P[:, 1] < 0.0455)
np.savez_compressed(out, x=P[box, 0], r=P[box, 1], mach=np.array(d.PointData["Mach"])[box],
                    p=np.array(d.PointData["Pressure"])[box], u=np.array(d.PointData["Velocity"])[box, 0])
print("saved", out)
