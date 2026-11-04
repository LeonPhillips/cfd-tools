"""pvbatch --mesa probes_md.py <flow.vtu of the starting solution>: prints the CUSTOM_OUTPUTS value for six Mach probes
either side of the Mach disk (user default, 2026-11-04), for a stage config:
    CUSTOM_OUTPUTS= '<printed line>'      (and add CUSTOM to HISTORY_OUTPUT)
The Mach disk (normal jet shock, ~8.6 mm) is found in the starting solution with grid_compare/feature_positions.py
(0.025 mm box). At r = 0, 1 and 2 mm: up<r> = 0.3 mm upstream of the disk (supersonic side), dn<r> = 0.1 mm behind it
(subsonic side). On the axis dn0 is capped at half-way to the stagnation point, which can be < 0.1 mm behind the disk.
Probe Mach is |velocity| / sound speed (SU2 rejects Probe{MACH}). Positions are printed in m, as SU2 wants.
"""
import os
import sys
import numpy as np
from paraview.simple import XMLUnstructuredGridReader, ResampleToImage
from paraview import servermanager as sm
from vtk.util.numpy_support import vtk_to_numpy as v2n

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "grid_compare"))
import feature_positions as fp   # noqa: E402

UP, DN, RS = 0.3, 0.1, (0.0, 1.0, 2.0)
X0, X1, R0, R1, H = 0.0, 20.0, 0.0, 9.5, 0.025
NX, NR = int(round((X1 - X0) / H)) + 1, int(round((R1 - R0) / H)) + 1
rs = ResampleToImage(Input=XMLUnstructuredGridReader(FileName=[sys.argv[1]]))
rs.UseInputBounds = 0
rs.SamplingBounds = [X0 * 1e-3, X1 * 1e-3, R0 * 1e-3, R1 * 1e-3, 0.0, 0.0]
rs.SamplingDimensions = [NX, NR, 1]
pd = sm.Fetch(rs).GetPointData()
ok = v2n(pd.GetArray("vtkValidPointMask")).astype(bool)
A = {}
for v in ("Mach", "Density", "Velocity"):
    a = v2n(pd.GetArray(v)).astype(float)
    a = a[:, 0] if a.ndim > 1 else a
    a[~ok] = np.nan
    A[v] = a.reshape(NR, NX)
fp.set_grid(np.linspace(X0, X1, NX), np.linspace(R0, R1, NR))
o = fp.features(A["Mach"], A["Density"], A["Velocity"])
if "_mach_disk_track" not in o:
    sys.exit(f"no Mach disk found: {o['flags']}")
tr, tx = o["_mach_disk_track"]
M = "sqrt(pow(VELOCITY_X, 2) + pow(VELOCITY_Y, 2)) / SOUND_SPEED"
parts = []
for r in RS:
    xd = float(np.interp(r, tr, tx))
    dn = min(DN, 0.5 * (o["stagnation_x"] - xd)) if r == 0 else DN
    for name, xp in ((f"up{int(r)}", xd - UP), (f"dn{int(r)}", xd + dn)):
        parts.append(f"{name} : Probe{{{M}}}[{xp * 1e-3:.6f}, {r * 1e-3:.6f}]")
    sys.stderr.write(f"r = {r} mm: Mach disk x = {xd:.3f} mm, probes at {xd - UP:.3f} and {xd + dn:.3f} mm\n")
sys.stderr.write(f"axis: Mach disk {o['mach_disk_x']:.3f}, stagnation {o['stagnation_x']:.3f} mm\n")
print("; ".join(parts))
