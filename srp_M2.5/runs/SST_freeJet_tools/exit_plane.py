"""pvbatch exit_plane.py <flow.vtu> [x_exit_mm=1.876] [r_lip_mm=2.389] [pa_Pa=101325]
Averages of the flow across the nozzle exit plane of an srp_M2.5 solution (user, 2026-11-23), to set the ambient
pressure of the free-jet case: the area-weighted (2 pi r dr) and mass-weighted static pressure, Mach, and pe/pa for the
given pa. The plane is sampled at 2000 points from r = 0 to the lip (just inside it).
"""
import sys
import numpy as np
from paraview.simple import XMLUnstructuredGridReader, PlotOverLine
from paraview import servermanager as sm
from vtk.util.numpy_support import vtk_to_numpy as v2n

src = sys.argv[1]
xe = float(sys.argv[2]) * 1e-3 if len(sys.argv) > 2 else 1.876e-3
rl = float(sys.argv[3]) * 1e-3 if len(sys.argv) > 3 else 2.389e-3
pa = float(sys.argv[4]) if len(sys.argv) > 4 else 101325.0
line = PlotOverLine(Input=XMLUnstructuredGridReader(FileName=[src]))
line.Point1, line.Point2, line.Resolution = [xe, 0.0, 0.0], [xe, rl * 0.999, 0.0], 1999
pd = sm.Fetch(line).GetPointData()
ok = v2n(pd.GetArray("vtkValidPointMask")).astype(bool)
r = np.linspace(0.0, rl * 0.999, 2000)[ok]
p, M, rho = (v2n(pd.GetArray(n)).astype(float)[ok] for n in ("Pressure", "Mach", "Density"))
u = v2n(pd.GetArray("Velocity")).astype(float)[ok, 0]
w_a, w_m = 2 * np.pi * r, 2 * np.pi * r * rho * u
avg = lambda a, w: np.trapz(a * w, r) / np.trapz(w, r)
pe_a, pe_m = avg(p, w_a), avg(p, w_m)
print(f"source {src}")
print(f"exit plane x = {xe*1e3:.3f} mm, r 0..{rl*0.999e3:.3f} mm, {ok.sum()} valid points")
print(f"static p: axis {p[0]:.0f} Pa, at lip {p[-1]:.0f} Pa, min {p.min():.0f}, max {p.max():.0f}")
print(f"Mach: axis {M[0]:.3f}, area-avg {avg(M, w_a):.3f}, mass-avg {avg(M, w_m):.3f}")
print(f"mass flow through the plane {np.trapz(w_m, r):.6e} kg/s")
print(f"pe area-avg {pe_a:.0f} Pa, mass-avg {pe_m:.0f} Pa; pa {pa:.0f} Pa -> pe/pa {pe_a/pa:.3f} (area), {pe_m/pa:.3f} (mass)")
print(f"pa for pe/pa = 1.10 / 1.15 / 1.20 (area-avg pe): {pe_a/1.1:.0f} / {pe_a/1.15:.0f} / {pe_a/1.2:.0f} Pa")
