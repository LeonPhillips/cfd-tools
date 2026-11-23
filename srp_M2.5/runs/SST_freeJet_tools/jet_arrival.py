"""pvbatch jet_arrival.py <flow.vtu> [x_mm=220] [u_min=20]
Has the free jet reached the outlet? (user, 2026-11-23) Probes the axial velocity just off the axis (r = 0.01 mm) at
x = x_mm and prints one line, starting "ARRIVED" if Ux > u_min [m/s], else "not yet". run_fj_chain.sh uses it to end
t1: once the jet has arrived, t1 runs one more 1000-iteration chunk. It does not stop SU2 itself.
"""
import sys
from paraview.simple import XMLUnstructuredGridReader, ProbeLocation
from paraview import servermanager as sm

vtu = sys.argv[1]
x_mm = float(sys.argv[2]) if len(sys.argv) > 2 else 220.0
u_min = float(sys.argv[3]) if len(sys.argv) > 3 else 20.0
r = XMLUnstructuredGridReader(FileName=[vtu])
p = ProbeLocation(Input=r, ProbeType="Fixed Radius Point Source")
p.ProbeType.Center = [x_mm * 1e-3, 1e-5, 0.0]
d = sm.Fetch(p)
pd = d.GetPointData()
mask = pd.GetArray("vtkValidPointMask")
if mask is not None and mask.GetTuple1(0) == 0:
    sys.exit(f"probe point x = {x_mm:g} mm, r = 0.01 mm is outside the mesh")
ux = pd.GetArray("Velocity").GetTuple3(0)[0]
mach = pd.GetArray("Mach").GetTuple1(0)
state = "ARRIVED" if ux > u_min else "not yet"
print(f"{state}: axis Ux {ux:.2f} m/s (Mach {mach:.3f}) at x = {x_mm:g} mm, threshold {u_min:g} m/s")
