"""Checks a finished Phase 1 before Phase 2 is started (pvbatch --mesa check_phase1.py <1_freestream folder>).
PASS needs: relrms[Rho] <= -2 (2-order drop), stopped before ITER, mass imbalance < 0.5 %, freestream inflow within
0.1 % of 1.86533 kg/s, and Mach >= 1 everywhere on the outlet (x = x min), because the outlet is a supersonic outlet.
Also prints the bow-shock position on the axis. Last line: RESULT PASS or RESULT FAIL.
"""
import csv
import sys

import numpy as np
from paraview.simple import XMLUnstructuredGridReader
from paraview import servermanager
from vtk.numpy_interface import dataset_adapter as dsa

D = sys.argv[1].rstrip("/")
rows = list(csv.reader(open(f"{D}/history.csv")))
head = [c.strip().strip('"') for c in rows[0]]
h = np.array([[float(v) for v in r] for r in rows[1:] if len(r) == len(head)])
c = {k: h[:, i] for i, k in enumerate(head)}
it = int(c["Inner_Iter"][-1])
rel = c["relrms[Rho]"][-1]
m_in, m_out = c["Avg_Massflow(inlet_freestream)"][-1], c["Avg_Massflow(outlet)"][-1]
imb = 100 * abs(m_in + m_out) / m_in
dev = 100 * abs(m_in / 1.86533 - 1)

r = XMLUnstructuredGridReader(FileName=[f"{D}/flow.vtu"]); r.UpdatePipeline()
d = dsa.WrapDataObject(servermanager.Fetch(r))
P, M = np.array(d.Points), np.array(d.PointData["Mach"])
out = P[:, 0] < P[:, 0].min() + 1e-7
ax = P[:, 1] < 1e-7
xa, ma = P[ax, 0], M[ax]
o = np.argsort(-xa)
xa, ma = xa[o], ma[o]
shock = xa[np.argmax(ma < 0.97 * ma[0])] * 1000

checks = [
    (f"iterations {it} (cap 20000), relrms[Rho] {rel:.3f} (need <= -2)", rel <= -2 and it < 19999),
    (f"mass: in {m_in:.5f} out {m_out:.5f} kg/s, imbalance {imb:.3f} % (need < 0.5)", imb < 0.5),
    (f"freestream inflow off its fixed value by {dev:.4f} % (need < 0.1)", dev < 0.1),
    (f"outlet x = {P[out,0].min()*1000:.2f} mm: {out.sum()} points, Mach {M[out].min():.3f}..{M[out].max():.3f} "
     f"(need >= 1 everywhere; subsonic points: {np.sum(M[out] < 1)})", M[out].min() >= 1.0),
]
for text, ok in checks:
    print(("ok    " if ok else "FAIL  ") + text)
print(f"info  bow shock on the axis (Mach < 97 % of inlet value) at x = {shock:.2f} mm; CD = {c['CD'][-1]:.5f}")
print("RESULT " + ("PASS" if all(ok for _, ok in checks) else "FAIL"))
