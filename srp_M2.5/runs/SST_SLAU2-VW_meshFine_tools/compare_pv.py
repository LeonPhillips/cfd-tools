"""python compare_pv.py <ref.npz> <new.npz>: check the Gadi ParaView 6.1 monitor sampling against M3's ParaView 5.13 one
(user, 2026-11-05). Both npz come from monitor/snap.py on the same solution (p0j 7.0 bar meshFine t3, iteration 19999).
Prints max / mean |difference| of line Mach (per line), of box Mach / Density / Velocity_x, the NaN-mask agreement, and the
monitor features (grid_compare/feature_positions.py) from each. PASS if every line and box difference is < 1e-4 in Mach
and the features agree within 1e-3 mm.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "grid_compare"))
import feature_positions as fp   # noqa: E402

a, b = np.load(sys.argv[1]), np.load(sys.argv[2])
ok = True
print(f"ref {sys.argv[1]}\nnew {sys.argv[2]}")
for i, n in enumerate(a["line_names"]):
    d = np.abs(a["line_Mach"][i] - b["line_Mach"][i])
    nan = int(np.sum(np.isnan(a["line_Mach"][i]) != np.isnan(b["line_Mach"][i])))
    print(f"line {n:12s} max |dM| {np.nanmax(d):.2e}  mean {np.nanmean(d):.2e}  NaN mismatches {nan}")
    ok &= np.nanmax(d) < 1e-4 and nan == 0
for k in ("Mach", "Density", "Velocity_x"):
    d = np.abs(a[k] - b[k])
    nan = int(np.sum(np.isnan(a[k]) != np.isnan(b[k])))
    rel = np.nanmax(d) / np.nanmax(np.abs(a[k]))
    print(f"box  {k:12s} max |d| {np.nanmax(d):.2e} ({rel:.1e} of max)  mean {np.nanmean(d):.2e}  NaN mismatches {nan}")
    ok &= rel < 1e-4 and nan == 0
fp.set_grid(a["x_mm"], a["r_mm"])
fa = fp.features(a["Mach"], a["Density"], a["Velocity_x"])
fb = fp.features(b["Mach"], b["Density"], b["Velocity_x"])
for k in sorted(fa):
    if k.startswith("_") or not isinstance(fa[k], (float, int, np.floating)):
        continue
    d = abs(fa[k] - fb.get(k, np.nan))
    print(f"feature {k:24s} ref {fa[k]:.4f}  new {fb.get(k, np.nan):.4f}  |d| {d:.1e}")
    ok &= not (d > 1e-3)
print("PASS" if ok else "FAIL")
