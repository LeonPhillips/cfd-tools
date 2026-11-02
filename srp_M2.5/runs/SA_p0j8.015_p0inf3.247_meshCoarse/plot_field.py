"""Mach contours and axis profiles from inspect_field.py exports. Usage: python plot_field.py out.png a.npz [b.npz ...]"""
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

out, files = sys.argv[1], sys.argv[2:]
fig, ax = plt.subplots(len(files) + 1, 1, figsize=(13, 3.6 * len(files) + 3.5))
for a, fn in zip(ax, files):
    d = np.load(fn)
    t = mtri.Triangulation(d["x"] * 1000, d["r"] * 1000)
    c = a.tricontourf(t, np.clip(d["mach"], 0, 6), levels=np.linspace(0, 6, 31), cmap="jet")
    a.set_aspect("equal"); a.set_xlim(-65, 74.1); a.set_ylim(0, 45)
    a.set_title(f"Mach: {fn}"); a.set_ylabel("r [mm]")
    fig.colorbar(c, ax=a, shrink=0.8, label="Mach (clipped at 6)")
    on_axis = d["r"] < 2e-4
    o = np.argsort(d["x"][on_axis])
    ax[-1].plot(d["x"][on_axis][o] * 1000, d["mach"][on_axis][o], lw=1, label=fn.split("/")[0])
ax[-1].set_xlim(-10, 74.1); ax[-1].set_ylim(0, 8); ax[-1].grid(alpha=0.3)
ax[-1].set_xlabel("x [mm] (freestream from the right, jet exits at x ~ 1.9 mm)"); ax[-1].set_ylabel("Mach on axis")
ax[-1].legend()
fig.tight_layout(); fig.savefig(out, dpi=100); print("saved", out)
