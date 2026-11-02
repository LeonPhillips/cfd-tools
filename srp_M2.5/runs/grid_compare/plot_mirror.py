"""python plot_mirror.py <box.npz> <labelTop> <labelBottom> <out.png> [xmax_mm]: Mach contours, top label above the
axis and bottom label mirrored below it. Box from resample_box.py; NaN (solid) drawn grey."""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
f, top, bot, out = sys.argv[1:5]
xmax = float(sys.argv[5]) if len(sys.argv) > 5 else 16.0
d = np.load(f)
x, r = d["x_mm"], d["r_mm"]
ix = x <= xmax
x = x[ix]
A, B = d[f"{top}_Mach"][:, ix], d[f"{bot}_Mach"][:, ix]
R = np.concatenate([-r[::-1], r[1:]])
M = np.vstack([B[::-1], A[1:]])
levels = np.linspace(0, 4, 81)
fig, ax = plt.subplots(figsize=(12, 9.6))
ax.set_facecolor("0.55")
cf = ax.contourf(x, R, np.clip(M, 0, 4), levels=levels, cmap="turbo")
ax.contour(x, R, M, levels=[1.0], colors="k", linewidths=0.8)
ax.axhline(0, color="w", lw=0.8, ls="--")
ax.set_aspect("equal")
ax.set_xlim(x[0], xmax); ax.set_ylim(R[0], R[-1])
ax.invert_xaxis() if False else None
ax.set_xlabel("x [mm] (freestream from the right)"); ax.set_ylabel("r [mm]")
kw = dict(transform=ax.transAxes, ha="left", fontsize=13, color="w", fontweight="bold",
          bbox=dict(boxstyle="round", fc="k", alpha=0.5))
ax.text(0.01, 0.97, f"{top.split('_')[0]} (above axis)", va="top", **kw)
ax.text(0.01, 0.03, f"{bot.split('_')[0]} (below axis)", va="bottom", **kw)
cb = fig.colorbar(cf, ax=ax, shrink=0.8, ticks=np.arange(0, 4.01, 0.5)); cb.set_label("Mach")
ax.set_title(f"Mach number, meshFineJet2: {top} (top) vs {bot} (bottom); black line: Mach 1; grey: body/nozzle")
fig.tight_layout(); fig.savefig(out, dpi=200)
print("saved", out)
