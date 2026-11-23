"""python plot_field_fj.py <in.npz> <out prefix> "<title>"
Plots field_fj.py output (2026-11-23): for each solution, the whole domain coloured by Ux with streamlines, and
p - pa; then the near jet coloured by Mach with streamlines. Writes <prefix>_full.png and <prefix>_near.png.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = np.load(sys.argv[1])
pre, title = sys.argv[2], sys.argv[3]
labs = list(d["labels"])
PA = 101325.0


def stream(ax, x, r, ux, ur, density):
    u, w = np.nan_to_num(ux), np.nan_to_num(ur)
    ax.streamplot(x, r, u, w, color="k", linewidth=0.4, density=density, arrowsize=0.5)


# whole domain: Ux + streamlines | p - pa
x, r = d["full_x"], d["full_r"]
fig, axs = plt.subplots(len(labs), 2, figsize=(16, 2.1 * len(labs) + 0.8), squeeze=False, constrained_layout=True)
for i, lab in enumerate(labs):
    ux, ur, p = d[f"full_ux_{lab}"], d[f"full_ur_{lab}"], d[f"full_p_{lab}"]
    a = axs[i, 0]
    im = a.pcolormesh(x, r, ux, cmap="RdBu_r", vmin=-30, vmax=30, shading="auto")
    stream(a, x, r, ux, ur, (3, 1))
    a.set_aspect("equal"); a.set_ylabel(f"{lab}\nr [mm]")
    fig.colorbar(im, ax=a, label="Ux [m/s]", shrink=0.9)
    b = axs[i, 1]
    im = b.pcolormesh(x, r, p - PA, cmap="RdBu_r", vmin=-1500, vmax=1500, shading="auto")
    b.set_aspect("equal")
    fig.colorbar(im, ax=b, label="p - pa [Pa]", shrink=0.9)
    sp = np.nanmax(np.hypot(ux, ur)[r > 60, :]) if np.any(r > 60) else np.nan
    b.set_title(f"max |U| for r > 60 mm: {sp:.1f} m/s", fontsize=8)
axs[-1, 0].set_xlabel("x [mm]"); axs[-1, 1].set_xlabel("x [mm]")
fig.suptitle(title + ": whole domain (Ux clipped at ±30 m/s, p - pa at ±1.5 kPa)")
fig.savefig(pre + "_full.png", dpi=110)

# near jet: Mach + streamlines
x, r = d["near_x"], d["near_r"]
fig, axs = plt.subplots(len(labs), 1, figsize=(12, 2.6 * len(labs) + 0.8), squeeze=False, constrained_layout=True)
for i, lab in enumerate(labs):
    m, ux, ur = d[f"near_mach_{lab}"], d[f"near_ux_{lab}"], d[f"near_ur_{lab}"]
    a = axs[i, 0]
    im = a.pcolormesh(x, r, m, cmap="viridis", vmin=0, vmax=2.5, shading="auto")
    stream(a, x, r, ux, ur, (3, 1))
    a.set_aspect("equal"); a.set_ylabel(f"{lab}\nr [mm]")
    fig.colorbar(im, ax=a, label="Mach", shrink=0.9)
axs[-1, 0].set_xlabel("x [mm]")
fig.suptitle(title + ": near jet")
fig.savefig(pre + "_near.png", dpi=110)
print("wrote", pre + "_full.png", pre + "_near.png")
