"""python plot_lines.py <out.png> <profiles.npz> [title]: Mach (left axis, solid) and density (right axis, dashed)
along each line from line_profiles.py, one colour per solution, and a table of differences against the last solution."""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
z = np.load(sys.argv[2]); labs = [str(l) for l in z["labels"]]
if "lines_mm" in z:
    L = z["lines_mm"]
else:  # files written before 2026-10-28 18:00 hold only the standard lines as r_mm and x_mm
    L = np.array([(z["x_mm"][0], r, z["x_mm"][-1], r) for r in z["r_mm"]])
n = len(L); nc = 2 if n > 1 else 1; nr = -(-n // nc)
cols = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]


def line_info(x0, r0, x1, r1):
    """Abscissa (the coordinate that varies), its axis label and the panel title."""
    t = np.linspace(0, 1, len(z[labs[0] + "_Mach"][0]))
    if r0 == r1:
        return x0 + t * (x1 - x0), "x [mm] (freestream from the right, jet exits at x ~ 1.9 mm)", f"r = {r0:g} mm, x {x0:g} to {x1:g} mm"
    if x0 == x1:
        return r0 + t * (r1 - r0), "r [mm]", f"x = {x0:g} mm, r {r0:g} to {r1:g} mm"
    return np.hypot(x1 - x0, r1 - r0) * t, "distance along the line [mm]", f"({x0:g}, {r0:g}) to ({x1:g}, {r1:g}) mm"


fig, axs = plt.subplots(nr, nc, figsize=(7 * nc, 4 * nr), squeeze=False)
for ax in axs.flat[n:]:
    ax.set_visible(False)
info = [line_info(*l) for l in L]
for i, (ax, (s, xlab, title)) in enumerate(zip(axs.flat, info)):
    ax2 = ax.twinx()
    for c, l in zip(cols, labs):
        ax.plot(s, z[l + "_Mach"][i], "-", color=c, lw=1.1, label=f"{l} Mach")
        ax2.plot(s, z[l + "_Density"][i], "--", color=c, lw=0.9, label=f"{l} density")
    ax.set_title(title); ax.set_ylabel("Mach (solid)"); ax2.set_ylabel("density [kg/m$^3$] (dashed)")
    ax.set_xlabel(xlab); ax.grid(alpha=0.3)
    if i == 0:
        h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper left")
fig.suptitle(sys.argv[3] if len(sys.argv) > 3 else "Mach and density along lines")
fig.tight_layout(); fig.savefig(sys.argv[1], dpi=110)
ref = labs[-1]
print(f"differences against {ref}: |dMach| and |drho|/rho_ref (rho_ref = local value of {ref})")
print(f"{'line':>28s} {'solution':>10s} {'mean|dM|':>9s} {'p99|dM|':>8s} {'max|dM|':>8s} {'at':>7s}"
      f" {'mean drho%':>10s} {'p99 drho%':>9s} {'max drho%':>9s}")
for l in labs[:-1]:
    for i, (s, _, title) in enumerate(info):
        dm = np.abs(z[l + "_Mach"][i] - z[ref + "_Mach"][i])
        dr = 100 * np.abs(z[l + "_Density"][i] - z[ref + "_Density"][i]) / z[ref + "_Density"][i]
        ok = np.isfinite(dm) & np.isfinite(dr)
        name = title.split(",")[0]
        if not ok.any():
            print(f"{name:>28s} {l:>10s}  no fluid points"); continue
        dm, dr, ss = dm[ok], dr[ok], s[ok]
        print(f"{name:>28s} {l:>10s} {dm.mean():9.4f} {np.percentile(dm, 99):8.4f} {dm.max():8.4f} {ss[dm.argmax()]:7.2f}"
              f" {dr.mean():10.3f} {np.percentile(dr, 99):9.3f} {dr.max():9.2f}")
print("saved", sys.argv[1])
