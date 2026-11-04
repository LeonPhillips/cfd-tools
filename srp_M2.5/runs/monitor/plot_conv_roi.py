"""Convergence plot for one stage (monitor copy, 2026-11-04, of the SST version made 2026-10-30), focused on the jet region
(user, 2026-10-29: Mach in the region of interest instead of CD; CD is printed in the title for information only).
Changes in this copy (user, 2026-11-04): the difference map covers x 0-16, r 0-6 mm (snapshots from monitor/snap.py, box
x 0-20, r 0-9.5 mm); box_change.csv is over the same region; the probes are read from the stage config's CUSTOM_OUTPUTS
(placed either side of the Mach disk by probes_md.py), not a fixed list.

Usage (on the desktop node): python plot_conv_roi.py <stage folder> [previous stage folder]
Writes <stage folder>/convergence.png with six panels:
  residuals  | probe Mach, as change from its current value | box change per 1000 iterations (box_change.csv)
  mass flows | average CFL                                  | Mach difference map: latest minus previous snapshot
The snapshots come one every 1000 iterations; conv_update.py keeps only the latest two. While a stage has fewer than two, the previous stage's latest snapshot stands in.
"""
import csv
import glob
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MAP = (0.0, 16.0, 0.0, 6.0)   # difference map region, mm (user, 2026-11-04)


def read_probes(folder):
    """{name: (x mm, r mm, label)} from CUSTOM_OUTPUTS in the stage's .cfg. Names from probes_md.py: up<r>/dn<r> =
    upstream/downstream of the Mach disk at r = <r> mm."""
    out = {}
    for cfg in glob.glob(f"{folder}/*.cfg"):
        for line in open(cfg):
            if line.startswith("CUSTOM_OUTPUTS="):
                for n, x, r in re.findall(r"(\w+)\s*:\s*Probe\{[^}]*\}\[\s*([-0-9.eE]+)\s*,\s*([-0-9.eE]+)\s*\]", line):
                    m = re.match(r"(up|dn)(\d+)$", n)
                    what = (f"{'upstream of' if m.group(1) == 'up' else 'behind'} the Mach disk" if m else "probe")
                    out[n] = (float(x) * 1e3, float(r) * 1e3, what)
    return out


folder = sys.argv[1].rstrip("/")
prev = sys.argv[2].rstrip("/") if len(sys.argv) > 2 and sys.argv[2] != "-" else None
PROBES = read_probes(folder)
rows = list(csv.reader(open(f"{folder}/history.csv")))
head = [c.strip().strip('"') for c in rows[0]]
data = np.array([[float(v) for v in r] for r in rows[1:] if len(r) == len(head)])
if data.ndim < 2 or len(data) < 2:
    sys.exit("history.csv has fewer than 2 rows yet")
col = {h: data[:, i] for i, h in enumerate(head)}
it = col["Inner_Iter"]
flows = [h for h in head if h.startswith("Avg_Massflow(")]
cfl = col["Avg CFL"] if "Avg CFL" in col else col["AvgCFL"]
probes = [p for p in PROBES if p in col]

fig, ax = plt.subplots(2, 3, figsize=(19, 9))
for k in [k for k in ["rms[Rho]", "rms[RhoU]", "rms[RhoE]", "rms[nu]", "rms[k]", "rms[w]"] if k in col]:
    ax[0, 0].plot(it, col[k], lw=0.8, label=k)
ax[0, 0].plot(it, col["relrms[Rho]"], "k--", lw=1, label="relrms[Rho] (drop from peak)")
ax[0, 0].set_ylabel("log10 residual"); ax[0, 0].legend(fontsize=8)

for p in probes:
    x, r, what = PROBES[p]
    ax[0, 1].plot(it, col[p] - col[p][-1], lw=0.9, label=f"{p} x {x:.2f} r {r:.2f}: {what}  (now {col[p][-1]:.4f})")
ax[0, 1].set_ylabel("probe Mach - current value"); ax[0, 1].legend(fontsize=7)
ax[0, 1].set_title("Mach at the probes (each shown relative to its latest value)", fontsize=9)

bc = f"{folder}/box_change.csv"
p = np.atleast_2d(np.genfromtxt(bc, delimiter=",", skip_header=1)) if os.path.exists(bc) and os.path.getsize(bc) > 60 else None
if p is not None and p.size:
    ax[0, 2].semilogy(p[:, 0], p[:, 2], "o-", ms=3, label="mean |dMach|")
    ax[0, 2].semilogy(p[:, 0], p[:, 3], "s-", ms=3, label="p99 |dMach| (99 % of points change less)")
    top = p[:, 5] if p.shape[1] > 5 else np.full(len(p), np.nan)
    if np.isfinite(top).any():
        ax[0, 2].semilogy(p[:, 0], top, "D-", ms=3, label="mean |dMach| of the top 1 % (>= p99)")
    ax[0, 2].semilogy(p[:, 0], p[:, 4], "^--", ms=3, label="mean |drho|/rho [%]")
    ax[0, 2].legend(fontsize=8)
    ax[0, 2].set_title(f"jet box x 1.5-8.5, r 0-3.6 mm: change per 1000 iterations\n"
                       f"latest: mean |dM| {p[-1, 2]:.2e}, p99 {p[-1, 3]:.2e}, top-1 % mean {top[-1]:.2e}, "
                       f"mean |drho|/rho {p[-1, 4]:.2e} %", fontsize=9)
else:
    ax[0, 2].text(0.5, 0.5, "no snapshot pair yet", ha="center", transform=ax[0, 2].transAxes)
ax[0, 2].set_ylabel("change per 1000 iterations")

for k in flows:
    ax[1, 0].plot(it, col[k], lw=0.8, label=k.replace("Avg_Massflow", ""))
ax[1, 0].set_ylabel("mass flow [kg/s, 360 deg]"); ax[1, 0].legend(fontsize=8)
ax[1, 1].semilogy(it, cfl, lw=0.8); ax[1, 1].set_ylabel("average CFL")
same = int(np.sum(np.cumprod((np.abs(cfl[::-1] - cfl[-1]) <= 1e-6 * cfl[-1]).astype(int))))
ax[1, 1].set_title(f"CFL {cfl[-1]:.3g}, unchanged for the last {same} iterations", fontsize=9)
for a in list(ax[0, :]) + list(ax[1, :2]):
    a.grid(alpha=0.3); a.set_xlabel("iteration")

key = lambda f: int(f.rsplit("_", 1)[-1][:-4])
snaps = sorted(glob.glob(f"{folder}/snaps/flow_*.npz"), key=key)
if len(snaps) < 2 and prev:
    snaps = sorted(glob.glob(f"{prev}/snaps/flow_*.npz"), key=key)[-1:] + snaps
if len(snaps) >= 2:
    za, zb = np.load(snaps[-2]), np.load(snaps[-1])
    ix = (zb["x_mm"] >= MAP[0]) & (zb["x_mm"] <= MAP[1]); ir = (zb["r_mm"] >= MAP[2]) & (zb["r_mm"] <= MAP[3])
    xm, rm = zb["x_mm"][ix], zb["r_mm"][ir]
    cm = lambda a: a[np.ix_(ir, ix)].astype(float)
    d = cm(zb["Mach"]) - cm(za["Mach"]) if za["Mach"].shape == zb["Mach"].shape else np.full((len(rm), len(xm)), np.nan)
    lim = max(np.nanpercentile(np.abs(d), 99.5), 1e-6) if np.isfinite(d).any() else 1.0
    ext = [xm[0], xm[-1], rm[0], rm[-1]]
    im = ax[1, 2].imshow(np.ma.masked_invalid(d), origin="lower", cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="equal",
                         extent=ext, interpolation="nearest")
    fig.colorbar(im, ax=ax[1, 2], shrink=0.8, label="dMach (clipped at p99.5)")
    ax[1, 2].contour(xm, rm, np.ma.masked_invalid(cm(zb["Mach"])), levels=[1.0, 2.0, 2.5, 3.0],
                     colors="k", linewidths=0.4, alpha=0.5)
    for q, (x, r, _) in PROBES.items():
        ax[1, 2].plot(x, r, "x", color="g", ms=7, mew=1.5, clip_on=False)
        ax[1, 2].annotate(q, (x, r), (3, 3), textcoords="offset points", color="g", fontsize=8)
    side = lambda f: ("prev. stage " if prev and f.startswith(prev) else "") + str(key(f))
    ax[1, 2].set_title(f"Mach, snapshot {side(snaps[-1])} minus {side(snaps[-2])}   max |dM| {np.nanmax(np.abs(d)):.2e}\n"
                       f"(grey: current Mach 1, 2, 2.5, 3)", fontsize=9)
    ax[1, 2].set_xlabel("x [mm] (freestream from the right)"); ax[1, 2].set_ylabel("r [mm]")
else:
    ax[1, 2].text(0.5, 0.5, "fewer than two snapshots yet", ha="center", transform=ax[1, 2].transAxes)

mass_in = col["Avg_Massflow(inlet_freestream)"][-1] + col.get("Avg_Massflow(inlet_jet)", np.zeros(1))[-1]
imb = 100 * abs(mass_in + col["Avg_Massflow(outlet)"][-1]) / mass_in
fig.suptitle(f"{([q for q in folder.split('/') if '_mesh' in q] or ['run'])[-1]} / {folder.split('/')[-1]}   iteration {int(it[-1])}   "
             f"rms[Rho] {col['rms[Rho]'][-1]:.2f}   mass imbalance {imb:.3f} %   jet {col.get('Avg_Massflow(inlet_jet)', [np.nan])[-1]:.6f} kg/s   "
             f"CD {col['CD'][-1]:.5f} (info only)", fontsize=10)
fig.tight_layout()
fig.savefig(f"{folder}/convergence.png", dpi=110)
print("saved", f"{folder}/convergence.png")
