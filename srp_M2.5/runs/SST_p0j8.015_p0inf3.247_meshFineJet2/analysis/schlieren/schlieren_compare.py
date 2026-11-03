"""python schlieren_compare.py <cfd.npz> <exp_mean.mat> <prefix> [D_mm=4.6]
Synthetic d(rho)/dx schlieren from the axisymmetric CFD, compared with the mean experimental image.
- LOS: line-of-sight integral 2*int_0^sqrt(R^2-y^2) drho/dx(x, sqrt(u^2+y^2)) du (forward Abel, u substitution so the
  kernel is not singular), out to the domain radius R. Outside the fluid drho/dx = 0.
- planar: drho/dx on the axis plane at r = |y|.
Registration: x/D = (x_cfd - x_exit)/D, y/D = +-r/D. One least-squares fit E = s*C + c per synthetic image, over pixels
where the beam does not hit the model in the CFD and the background is not dark in the experiment.
Writes <prefix>_mirror_LOS.png, <prefix>_mirror_planar.png, <prefix>_profiles.png, <prefix>_features.csv, <prefix>.npz
Env for the plot labels: CFD_LABEL (default "SST meshFineJet2"), CASE_LABEL (default "p0j 8.015 bar, p0inf 3.247 bar")."""
import os
import sys
import numpy as np
import scipy.io as sio
from scipy.interpolate import RegularGridInterpolator
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fc, fe, pre = sys.argv[1:4]
D = float(sys.argv[4]) if len(sys.argv) > 4 else 4.6
CFD_LABEL = os.environ.get("CFD_LABEL", "SST meshFineJet2")
CASE_LABEL = os.environ.get("CASE_LABEL", "p0j 8.015 bar, p0inf 3.247 bar")
c = np.load(fc)
xg, rg, xe = c["x_mm"], c["r_mm"], float(c["xe_mm"])
h = rg[1] - rg[0]
f = np.nan_to_num(c["drhodx"].astype(np.float64))         # (nr, nx), 0 outside the fluid
R = rg[-1]
e = sio.loadmat(fe)
xD, yD, E, bg = e["x"], e["y"], e["im_mean"].astype(np.float64), e["bg"].astype(np.float64)
xcol, yrow = xD[0, :], yD[:, 0]
xmm_pix = xe + xcol * D                                     # CFD x of each pixel column

# beam blocked by the model: solid (invalid, behind the exit, r < 12 mm) anywhere at r >= |y|
X, Rr = np.meshgrid(xg, rg)
solid = (~c["valid"]) & (X < xe + 0.05) & (Rr < 12.0)
rsol = np.where(solid.any(0), np.where(solid, Rr, -1).max(0), -1.0)
rsol_pix = np.interp(xmm_pix, xg, rsol)
blocked = np.abs(yrow)[:, None] * D <= rsol_pix[None, :]

# LOS projection, one row per unique |y|
los = np.empty_like(E)
cache = {}
for i, y in enumerate(np.abs(yrow) * D):
    key = round(y, 6)
    if key not in cache:
        u = np.arange(0.0, np.sqrt(max(R * R - y * y, 0.0)), h)
        fr = np.sqrt(u * u + y * y) / h
        i0 = np.minimum(fr.astype(int), len(rg) - 2)
        w = (fr - i0)[:, None]
        cache[key] = np.interp(xmm_pix, xg, 2.0 * h * ((1 - w) * f[i0] + w * f[i0 + 1]).sum(0))
    los[i] = cache[key]
fpl = c["drhodx"].astype(np.float64)
ip = RegularGridInterpolator((rg, xg), fpl, bounds_error=False, fill_value=np.nan)
P, Q = np.meshgrid(np.abs(yrow) * D, xmm_pix, indexing="ij")
plan = ip(np.stack([P, Q], -1))

expbody = bg < 30
use = ~blocked & ~expbody & np.isfinite(E)


def fit(C):
    m = use & np.isfinite(C)
    A = np.stack([C[m], np.ones(m.sum())], 1)
    (s, k), *_ = np.linalg.lstsq(A, E[m], rcond=None)
    cc = np.corrcoef(E[m], C[m])[0, 1]
    return s, k, cc


fits = {}
for name, C in (("LOS", los), ("planar", plan)):
    s, k, cc = fit(C)
    fits[name] = (s, k, cc)
    print(f"{name}: E = {s:.4g} * C + {k:.3g}   correlation {cc:.3f}  ({use.sum()} pixels)")

lim = np.nanpercentile(np.abs(E[use]), 99.5)
ext = [xcol[0], xcol[-1], yrow[0], yrow[-1]]
for name, C in (("LOS", los), ("planar", plan)):
    s, k, cc = fits[name]
    S = s * C + k
    comp = np.where(yrow[:, None] >= 0, E, S)
    comp = np.where(blocked | expbody, np.nan, comp)
    fig, ax = plt.subplots(figsize=(13, 9.5))
    cm = plt.get_cmap("gray").copy(); cm.set_bad("0.75")
    im = ax.imshow(comp, extent=ext, origin="lower" if yrow[0] < yrow[-1] else "upper", cmap=cm,
                   vmin=-lim, vmax=lim, aspect="equal", interpolation="nearest")
    ax.axhline(0, color="y", lw=0.6, ls="--")
    kw = dict(transform=ax.transAxes, ha="left", fontsize=12, color="w", fontweight="bold",
              bbox=dict(boxstyle="round", fc="k", alpha=0.6))
    ax.text(0.01, 0.97, "experiment, mean (above axis)", va="top", **kw)
    lab = f"CFD {CFD_LABEL}, line-of-sight integrated" if name == "LOS" else f"CFD {CFD_LABEL}, planar (axis plane)"
    ax.text(0.01, 0.03, f"{lab} (below axis)", va="bottom", **kw)
    ax.set_xlabel("x / D from nozzle exit (freestream from the right)"); ax.set_ylabel("y / D")
    cb = fig.colorbar(im, ax=ax, shrink=0.75); cb.set_label("schlieren signal [exp grey levels]")
    ax.set_title(f"d(rho)/dx schlieren, NPR 501 case ({CASE_LABEL}), D = {D} mm\n"
                 f"CFD scaled by one fit: E = {s:.3g} C + {k:.3g}, correlation {cc:.2f}; grey: model / masked")
    fig.tight_layout(); fig.savefig(f"{pre}_mirror_{name}.png", dpi=150); plt.close(fig)
    print("saved", f"{pre}_mirror_{name}.png")

# profiles along y = const (both halves of the experiment averaged)
lines_mm = [0.01, 1.3, 2.3, 3.0, 6.0, 9.0]


def rowavg(A, ymm):
    iu = np.argmin(np.abs(yrow * D - ymm)); il = np.argmin(np.abs(yrow * D + ymm))
    return 0.5 * (A[iu] + A[il]) if iu != il else A[iu]


SL = fits["LOS"][0] * los + fits["LOS"][1]
SP = fits["planar"][0] * plan + fits["planar"][1]
fig, axs = plt.subplots(3, 2, figsize=(15, 12), sharex=True)
for ax, ymm in zip(axs.flat, lines_mm):
    bm = rowavg((blocked | expbody).astype(float), ymm) > 0
    for A, lab, sty in ((E, "experiment", "k-"), (SL, "CFD LOS", "C0-"), (SP, "CFD planar", "C1--")):
        p = rowavg(A, ymm).copy(); p[bm] = np.nan
        ax.plot(xcol, p, sty, lw=0.9 if lab != "experiment" else 1.2, label=lab)
    ax.set_title(f"y = {ymm} mm ({ymm / D:.3f} D)"); ax.grid(alpha=0.3); ax.set_ylabel("signal [grey levels]")
axs[0, 0].legend()
for ax in axs[-1]:
    ax.set_xlabel("x / D from nozzle exit")
fig.suptitle(f"Schlieren d(rho)/dx along lines of constant y: experiment (mean of +y and -y rows) vs CFD {CFD_LABEL} (fitted scale)")
fig.tight_layout(); fig.savefig(f"{pre}_profiles.png", dpi=130); plt.close(fig)
print("saved", f"{pre}_profiles.png")

# feature positions: bow shock = most negative signal, terminal shock = most positive, in windows of x/D
rows = []
def pick(p, lo, hi, sgn):
    m = (xcol >= lo) & (xcol <= hi) & np.isfinite(p)
    j = np.where(m)[0][np.argmax(sgn * p[m])]
    return xcol[j]
for ymm, feat, lo, hi, sgn in ((0.01, "terminal shock (max)", 0.8, 2.2, 1), (0.01, "bow shock (min)", 1.6, 4.5, -1),
                               (4.6, "bow shock (min)", 1.0, 4.5, -1), (9.2, "bow shock (min)", 0.5, 4.5, -1)):
    vals = [pick(rowavg(A, ymm), lo, hi, sgn) for A in (E, SL, SP)]
    rows.append((ymm, feat, *vals))


# the experimental bow-shock band is saturated (flat bottom), so also give its upstream edge: first x, scanning from
# upstream, where the signal falls below half its minimum in the window
def edge(p, lo, hi):
    m = (xcol >= lo) & (xcol <= hi) & np.isfinite(p)
    xs, ps = xcol[m][::-1], p[m][::-1]
    return xs[np.argmax(ps < 0.5 * ps.min())]
for ymm, lo, hi in ((0.01, 1.6, 4.5), (4.6, 1.0, 4.5), (9.2, 0.5, 4.5)):
    rows.append((ymm, "bow shock upstream edge", *[edge(rowavg(A, ymm), lo, hi) for A in (E, SL, SP)]))
with open(f"{pre}_features.csv", "w") as fo:
    fo.write("y_mm,feature,exp_xD,cfdLOS_xD,cfdPlanar_xD,exp_mm,cfdLOS_mm,LOS_minus_exp_mm\n")
    print(f"{'y mm':>6} {'feature':<22} {'exp x/D':>8} {'LOS x/D':>8} {'plan x/D':>8} {'LOS-exp mm':>10}")
    for ymm, feat, a, b, p in rows:
        fo.write(f"{ymm},{feat},{a:.3f},{b:.3f},{p:.3f},{a * D:.2f},{b * D:.2f},{(b - a) * D:.2f}\n")
        print(f"{ymm:6.2f} {feat:<22} {a:8.3f} {b:8.3f} {p:8.3f} {(b - a) * D:10.2f}")
np.savez_compressed(f"{pre}.npz", xD=xcol, yD=yrow, exp=E, los=los, planar=plan, blocked=blocked, expbody=expbody,
                    fit_los=np.array(fits["LOS"]), fit_planar=np.array(fits["planar"]), xe_mm=xe, D_mm=D)
print("saved", f"{pre}.npz")
