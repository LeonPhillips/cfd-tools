"""python diff_contours.py <box.npz> <outdir>: contour plots of very fine - fine over the refinement box, one figure for
Mach and one for density (% of very fine), plus region statistics and a 0.02 mm copy for the interactive html."""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
z = np.load(sys.argv[1]); out = sys.argv[2]
x, r = z["x_mm"], z["r_mm"]; X, R = np.meshgrid(x, r)
mf, mv = z["fine_Mach"].astype(float), z["veryFine_Mach"].astype(float)
df, dv = z["fine_Density"].astype(float), z["veryFine_Density"].astype(float)
dM = mv - mf; dRho = 100 * (dv - df) / dv
# bow shock per row: first x from upstream (x = 20 mm) where very-fine Mach drops below 2.425
bow = np.full(len(r), np.nan)
for j in range(len(r)):
    below = np.flatnonzero(~(mv[j] >= 2.425))
    if below.size: bow[j] = x[below[-1]]
B = bow[:, None]
regions = {
    "jet core x 2-8, r<2": (X >= 2) & (X < 8) & (R < 2),
    "jet shear layer x 2-8, r 2-4": (X >= 2) & (X < 8) & (R >= 2) & (R < 4),
    "subsonic cap x>8.5 to bow-0.2": (X > 8.5) & (X < B - 0.2),
    "bow shock layer |x-bow|<=0.2": np.abs(X - B) <= 0.2,
    "upstream x>bow+0.2": X > B + 0.2,
    "whole box": np.ones_like(X, bool),
}
print("very fine - fine; bow shock found at x =", f"{np.nanmin(bow):.2f}..{np.nanmax(bow):.2f} mm")
print(f"{'region':32s}{'points':>9s}{'mean|dM|':>10s}{'p99|dM|':>9s}{'mean dM':>10s}{'mean|drho|%':>12s}{'p99|drho|%':>11s}{'mean drho%':>11s}")
for name, m in regions.items():
    m = m & np.isfinite(dM) & np.isfinite(dRho)
    a, b = dM[m], dRho[m]
    print(f"{name:32s}{m.sum():9d}{np.abs(a).mean():10.4f}{np.percentile(np.abs(a), 99):9.4f}{a.mean():10.4f}"
          f"{np.abs(b).mean():12.3f}{np.percentile(np.abs(b), 99):11.3f}{b.mean():11.3f}")
ext = [x[0], x[-1], r[0], r[-1]]
LINES = [0.01, 1.3, 2.3, 3.0, 6.0, 9.0]
ZOOMS = [((2, 9), (0, 4), "jet core, Mach disk and shear layer"), ((13, 17), (0, 9.5), "bow shock layer")]


def figure(field, diff, fname, flabel, dlabel, clip, title, fvmin=None, fvmax=None):
    fig, ax = plt.subplots(2, 2, figsize=(18, 11), constrained_layout=True)
    panels = [(ax[0, 0], field, "viridis", (fvmin, fvmax), f"very fine: {flabel}", None),
              (ax[0, 1], diff, "RdBu_r", (-clip, clip), f"very fine - fine: {dlabel} (clipped at +/-{clip:g})", None)]
    panels += [(a, diff, "RdBu_r", (-clip, clip), f"very fine - fine, zoom: {t}", (xl, yl)) for a, (xl, yl, t) in zip(ax[1], ZOOMS)]
    for a, data, cmap, (lo, hi), t, zoom in panels:
        im = a.imshow(np.ma.masked_invalid(data), origin="lower", extent=ext, cmap=cmap, vmin=lo, vmax=hi,
                      interpolation="nearest", aspect="equal")
        for rl in LINES:
            a.axhline(rl, color="k", lw=0.5, ls="--", alpha=0.35)
        if zoom:
            a.set_xlim(zoom[0]); a.set_ylim(zoom[1])
        a.set_xlabel("x [mm] (freestream from the right)"); a.set_ylabel("r [mm]"); a.set_title(t, fontsize=10)
        fig.colorbar(im, ax=a, shrink=0.85)
    fig.suptitle(title + "  (dashed: standard line positions r = 0.01, 1.3, 2.3, 3, 6, 9 mm)")
    fig.savefig(f"{out}/{fname}", dpi=130); print("saved", f"{out}/{fname}")


figure(mv, dM, "diff_contours_Mach.png", "Mach", "dMach", 0.05,
       "Mach: very fine (0.0125 mm) minus fine (0.0177 mm), refinement box x 0.5-20 mm, r 0-9.5 mm", 0, 3.4)
figure(dv, dRho, "diff_contours_density.png", "density [kg/m^3]", "drho [% of very fine]", 3.0,
       "Density: very fine (0.0125 mm) minus fine (0.0177 mm), % of very-fine density, refinement box")
np.savez_compressed(f"{out}/diff_box_0.02mm.npz", x_mm=x[::2], r_mm=r[::2], veryFine_Mach=mv[::2, ::2].astype(np.float32),
                    dMach=dM[::2, ::2].astype(np.float32), veryFine_Density=dv[::2, ::2].astype(np.float32),
                    dRho_pct=dRho[::2, ::2].astype(np.float32))
print("saved", f"{out}/diff_box_0.02mm.npz")
