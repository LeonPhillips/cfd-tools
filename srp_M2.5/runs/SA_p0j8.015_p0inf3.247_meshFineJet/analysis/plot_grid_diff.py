"""python plot_grid_diff.py <grid_box.npz> <out.png>: very fine minus fineJet, Mach and density (% of very fine), over the
whole box and zoomed on the jet box (x 1.5-8.5, r 0-3.6 mm). Colour clipped at the 99.5th percentile of |difference|;
grey: very-fine Mach 1, 2, 2.5, 3; dotted: the line-profile radii."""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
z = np.load(sys.argv[1]); x, r = z["x_mm"], z["r_mm"]
vf, fj = "veryFine", "fineJet"
dM = z[f"{vf}_Mach"].astype(float) - z[f"{fj}_Mach"].astype(float)
dR = 100 * (z[f"{vf}_Density"].astype(float) - z[f"{fj}_Density"].astype(float)) / z[f"{vf}_Density"]
M = np.ma.masked_invalid(z[f"{vf}_Mach"])
LINES = [0.01, 1.3, 2.3, 3, 6, 9]
fig, ax = plt.subplots(2, 2, figsize=(17, 10), gridspec_kw={"width_ratios": [1.6, 1]})
for i, (d, name, unit) in enumerate([(dM, "Mach", "dMach"), (dR, "density", "drho / rho_veryFine [%]")]):
    for j, (xa, xb, ra, rb) in enumerate([(x[0], x[-1], r[0], r[-1]), (1.5, 8.5, 0.0, 3.6)]):
        a = ax[i, j]
        sx, sr = (x >= xa) & (x <= xb), (r >= ra) & (r <= rb)
        dd = d[np.ix_(sr, sx)]
        lim = max(np.nanpercentile(np.abs(dd), 99.5), 1e-9)
        im = a.imshow(np.ma.masked_invalid(dd), origin="lower", cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="equal",
                      extent=[x[sx][0], x[sx][-1], r[sr][0], r[sr][-1]], interpolation="nearest")
        a.contour(x[sx], r[sr], M[np.ix_(sr, sx)], levels=[1.0, 2.0, 2.5, 3.0], colors="k", linewidths=0.4, alpha=0.5)
        for rl in LINES:
            if ra <= rl <= rb:
                a.axhline(rl, color="g", ls=":", lw=0.8)
        fig.colorbar(im, ax=a, shrink=0.85, label=f"{unit} (clipped at p99.5)")
        a.set_title(f"{name}: very fine minus fineJet   mean |d| {np.nanmean(np.abs(dd)):.2e}, p99 "
                    f"{np.nanpercentile(np.abs(dd), 99):.2e}, max {np.nanmax(np.abs(dd)):.2e}", fontsize=9)
        a.set_xlabel("x [mm] (freestream from the right)"); a.set_ylabel("r [mm]")
fig.suptitle("SA_p0j8.015_p0inf3.247: very fine v3 minus fineJet j3 (end of each run), resampled at 0.01 mm\n"
             "grey: very-fine Mach 1, 2, 2.5, 3   green dotted: line-profile radii", fontsize=11)
fig.tight_layout()
fig.savefig(sys.argv[2], dpi=110)
print("saved", sys.argv[2])
