"""python machdisk_fine.py <box.npz> <out.csv>
Mach-disk position and radius from a fine box around it (run_features.sh MODE=machdisk: x 4-5 mm, r 0-0.5 mm,
0.005 mm), because on the 0.025 mm box the disk is only ~2 cells wide (user default, 2026-11-04).
Same definitions as feature_positions.py:
- mach_disk_x: largest positive d(rho)/dx peak on the axis (parabolic sub-cell fit).
- Track the shock outward row by row (peak within +-WIN mm of the previous row).
- mach_disk_subsonic_edge_r: r where the minimum Mach within BEHIND mm downstream of the shock reaches 1 (0 if the flow
  behind the shock is supersonic on the axis already, i.e. no Mach disk).
- mach_disk_edge_r: r where the peak falls to 50 % of its axis value.
Sampling finer than the mesh only interpolates, so on coarse meshes the value is still mesh-limited.
"""
import sys
import numpy as np

WIN, BEHIND = 0.1, 0.3
z = np.load(sys.argv[1])
x, r = z["x_mm"], z["r_mm"]
h = x[1] - x[0]


def vertex(g, i):
    a, b, c = g[i - 1], g[i], g[i + 1]
    den = a - 2 * b + c
    return x[i] + (0.5 * (a - c) / den * h if den != 0 else 0.0)


def edge(rs, v, level, rising):
    hit = np.where(v >= level)[0] if rising else np.where(v < level)[0]
    hit = hit[hit > 0]
    if len(hit) == 0:
        return np.nan
    k = hit[0]
    return rs[k - 1] + (level - v[k - 1]) / (v[k] - v[k - 1]) * (rs[k] - rs[k - 1])


rows = []
for lab in [str(s) for s in z["labels"]]:
    M, G = z[f"{lab}_Mach"], np.gradient(z[f"{lab}_Density"], h, axis=1)
    i0 = int(np.nanargmax(G[0][1:-1])) + 1
    xp = x0 = vertex(G[0], i0)
    rs, gs, mb = [], [], []
    for j in range(len(r)):
        w = np.where((np.abs(x - xp) <= WIN) & np.isfinite(G[j]))[0]
        w = w[(w > 0) & (w < len(x) - 1)]
        if len(w) == 0:
            break
        i = w[np.argmax(G[j][w])]
        xp = vertex(G[j], i)
        b = (x >= x[i]) & (x <= x[i] + BEHIND)
        rs.append(r[j]); gs.append(G[j][i]); mb.append(np.nanmin(M[j][b]))
    rs, gs, mb = map(np.array, (rs, gs, mb))
    sub = 0.0 if mb[0] >= 1 else edge(rs, mb, 1.0, True)
    g50 = edge(rs, gs, 0.5 * gs[0], False)
    rows.append((lab, x0, sub, g50, mb[0]))
    print(f"{lab}: x={x0:.4f} subsonic_edge_r={sub:.4f} edge50_r={g50:.4f} M_behind_axis={mb[0]:.3f}")

with open(sys.argv[2], "w") as f:
    f.write(f"# Mach disk from machdisk_fine.py (cfd-tools); box x {x[0]}-{x[-1]} mm, r {r[0]}-{r[-1]} mm at {h:.4f} mm\n")
    f.write("model,inlet,mesh,stage,mach_disk_x_fine,mach_disk_subsonic_edge_r_fine,mach_disk_edge_r_fine,"
            "M_behind_axis\n")
    for lab, *v in rows:
        f.write(",".join(lab.split("__") + [f"{a:.4f}" for a in v]) + "\n")
print("saved", sys.argv[2])
