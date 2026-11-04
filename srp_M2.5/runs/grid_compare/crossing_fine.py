"""python crossing_fine.py <box.npz> <out.csv>
Position and normal-stem radius of the barrel-shock crossing on the axis (~4.5 mm) from a fine box around it
(run_features.sh MODE=crossing: x 4-5 mm, r 0-0.5 mm, 0.005 mm). Written as machdisk_fine.py; renamed 2026-11-04 when
the user pointed out that this feature is the shock crossing (regular reflection) and the Mach disk is the normal
shock at ~8.6 mm (feature_positions.py mach_disk_*). The 2026-11-04 run output machdisk_fine_sst.csv uses the old names.
Same definitions as feature_positions.py:
- shock_crossing_x: largest positive d(rho)/dx peak on the axis (parabolic sub-cell fit).
- Track the shock outward row by row (peak within +-WIN mm of the previous row).
- shock_crossing_subsonic_edge_r: r where the minimum Mach within BEHIND mm downstream of the shock reaches 1 (0 if the flow
  behind the shock is supersonic on the axis already, i.e. no normal stem).
- shock_crossing_edge_r: r where the peak falls to 50 % of its axis value.
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
    f.write(f"# Shock crossing from crossing_fine.py (cfd-tools); box x {x[0]}-{x[-1]} mm, r {r[0]}-{r[-1]} mm at {h:.4f} mm\n")
    f.write("model,inlet,mesh,stage,shock_crossing_x_fine,shock_crossing_subsonic_edge_r_fine,shock_crossing_edge_r_fine,"
            "M_behind_axis\n")
    for lab, *v in rows:
        f.write(",".join(lab.split("__") + [f"{a:.4f}" for a in v]) + "\n")
print("saved", sys.argv[2])
