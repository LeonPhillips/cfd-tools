"""python feature_positions.py <box.npz> <out.csv> [label ...]
Key flow-feature positions for srp_M2.5 (counterflowing jet at M 2.5) from a resample_box.py box that holds
Mach, Density and Velocity_x (VARS=Mach,Density,Velocity). One CSV row per label (default: every label in the box).

Coordinates: x along the axis (jet flows +x, freestream flows -x), r radial, both mm. Nozzle exit X_EXIT = 1.875 mm,
D = 4.6 mm. Definitions (user-approved plan, 2026-11-04):
- Shocks sit at the peak of d(rho)/dx along a line, with a parabolic sub-cell fit. Jet shocks raise rho going +x
  (positive peak); the bow shock raises rho going -x (negative peak).
- Axis: Mach disk = first and terminal shock = last positive peak between the exit and the stagnation point, counting
  peaks of at least 25 % of the strongest. Stagnation point = last u_x zero crossing (+ to -) before the bow shock.
- Bow shock per row: the "edge" is the first x from upstream (x = 20 mm) where Mach < 97 % of that row's upstream
  Mach (the old threshold method, kept as a cross-check); the shock is the d(rho)/dx minimum within -0.5/+0.25 mm of it.
- Radial extent of the Mach disk (triple-point estimate) and of the terminal-shock arc: follow the shock outward row by
  row (peak within +-0.3 mm of the previous row); the edge is the first r where the peak is below 50 % of its axis value.
- Jet sonic-line max radius: the supersonic region connected to the nozzle exit on the axis; its largest r for x below
  the stagnation point.
"""
import sys
import numpy as np

X_EXIT, D = 1.875, 4.6
z = np.load(sys.argv[1])
x, r = z["x_mm"], z["r_mm"]
h = x[1] - x[0]
labels = sys.argv[3:] or [str(s) for s in z["labels"]]


def vertex(xs, g, i):
    """Sub-cell position of an extremum of g at index i (parabola through i-1, i, i+1)."""
    if i <= 0 or i >= len(g) - 1 or not np.all(np.isfinite(g[i - 1:i + 2])):
        return xs[i]
    a, b, c = g[i - 1], g[i], g[i + 1]
    den = a - 2 * b + c
    return xs[i] + (0.5 * (a - c) / den * h if den != 0 else 0.0)


def peaks(g, lo, hi):
    """Indices of local maxima of g with lo <= x <= hi."""
    i = np.arange(1, len(g) - 1)
    m = (g[1:-1] > g[:-2]) & (g[1:-1] >= g[2:]) & np.isfinite(g[1:-1]) & (x[i] >= lo) & (x[i] <= hi)
    return i[m]


def bow(Mrow, grow):
    """Bow shock on one row: (edge x by Mach threshold, d(rho)/dx-minimum x, upstream Mach)."""
    ok = np.where(np.isfinite(Mrow))[0]
    Mup = Mrow[ok[-1]]
    below = np.where(np.isfinite(Mrow) & (Mrow < 0.97 * Mup))[0]
    if len(below) == 0:
        return np.nan, np.nan, Mup
    k = below[-1]   # first point below the threshold seen from upstream
    xe = x[k] + (x[k + 1] - x[k]) * (0.97 * Mup - Mrow[k]) / (Mrow[k + 1] - Mrow[k])
    w = np.where((x >= xe - 0.5) & (x <= xe + 0.25) & np.isfinite(grow))[0]
    i = w[np.argmin(grow[w])]
    return xe, vertex(x, -grow, i), Mup


def track(G, x0, s0, sign, win=0.3):
    """Follow a shock outward from the axis. Returns (edge r, x at the edge, r and x of the track)."""
    rs, xs, xp = [0.0], [x0], x0
    for j in range(1, len(r)):
        g = sign * G[j]
        w = np.where((np.abs(x - xp) <= win) & np.isfinite(g))[0]
        if len(w) == 0:
            break
        i = w[np.argmax(g[w])]
        if g[i] < 0.5 * s0:
            return r[j], xs[-1], rs, xs
        xp = vertex(x, g, i); rs.append(r[j]); xs.append(xp)
    return np.nan, np.nan, rs, xs   # never dropped below 50 % inside the box


def connected(mask, seed):
    """Region of mask 4-connected to seed (iterative dilation)."""
    reg = np.zeros_like(mask); reg[seed] = mask[seed]
    while True:
        n = reg.copy()
        n[1:] |= reg[:-1]; n[:-1] |= reg[1:]; n[:, 1:] |= reg[:, :-1]; n[:, :-1] |= reg[:, 1:]
        n &= mask
        if (n == reg).all():
            return reg
        reg = n


RB = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]   # bow shock shape, r in mm
cols = ["label", "M_upstream_axis",
        "mach_disk_x", "mach_disk_M_before", "mach_disk_edge_r", "mach_disk_edge_x",
        "terminal_shock_x", "terminal_shock_M_before", "terminal_arc_edge_r", "terminal_arc_edge_x",
        "stagnation_x",
        "bow_x_axis", "bow_edge_x_axis", "bow_x_1D", "bow_edge_x_1D", "bow_x_2D", "bow_edge_x_2D",
        "standoff_bow_minus_exit_axis", "shock_layer_bow_minus_stag_axis",
        "sonic_rmax", "sonic_rmax_x"] + [f"bow_x_r{v}" for v in RB] + ["n_jet_peaks", "flags"]
rows = []
for lab in labels:
    M, rho, u = z[f"{lab}_Mach"], z[f"{lab}_Density"], z[f"{lab}_Velocity_x"]
    G = np.gradient(rho, h, axis=1)   # d(rho)/dx per mm, rows r
    flags = []
    out = {"label": lab}
    # axis
    xb_edge, xb, Mup = bow(M[0], G[0])
    out.update(M_upstream_axis=Mup, bow_x_axis=xb, bow_edge_x_axis=xb_edge)
    c = np.where((x > X_EXIT) & (x < xb) & np.isfinite(u[0]) & (u[0] > 0) & (np.roll(u[0], -1) <= 0))[0]
    xs = x[c[-1]] + h * u[0][c[-1]] / (u[0][c[-1]] - u[0][c[-1] + 1]) if len(c) else np.nan
    if len(c) > 1:
        flags.append(f"{len(c)} u=0 crossings on axis")
    out["stagnation_x"] = xs
    p = peaks(G[0], X_EXIT + 0.1, xs)
    p = p[G[0][p] >= 0.25 * G[0][p].max()] if len(p) else p
    out["n_jet_peaks"] = len(p)
    if len(p) >= 2:
        for name, i, lo in (("mach_disk", p[0], X_EXIT), ("terminal_shock", p[-1], x[p[-2]])):
            out[f"{name}_x"] = xi = vertex(x, G[0], i)
            w = (x >= lo) & (x < xi)
            out[f"{name}_M_before"] = np.nanmax(M[0][w])
        er, ex, _, _ = track(G, out["mach_disk_x"], G[0][p[0]], +1)
        out.update(mach_disk_edge_r=er, mach_disk_edge_x=ex)
        er, ex, _, _ = track(G, out["terminal_shock_x"], G[0][p[-1]], +1)
        out.update(terminal_arc_edge_r=er, terminal_arc_edge_x=ex)
    else:
        flags.append(f"only {len(p)} jet shock peak(s) on axis")
    out["standoff_bow_minus_exit_axis"] = xb - X_EXIT
    out["shock_layer_bow_minus_stag_axis"] = xb - xs
    # bow shock off the axis
    for key, rv in (("1D", D), ("2D", 2 * D)):
        j = int(round(rv / (r[1] - r[0])))
        e, b, _ = bow(M[j], G[j]); out[f"bow_edge_x_{key}"] = e; out[f"bow_x_{key}"] = b
    for rv in RB:
        j = int(round(rv / (r[1] - r[0])))
        out[f"bow_x_r{rv}"] = bow(M[j], G[j])[1]
    # jet sonic line
    sup = np.isfinite(M) & (M > 1.0)
    i0 = int(np.searchsorted(x, X_EXIT + 0.05))
    reg = connected(sup, (0, i0))
    if (reg & (x[None, :] > xb)).any():
        flags.append("jet supersonic region reaches the freestream")
    js, jx = np.where(reg & (x[None, :] < xs))
    if len(js):
        k = np.argmax(js); out.update(sonic_rmax=r[js[k]], sonic_rmax_x=x[jx[k]])
    out["flags"] = "; ".join(flags)
    rows.append(out)
    print(lab, " ".join(f"{k}={out.get(k, np.nan):.3f}" for k in
                        ("mach_disk_x", "mach_disk_edge_r", "terminal_shock_x", "terminal_arc_edge_r", "stagnation_x",
                         "bow_x_axis", "bow_x_1D", "bow_x_2D", "sonic_rmax")), out["flags"])

with open(sys.argv[2], "w") as f:
    f.write("# srp_M2.5 key-feature positions from feature_positions.py (cfd-tools); box 0.025 mm, x/r in mm, "
            f"nozzle exit x = {X_EXIT} mm, D = {D} mm; definitions in the script header\n")
    # labels are model__inlet__mesh__stage (run_features.sh); written as four columns
    f.write(",".join(["model", "inlet", "mesh", "stage"] + cols[1:]) + "\n")
    for o in rows:
        parts = (o["label"].split("__") + ["", "", "", ""])[:4]
        f.write(",".join(parts + [o[k] if k == "flags" else
                         (str(o[k]) if k == "n_jet_peaks" else f"{o.get(k, np.nan):.4f}") for k in cols[1:]]) + "\n")
print("saved", sys.argv[2])
