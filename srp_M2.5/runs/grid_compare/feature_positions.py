"""python feature_positions.py <box.npz> <out.csv> [label ...]
Key flow-feature positions for srp_M2.5 (counterflowing jet at M 2.5) from a resample_box.py box that holds
Mach, Density and Velocity_x (VARS=Mach,Density,Velocity). One CSV row per label (default: every label in the box).

Coordinates: x along the axis (jet flows +x, freestream flows -x), r radial, both mm. Nozzle exit X_EXIT = 1.875 mm,
D = 4.6 mm. Definitions (user-approved plan, 2026-11-04):
- Shocks sit at the peak of d(rho)/dx along a line, with a parabolic sub-cell fit. Jet shocks raise rho going +x
  (positive peak); the bow shock raises rho going -x (negative peak).
- Axis: Mach disk = first and terminal shock = last positive peak between the exit and the stagnation point, counting
  peaks of at least 25 % of the strongest. If only the Mach disk passes (SA: smeared terminal compression), the
  strongest peak more than 0.5 mm after it is used and flagged. Stagnation point = first u_x zero crossing (+ to -)
  after the exit (coarse meshes have further crossings from a recirculation pocket; all are in axis_u0_x).
- *_subsonic_edge_r (second criterion, added 2026-11-04 after the test): along the same track, the r where the minimum
  Mach within 0.3 mm behind the shock reaches 1, i.e. the end of the normal (subsonic-behind) part. The Mach disk is
  only ~2 box cells (0.05 mm) wide on this criterion, so its value is resolution-limited.
- Bow shock per row: the "edge" is the first x from upstream (x = 20 mm) where Mach < 97 % of that row's upstream
  Mach (the old threshold method, kept as a cross-check); the shock is the d(rho)/dx minimum within -0.5/+0.25 mm of it.
- Radial extent of the Mach disk (triple-point estimate) and of the terminal-shock arc: follow the shock outward row by
  row (peak within +-0.3 mm of the previous row); *_edge_r is the r where the peak falls to 50 % of its axis value
  (linear between rows). Fragile: the terminal-shock strength off the axis hovers near 55-60 % of the axis value.
- Jet sonic-line max radius: per x column from the exit to the terminal shock, the outer M = 1 edge of the first
  supersonic band out from the axis; the largest r and its x.
- axis_u0_x: every u_x sign change on the axis between the exit and the bow shock; (+/-) = +x flow ends there.
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


def track(G, M, x0, sign, win=0.3, behind=0.3):
    """Follow a jet shock outward from the axis, row by row (peak within +-win of the previous row).
    Returns per row: r, x, peak strength, min Mach within `behind` mm downstream (+x) of the shock."""
    rs, xs, gs, mb, xp = [], [], [], [], x0
    for j in range(len(r)):
        g = sign * G[j]
        w = np.where((np.abs(x - xp) <= win) & np.isfinite(g))[0]
        if len(w) == 0:
            break
        i = w[np.argmax(g[w])]
        xp = vertex(x, g, i)
        b = (x >= x[i]) & (x <= x[i] + behind) & np.isfinite(M[j])
        rs.append(r[j]); xs.append(xp); gs.append(g[i]); mb.append(M[j][b].min() if b.any() else np.nan)
    return np.array(rs), np.array(xs), np.array(gs), np.array(mb)


def edge(rs, xs, v, level, rising):
    """First r where v crosses `level` (rising: from below to >= level; else falling below it), linear between rows.
    Returns (r, x of the last row before it), or NaN if it never crosses inside the box."""
    hit = np.where(v >= level)[0] if rising else np.where(v < level)[0]
    hit = hit[hit > 0]
    if len(hit) == 0:
        return np.nan, np.nan
    k = hit[0]
    f = (level - v[k - 1]) / (v[k] - v[k - 1]) if v[k] != v[k - 1] else 0.0
    return rs[k - 1] + f * (rs[k] - rs[k - 1]), xs[k - 1]


RB = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]   # bow shock shape, r in mm
cols = ["label", "M_upstream_axis",
        "mach_disk_x", "mach_disk_M_before", "mach_disk_edge_r", "mach_disk_edge_x",
        "mach_disk_subsonic_edge_r", "mach_disk_subsonic_edge_x",
        "terminal_shock_x", "terminal_shock_M_before", "terminal_shock_edge_r", "terminal_shock_edge_x",
        "terminal_shock_subsonic_edge_r", "terminal_shock_subsonic_edge_x",
        "stagnation_x",
        "bow_x_axis", "bow_edge_x_axis", "bow_x_1D", "bow_edge_x_1D", "bow_x_2D", "bow_edge_x_2D",
        "standoff_bow_minus_exit_axis", "shock_layer_bow_minus_stag_axis",
        "sonic_rmax", "sonic_rmax_x"] + [f"bow_x_r{v}" for v in RB] + ["n_jet_peaks", "axis_u0_x", "flags"]
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
    xs = x[c[0]] + h * u[0][c[0]] / (u[0][c[0]] - u[0][c[0] + 1]) if len(c) else np.nan   # first +x-flow end
    out["stagnation_x"] = xs
    pa = peaks(G[0], X_EXIT + 0.1, xs)
    p = pa[G[0][pa] >= 0.25 * G[0][pa].max()] if len(pa) else pa
    out["n_jet_peaks"] = len(p)
    if len(p) == 1:   # smeared terminal compression (SA): strongest peak after the Mach disk, flagged
        q = pa[x[pa] > x[p[0]] + 0.5]
        if len(q):
            k = q[np.argmax(G[0][q])]
            p = np.array([p[0], k])
            flags.append(f"terminal compression weak/smeared ({100 * G[0][k] / G[0][p[0]]:.0f} % of Mach disk peak)")
    if len(p) >= 2:
        for name, i, lo in (("mach_disk", p[0], X_EXIT), ("terminal_shock", p[-1], x[p[-2]])):
            out[f"{name}_x"] = xi = vertex(x, G[0], i)
            w = (x >= lo) & (x < xi)
            out[f"{name}_M_before"] = np.nanmax(M[0][w])
            rs_, xs_, gs_, mb_ = track(G, M, xi, +1)
            out[f"{name}_edge_r"], out[f"{name}_edge_x"] = edge(rs_, xs_, gs_, 0.5 * gs_[0], rising=False)
            out[f"{name}_subsonic_edge_r"], out[f"{name}_subsonic_edge_x"] = (
                edge(rs_, xs_, mb_, 1.0, rising=True) if mb_[0] < 1 else (0.0, xi))
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
    # jet sonic line: per column, outer edge of the first supersonic band out from the axis, exit to terminal shock
    # (the jet core joins an outer supersonic region at large r, so a connected-region search leaks; 2026-11-04 test)
    best = (np.nan, np.nan)
    for c in np.where((x >= X_EXIT + 0.05) & (x <= out.get("terminal_shock_x", xs)))[0]:
        s = np.isfinite(M[:, c]) & (M[:, c] > 1.0)
        on = np.where(s)[0]
        if len(on) == 0:
            continue
        off = np.where(~s[on[0]:])[0]
        if len(off) == 0:
            flags.append(f"sonic band reaches the box top at x={x[c]:.2f}"); break
        j = on[0] + off[0]   # first subsonic row
        rj = r[j - 1] + (r[j] - r[j - 1]) * (M[j - 1, c] - 1.0) / (M[j - 1, c] - M[j, c])
        if not rj <= best[0]:
            best = (rj, x[c])
    out.update(sonic_rmax=best[0], sonic_rmax_x=best[1])
    # every u_x sign change on the axis between the exit and the bow shock
    w = np.where((x > X_EXIT) & (x < xb))[0]
    out["axis_u0_x"] = " ".join(f"{x[i] + h * u[0][i] / (u[0][i] - u[0][i + 1]):.3f}({'+/-' if u[0][i] > 0 else '-/+'})"
                                for i in w[:-1] if np.isfinite(u[0][i]) and np.isfinite(u[0][i + 1])
                                and np.sign(u[0][i]) != np.sign(u[0][i + 1]))
    out["flags"] = "; ".join(flags)
    rows.append(out)
    print(lab, " ".join(f"{k}={out.get(k, np.nan):.3f}" for k in
                        ("mach_disk_x", "mach_disk_subsonic_edge_r", "terminal_shock_x", "terminal_shock_edge_r",
                         "terminal_shock_subsonic_edge_r", "stagnation_x",
                         "bow_x_axis", "bow_x_1D", "bow_x_2D", "sonic_rmax")), out["flags"])

with open(sys.argv[2], "w") as f:
    f.write("# srp_M2.5 key-feature positions from feature_positions.py (cfd-tools); box 0.025 mm, x/r in mm, "
            f"nozzle exit x = {X_EXIT} mm, D = {D} mm; definitions in the script header\n")
    # labels are model__inlet__mesh__stage (run_features.sh); written as four columns
    f.write(",".join(["model", "inlet", "mesh", "stage"] + cols[1:]) + "\n")
    for o in rows:
        parts = (o["label"].split("__") + ["", "", "", ""])[:4]
        f.write(",".join(parts + [o[k] if k in ("flags", "axis_u0_x") else
                         (str(o[k]) if k == "n_jet_peaks" else f"{o.get(k, np.nan):.4f}") for k in cols[1:]]) + "\n")
print("saved", sys.argv[2])
