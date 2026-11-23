"""python seed_jet.py <base restart.dat> <out restart.dat>
Seeded starting field for the coarse free jet (user, 2026-11-23, option 1). Started from uniform ambient, the steady
solver never carries the jet to the outlet: the jet front acts like a wall, and the jet turns outward into a hollow
cone around reversed flow on the axis. Here the starting field already has a decaying round jet along the axis out
to the outlet, so there is no front to push through. Only the starting field changes, not the physics.

- Kept from the base solution (a coarse-mesh restart with the nozzle running): the plenum and nozzle interior,
  -15 mm < x < 0 with r < 9 mm (inside the body there is only plenum/nozzle fluid), and 0 <= x < XE with r < 2.45 mm.
- Everywhere else: ambient at pa, T0 = 288.15 K with the Mach 0.01 coflow (3.40 m/s), plus for x' = x - XE > 0 a
  self-similar round jet:
    Uc(x') = Uj * min(1, LC / x'),  r_half(x') = DJ/2 + 0.094 x',  U = Ucf + (Uc - Ucf) exp(-ln2 (r/r_half)^2), Ur = 0
  Uj 490 m/s and DJ 4.92 mm are the fully expanded jet at NPR 6.51 (Mj 1.88, Tj 168.7 K, mdot 0.0195 kg/s).
  LC = 60 mm (between 6 momentum diameters and the ~70 mm core of the first-order runs).
  T from constant total temperature (T = T0 - U^2 / 2cp), p = pa, rho = p / (R T).
  k = max(k_inf, 0.05 (Uc - Ucf)^2 min(1, x'/LC) g),  nu_t = 0.028 (Uc - Ucf) r_half g',  omega = max(omega_inf, k / nu_t)
  with g, g' Gaussian in r/r_half (widths 1 and 1.4), k_inf, omega_inf the freestream values (TI 1 %, mu_t/mu 1).
The output has the base file's header and point order; x and y are copied unchanged.
"""
import sys
import numpy as np

GAM, R = 1.4, 287.058
CP, CV = GAM * R / (GAM - 1), R / (GAM - 1)
PA, T0 = 101325.0, 288.15
XE, RLIP = 1.876e-3, 2.45e-3
UJ, DJ, LC, SPREAD = 490.0, 4.92e-3, 60e-3, 0.094
LN2 = np.log(2.0)

base, out = sys.argv[1], sys.argv[2]
with open(base, "rb") as f:
    head = f.read(20)
    h = np.frombuffer(head, dtype=np.int32)
    assert h[0] == 535532, "not an SU2 binary restart"
    nv, n = int(h[1]), int(h[2])
    names_raw = f.read(33 * nv)
    names = [names_raw[33 * i:33 * (i + 1)].split(b"\0")[0].decode() for i in range(nv)]
    d = np.frombuffer(f.read(8 * nv * n), dtype=np.float64).reshape(n, nv).copy()
    assert f.read() == b"", "unexpected trailing data"
want = ["x", "y", "Density", "Momentum_x", "Momentum_y", "Energy", "Turb_Kin_Energy", "Omega"]
assert names[:8] == want, f"unexpected fields {names}"
x, r = d[:, 0], d[:, 1]

# freestream
Tinf = T0 / (1 + 0.5 * (GAM - 1) * 0.01 ** 2)
Ucf = 0.01 * np.sqrt(GAM * R * Tinf)
mu = lambda T: 1.716e-5 * (T / 273.15) ** 1.5 * (273.15 + 110.4) / (T + 110.4)
rho_inf = PA / (R * Tinf)
k_inf = 1.5 * (0.01 * Ucf) ** 2
om_inf = rho_inf * k_inf / mu(Tinf)

keep = ((x > -15e-3) & (x < 0.0) & (r < 9e-3)) | ((x >= 0.0) & (x < XE) & (r < RLIP))
xp = x - XE
jet = xp > 0
xj = np.where(jet, xp, 0.0)
Uc = UJ * np.minimum(1.0, LC / np.maximum(xj, 1e-12))
rh = DJ / 2 + SPREAD * xj
eta = r / rh
dU = np.where(jet, Uc - Ucf, 0.0)
U = Ucf + dU * np.exp(-LN2 * eta ** 2)
T = T0 - U ** 2 / (2 * CP)
rho = PA / (R * T)
k = np.maximum(k_inf, 0.05 * dU ** 2 * np.minimum(1.0, xj / LC) * np.exp(-LN2 * eta ** 2))
nut = np.maximum(mu(T) / rho, 0.028 * dU * rh * np.exp(-LN2 * (eta / 1.4) ** 2))
om = np.maximum(om_inf, k / nut)

s = ~keep
d[s, 2] = rho[s]
d[s, 3] = rho[s] * U[s]
d[s, 4] = 0.0
d[s, 5] = rho[s] * (CV * T[s] + 0.5 * U[s] ** 2)
d[s, 6] = k[s]
d[s, 7] = om[s]
with open(out, "wb") as f:
    f.write(head)
    f.write(names_raw)
    f.write(d.tobytes())

print(f"base {base}: {n} points, {nv} fields")
print(f"kept from base (plenum and nozzle): {keep.sum()} points; seeded: {s.sum()} points")
print(f"freestream: T {Tinf:.2f} K, U {Ucf:.3f} m/s, rho {rho_inf:.4f}, k {k_inf:.3e}, omega {om_inf:.1f}")
for xm in (5, 20, 60, 100, 160, 220, 238):
    sel = s & (np.abs(x - xm * 1e-3) < 1.0e-3) & (r < 1.0e-3)
    if sel.any():
        i = np.where(sel)[0][np.argmin(r[sel])]
        print(f"  x {xm:3d} mm, r {r[i]*1e3:5.2f} mm: Ux {d[i,3]/d[i,2]:7.1f} m/s, T {T[i]:6.1f} K, "
              f"k {k[i]:.3e}, nu_t/nu {k[i]/om[i]/(mu(T[i])/rho[i]):.0f}, r_half {rh[i]*1e3:.2f} mm")
print(f"wrote {out}")
