# lsa: local spatial linear stability analysis of the srp_M2.5 jet

`lsa_local.m` is one MATLAB file, commented line by line. At each axial station it treats the RANS mean flow as
locally parallel and finds the complex axial wavenumbers alpha of waves exp(i(alpha x + m theta - omega t)) for a
given real frequency. The method follows Tran, Nogueira, Weightman & Edgington-Mitchell 2026 (JFM 1046 A14, §2.3 and
App. D): inviscid linearised Euler equations in specific volume, mean radial velocity dropped, Chebyshev collocation
with mapping, `eig` on L q = alpha R q.

Work folder on M3: `/scratch/qk75/LeonPhillips/PhD/LSA/` (a copy of `lsa_local.m` stamped with the commit hash, plus
one output folder per mode).

## Confirmed assumptions (user, 2026-11-10)

1. Base flow: `SST_p0j6.6_p0inf3.247_meshFine/2_jet_staged/t3_cap3` (the GSA base flow), read from
   `/scratch/qk75/LeonPhillips/PhD/GSA/srp_p0j6.6_meshFine_t3/baseflow/su2L_flow.{bin,json}`.
2. Inviscid; no eddy viscosity; mean u_r (and u_theta) dropped, as in Tran 2026.
3. Scales: D = 18 mm, a_inf = 226.865 m/s, rho_inf = 0.51693 kg/m3, p by rho_inf a_inf^2 (free-stream p = 1/gamma).
   omega (code units) = 2 pi f D / a_inf = 2 pi St_D M_inf, with St_D = f D / U_inf, U_inf = 567.162 m/s,
   so f [Hz] = 31,509 St_D.
4. Stations x = 1.9 to 8.0 mm every 0.1 mm. Nozzle exit x_e = 1.876 mm, inner lip radius 2.389 mm (from the
   exported wall; throat r = 2.000 mm at x = 0). The axis is subsonic by x = 8.5 mm (Mach disk). Stations inside the
   nozzle are not included.
5. r_max just inside the bow shock at each station (bow-shock radius - 0.5 mm - 3 smoothing lengths), Dirichlet there.
6. Base-flow profile is a switch, `PROFILE` (user, 2026-11-10: Tran's tanh fit, left as a switch):
   - `'tanh'` (default): Tran eqs. 3.1-3.3, two tanh layers per variable (U, rho, p), fitted by variable projection
     (linear least squares for the amplitudes, `fminsearch` for R and thickness). Differences from Tran: the outer
     constant is fitted (the outer gas is the moving shock layer, not still air), and no layer may be thinner than
     0.05 mm half-thickness (one RANS cell; without this the fit snapped onto the lip and barrel-shock jumps).
   - `'smooth'`: Whittaker (penalised least-squares) smoother, cut-off length 0.05 mm.
7. m = 0 and 1; blind scan St_D 0.05-2 (experimental tones withheld).
8. N = 250 Chebyshev points, clustered at the jet shear layer (r1) and the outer shear layer (r2). An eigenvalue
   counts only if the N = 300 solve reproduces it to 1e-3 (relative).
9. Spectra of KH and upstream (guided-jet-like) waves only; no closure (STIM) frequency prediction.

**Caveat found on the profiles:** the mean radial velocity is large near the jet boundary (|V|/U_max about 0.6-0.7),
so dropping it is a strong approximation in this flow. The flow also has two shear layers: the jet boundary
(r about 2.3-4.7 mm) and an outer one (r about 7-9 mm) between the reverse-flow pocket and the deflected free stream.

## Running

On M3, on a compute node (your desktop job or an sbatch job), never on the login node:

```bash
module load matlab/r2023b
cd /scratch/qk75/LeonPhillips/PhD/LSA
matlab -batch "MODE='selfcheck'; lsa_local"
matlab -batch "MODE='single'; lsa_local"
matlab -batch "MODE='sensitivity'; lsa_local"
matlab -batch "MODE='scan'; NWORKERS=20; lsa_local"
```

| MODE | What it does | Output (in `LSA/<MODE>/`) |
|---|---|---|
| `selfcheck` | Thin tanh jet (M_j 1.5, cold and nu_j = 1.5) against the exact compressible vortex sheet, R/theta = 20, 40, 80 | `selfcheck.csv`, PASSED/FAILED line in the log |
| `single` | One station and St_D (`SINGLE_X`, `SINGLE_ST`; default 3.0 mm, 1.271), m = 0 and 1: full spectrum, labels, eigenfunctions | `x_<x>_St_<St>_lsa_single.mat`, `x_<x>_St_<St>_lsa_single_spectrum.png` |
| `sensitivity` | x = 2.5, 4, 6, 7.5 mm, St_D 0.25, 1, 2: N 150/250/350, r_max 2.5 mm smaller, smoothing 0.025/0.1 mm | `lsa_sensitivity.csv` |
| `scan` | All stations x St x m | `lsa_scan.mat`, `lsa_scan_m0.png`, `lsa_scan_m1.png`, `lsa_scan_radii.png` |

Every mode also writes `lsa_profiles.png` (raw and smoothed profiles) and `lsa_local.log`.

## Mode labels (heuristics, not a Briggs-Bers analysis)

- `continuum`: alpha = omega/U(r_j) at a collocation point, or a neutral non-dispersive mode (group velocity = phase
  speed) moving at a speed the mean flow has. These are disturbances carried by the mean flow (the inviscid
  continuous spectrum), discretised as dense clouds of nearly neutral eigenvalues. Never counted as waves.
- `KH_jet`: clearly growing (alpha_i < -1e-3 |alpha|), group velocity > 0, |u| peak within 0.6 mm of the jet shear
  layer r1. (|u| rather than |p|: the m = 0 KH pressure can peak inside the jet.)
- `KH_outer`: the same near the outer shear layer r2.
- `upstream_jet`: |alpha_i| < 0.02 |alpha_r|, group velocity < 0, more than half of the |p|^2 r energy inside r1.
- `other`: everything else (including evanescent duct-like pairs at large |alpha_i|).

**Convergence test.** An eigenvalue counts only if the N = 300 solve has one within 1e-3 (relative) AND that match is
10 times closer than the eigenvalue's nearest neighbour in its own spectrum. The second condition stops a dense
continuum cloud from passing (found on the first x = 3 mm run, where most "converged" modes were continuum).

**Self-check result (local, 2026-11-10).** St_j 0.4 and 0.6 pass in all four cases (cold and nu_j 1.5, m = 0
and 1): the error to the exact vortex sheet falls at first order (err(40)/err(80) = 2.11-2.14), and the two-point
extrapolation lands 1.8-3.0 % (St_j 0.4) and 4.5-7.4 % (St_j 0.6) from the root, the O((alpha theta)^2) residual.
**St_j 0.2 fails:** that KH wave has a phase speed above the ambient sound speed and grows slowly, so it radiates
and its eigenfunction hardly decays in r; the Dirichlet condition at r_max reflects the radiation and the mode is
lost among box modes. The same can affect weakly growing KH modes of the srp flow, whose phase speed is supersonic
relative to the outer gas. Check with the r_max sensitivity; a radiation (PML / complex-mapping) condition would fix it.

**Tanh fits (local, 2026-11-10, `MODE='profiles'`).** U: R^2 0.989-0.999 up to x = 7.0 mm. rho: 0.90-0.98, p:
0.84-0.96 (Tran: > 0.97 for all). Two layers cannot follow both the jet boundary and the barrel shock, so the rho
and p fits put their inner layer on the barrel shock (x = 3 mm: R = 2.3 and 2.1 mm) while U's is at the jet
boundary (3.04 mm, half-thickness 0.064 mm). Beyond x = 7.1 mm (Mach disk / triple point) the fits break down
(R^2 for rho down to 0.71, r1 jumping); do not use tanh results there.
At x = 3 mm, St_D 0.5 the KH candidates are close to the smoothed-profile values (m = 0: 4.833 - 1.341i, m = 1:
4.253 - 1.987i) but are not converged: they sit in a string of nearby eigenvalues and fail the isolation test.

**First look, x = 3.0 mm, St_D 0.5 (local test, 2026-11-10):** one KH_jet mode per m, m = 0: alpha D = 4.707 - 1.291i,
m = 1: 4.484 - 2.287i (both with downstream group velocity); no upstream_jet wave at this St.

The group velocity comes from first-order perturbation theory with the left eigenvectors
(d alpha / d omega = y^H (dL/d omega) q / y^H R q), so no second solve is needed.
