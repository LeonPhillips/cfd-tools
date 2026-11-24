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
3. Scales (changed by the user on 2026-11-13): length **D_e = 4.6 mm**, the CFD nozzle exit diameter (`S.D_ref`), so
   alpha, r and x are in units of D_e; a_inf = 226.865 m/s, rho_inf = 0.51693 kg/m3, p by rho_inf a_inf^2 (free-stream
   p = 1/gamma). **St keeps the body diameter:** St_D = f D_body / U_inf, D_body = 18 mm (`S.D_body`),
   U_inf = 567.162 m/s, so f [Hz] = 31,509 St_D and omega (code units) = 2 pi f D_e / a_inf = 2 pi St_D (D_e/D_body) M_inf.
   The eigenvalue window keeps its physical size: |alpha_r| < 5 omega + 5.11 and |alpha_i| < 12.8 (were 20 and 50 per
   18 mm). Results before 2026-11-13 were per 18 mm: multiply their alpha by 4.6/18 = 0.2556 to compare.
   Note: the exported wall gives an inner lip radius of 2.389 mm (diameter 4.78 mm), 4 % more than 4.6 mm; D_e is
   only the reference length, so this does not change the physics.
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
   - `'uniformp'` (user, 2026-11-17): uniform mean pressure p̄ (median RANS p from r1 + 0.3 mm to r2 − 0.3 mm, the
     pressure the jet shear layer sees), U and the RANS temperature θ = pν (= RT/a∞²) as two-tanh fits, ν = θ/p̄.
     A parallel flow with no mean radial velocity can only carry a uniform p̄; the RANS core pressure (expansion,
     barrel shock) is dropped. This is the usual choice in local jet LSA (Michalke 1984, Mancinelli 2021, Maia 2024).
   - `'cb'`: as `'uniformp'` but θ from Crocco–Busemann with one T0 (288.15 K, jet and free stream):
     T/T∞ = 2.25 − 0.2 U². Within 7–9 % of the RANS T (median) but 30–70 % too hot just inside r1.
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
matlab -batch "CASE='freejet'; MODE='single'; OUT_SUBDIR='single_freeJet'; NWORKERS=20; lsa_local"
```

### Free jet (`CASE='freejet'`, user 2026-11-24)

The srp_M2.5 nozzle without counterflow (`cases/srp_M2.5_freeJet`, p0j 6.6 bar into still air at 1 atm, 288.15 K).
The base flow is the fine-mesh final stage `t4_cap3`, exported with `gsa/baseflow/export_baseflow.py` to
`LSA/baseflow/srp_M2.5_freeJet_meshFine_t4/su2L_*`. Confirmed assumptions:

- Station x = 3.0 mm (same frame as the counterflow case). Lengths and alpha by D_e = 4.6 mm; velocity, density and
  pressure by the still ambient (a = 340.30 m/s, rho = 1.2250 kg/m3, ambient p = 1/gamma).
- **St_j = f D_j/U_j**, ideally expanded jet: M_j = 1.882 (isentropic from NPR 6.514), T0 = 288.15 K, U_j = 489.9 m/s,
  D_j = 4.952 mm from the throat D* = D_e/sqrt(1.3225) = 4.0 mm. omega = 2 pi St_j (U_j/D_j) D_e/a = 8.40 St_j.
  Default scan St_j 0.005-0.455 in 0.005 steps (up to 45 kHz); `SINGLE_ST` overrides.
- r_max fixed at 20 mm (`FJ_RMAX` overrides; a non-default value is tagged `_rmax<mm>`), Dirichlet there.
- Tran eqs. 3.1-3.3 two-tanh fits with the far-field constant **fixed** at ambient (U 0, rho 1, p 1/gamma), as Tran
  do. r1 = the U-fit layer at the jet edge; r2 (second grid cluster, and `KH_outer` in `label_modes`) = the rho-fit
  layer inside the jet, the barrel shock (so here r2 < r1 is possible).
- Everything else as the counterflow `single` (inviscid, N 250/300, filters, labels). Output tag `_freeJet`.

| MODE | What it does | Output (in `LSA/<MODE>/`) |
|---|---|---|
| `selfcheck` | Thin tanh jet (M_j 1.5, cold and nu_j = 1.5) against the exact compressible vortex sheet, R/theta = 20, 40, 80 | `selfcheck.csv`, PASSED/FAILED line in the log |
| `single` | One station and St_D (`SINGLE_X`, `SINGLE_ST`; default 3.0 mm, 1.271; optional `SINGLE_N`, default 250, N_check = N + 50), m = 0 and 1: full spectrum, labels, eigenfunctions. **St scan (2026-11-18):** a vector of St (set `S.single_St` near the top of `lsa_local.m`, or `SINGLE_ST`) solves every St at the one station (parfor over St with `NWORKERS`) into one file. `SAVE_V`: `'converged'` (default) saves eigenvectors only of converged, non-continuum modes; `'all'` saves every one (see below) | `x_<x>_St_<St>[_N<N>][_<PROFILE>]_lsa_single.mat`, `..._lsa_single_spectrum.png` (`_N<N>` only when N ≠ 250; `_<PROFILE>` only when PROFILE ≠ tanh) |
| `sensitivity` | x = 2.5, 4, 6, 7.5 mm, St_D 0.25, 1, 2: N 150/250/350, r_max 2.5 mm smaller, smoothing 0.025/0.1 mm | `lsa_sensitivity.csv` |
| `scan` | All stations x St x m | `lsa_scan.mat`, `lsa_scan_m0.png`, `lsa_scan_m1.png`, `lsa_scan_radii.png` |
| `thickness` | Test T1 (KH û shape against jet-layer thickness, user 2026-11-17): analytic Tran two-tanh base flows, R/θ 5-100, families b / cross / a (below) | `T1_kh_modes.csv`, `T1_thickness.mat`, `<date>_lsa_T1_*.png`, `viewer/T1_<fam>_Rtheta<NNN|ref>_lsa_single.mat` |

**St scan file (`MODE='single'` with several St; cfd-tools bcbbdaf, 07860ea).** `out.sol{iSt, im}` (row = St, column = m;
a one-St file is 1 x nm, as before), named `x_<x>_St_<min>-<max>_n<nSt>[tags]_lsa_single.mat`, with one
`..._spectrum.png` (converged non-continuum eigenvalues coloured by St, and the most unstable KH_jet against St).
Every eigenvalue, label, cg, conv and peak radius is saved for every mode; `sol.labid` (uint8, 1-5 = KH_jet, KH_outer,
upstream_jet, other, continuum) repeats the labels for fast reading; `sol.Vidx` lists which modes the columns of
`sol.V` are. With `SAVE_V='converged'` the |p|/|u| arrays `pabs`/`uabs` are not saved. Size (x = 3 mm, N 250, m 0/1):
about 3.2 MB per St with `'converged'` (test of 3 St: 9.5 MB, 83 s with 3 workers), about 31 MB per St with `'all'`.
Opened by LSA Viewer v1.0.

Every mode except `selfcheck` and `thickness` also writes `lsa_profiles.png` (raw and smoothed profiles); all write
`lsa_local.log`. Optional `OUT_SUBDIR` names the output folder under `OUT_DIR` (default: the mode name).

### MODE='michalke' (analytic base flow, user 2026-11-24)

Michalke (1971, Z. Flugwiss. 19; NASA TM 75190) base flow, with **no RANS input**:
- velocity: profile 2 (eq. 2.8), U/U1 = ½{1 + tanh[b(R/r − r/R)]}, U(R) = U1/2;
- temperature: Crocco–Busemann (eq. 2.5), T/T1 = T* + (1 − T*)y + ½(γ−1)M²y(1 − y), with T* = T0/T1 and M = U1/a1;
- uniform pressure, so ρ ∝ 1/T (eq. 2.2);
- b from R/θ, with θ the density-weighted momentum thickness (eq. 2.6; `MK_THETA_DEF='vel'` for the velocity θ),
  solved numerically (eq. 2.11's closed form holds only for T* = 1).

`MK_PART='validate'` runs a unit jet (T* = 1, R/θ = `MK_VAL_RTHETA`, default 6.25; M = `MK_VAL_M`, default 0.05, 0.4,
0.8 and 1.2, with M = 0 run as 0.05 because velocities are scaled by a sound speed) against the peak −α_iθ read
from Michalke's Figs 10–11 (library note); it passes if every peak is within 0.003. It also prints the solved b
against eq. 2.11.

`MK_PART='freejet'` (with `CASE='freejet'`) is the ideally expanded free jet (confirmed, user 2026-11-24): U1 = U_j
(489.9 m/s), M = M_j = 1.882, T* = 1 + 0.2M_j² = 1.708 (unheated), R = D_j/2 = 2.476 mm, p = p_amb, one St_j scan
(the free-jet list) per R/θ in `MK_RTHETA` (default 6.25, 12.5, 25 and 50). There is no x: the flow is parallel.
The KH pick is the most unstable converged `KH_jet` mode with α_r > 0 and c < 1.5 U1 (drops Michalke's long-wave
"mode II"). Writes one Viewer file per R/θ, `michalke_freeJet_KH_vs_Stj.{csv,png}` (with the two-tanh RANS-fit KH at
x = 3 and 4 mm overlaid, for reference only) and `michalke_freeJet_profiles.png`.

    matlab -batch "MODE='michalke'; MK_PART='validate'; OUT_SUBDIR='michalke_check'; NWORKERS=20; lsa_local"
    matlab -batch "MODE='michalke'; MK_PART='freejet'; CASE='freejet'; OUT_SUBDIR='single_freeJet_michalke'; NWORKERS=20; lsa_local"

### MODE='thickness' (test T1)

Builds each base flow as a `P` struct in the `radial_profile` format, from Tran two-tanh fits (the same as
`MODE='single'` with `PROFILE='tanh'`: exact derivatives, r1/r2 from the U fit), and solves it with the single-mode
path (`solve_case`, `mark_converged`, `label_modes`) and the single-mode defaults (N 250/300, grid, r_max, filters).
The KH mode is the most unstable converged `KH_jet` with alpha_r > 0; eigenvectors normalised to max|p̂| = 1.
One override: the `KH_jet` window (|û| peak within `sl_window_mm` of r1) is widened to max(0.6 mm, 4θ), because a
thick layer's |û| is double-peaked and its larger lobe can sit a few θ from r1 (family a, R/θ 8, m = 0 was labelled
`other` with the 0.6 mm window).
No RANS data are read.

- **b**: the stored fits of the x = 3 mm single run (`THICK_REF`, default `<OUT_DIR>/single/x_3.0_St_1.27_lsa_single.mat`).
  Only the inner-layer gradient parameter b1 of the U fit changes, to give R/θ = 5, 8, 12, 20, 40, 100; b1 of the ρ
  and p fits is scaled by the same factor. Amplitudes, radii, the outer layer and ω (St_D 1.271) unchanged. The
  unmodified fit is included and must reproduce the single-run KH eigenvalues.
- **cross**: b with U = U_j T(r; R1, b1) (outer stream at rest, U_j kept).
- **a**: one jet layer, Edgington-Mitchell 2020 regime: M_j 1.09, T0 = T∞, U₂ = 0, flat p = 1/γ,
  ρ_j/ρ∞ = 1 + 0.2 M_j², ρ layer as wide as U's, R = 0.5 D_e, ω = 2π · 0.67 · U_j (St_j on D_e). r2 = 3R is a dummy
  outer radius for `make_grid` only (amplitude 0).
- θ = ∫ f(1 − f) dr of the jet layer alone (f = T(r; R1, b1)), trapz on 200,001 points; b1 from `fzero`.
- Per mode: α, c_r/U_j, (c_r − U₂)/(U_j − U₂), α_i/α_r, α_rθ, |û| layer/inner/outer-lobe peaks over the max, the
  second peak (largest other local maximum separated by a dip), interior minimum < 0.5 × smaller peak, a visible
  double peak (that, and second peak ≥ 0.25), min-|û| radius, phase jumps, r_c (U = c_r) and the peak of
  |iαv̂ − dû/dr|. Lengths in D_e.
- Options: `THICK_FAMS` (cell), `THICK_RTHETA`, `THICK_M`, `FIG_DIR` (figure folder; default the output folder).
- `THICK_VIEWER` (default true, added 2026-11-18): also save one `MODE='single'`-format file per base flow in
  `<out_dir>/viewer/` (every eigenvalue and eigenvector of both m, `prof`, `S`), so the LSA Viewer can open it
  (about 30 MB each for b and cross, 10 MB for a). The dotted "RANS" lines are the x = 3 mm RANS line for b and cross
  and the analytic profile itself for a; the viewer's 2D field is always the x = 3 mm SRP RANS (irrelevant for a).
  `S.thick_info` holds the family, R/θ, θ, R, U_j, U₂, ω and the b factor.

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

**First look, x = 3.0 mm, St_D 0.5 (local test, 2026-11-10, per 18 mm):** one KH_jet mode per m, m = 0: alpha D = 4.707 - 1.291i,
m = 1: 4.484 - 2.287i (both with downstream group velocity); no upstream_jet wave at this St.

The group velocity comes from first-order perturbation theory with the left eigenvectors
(d alpha / d omega = y^H (dL/d omega) q / y^H R q), so no second solve is needed.
