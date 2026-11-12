# gsa: LinStab2D global stability analysis of the srp_M2.5 base flow

Work folder on M3: `/scratch/qk75/LeonPhillips/PhD/GSA/` (LinStab2D clone `LinStab2D/` at eefce30).
Goal: predict where the discrete pressure-fluctuation frequencies are, *without* being told the experimental
frequencies (user, 2026-11-09). The comparison with experiment comes after the prediction.

## Confirmed assumptions (user, 2026-11-09)

1. Experimental tone frequencies are withheld on purpose; the GSA is a blind prediction.
2. Base flow: `cases/srp_M2.5/runs/SST_p0j6.6_p0inf3.247_meshFine/2_jet_staged/t3_cap3/flow.vtu`
   (SST V2003m, AUSM+UP2, p0j 6.6 bar, p0inf 3.247 bar, meshFine, final stage).
3. Molecular viscosity only (Sutherland), no eddy viscosity in the linear operator.
4. Azimuthal wavenumbers m = 0 and m = 1.
5. Scales: length = body D = 0.018 m; velocity, Re and Ma from the freestream. Derived from the t3 config
   (p = 19003.9 Pa, T = 128.067 K, M = 2.5, air gamma 1.4, R 287.058, Sutherland 1.716e-5 @ 273.15 K, S = 110.4):
   - a = 226.865 m/s, U = 567.162 m/s, rho = 0.51693 kg/m3, mu = 8.8606e-6 Pa s
   - Re_D = rho U D / mu = 5.956e5; Ma = 2.5; Pr = 0.72; kappa = 1.4; T_0 = 128.067 K
   - omega = 2 pi f D / U, so f [Hz] = 5014.7 * omega, and St_D = omega / (2 pi)
   - LinStab2D's `sutherland_air` uses other reference constants (1.735e-5 @ 280 K) but the same S = 110.4, so the
     ratio mu/mu_inf it uses is identical to SU2's.
6. Domain and boundary conditions: match the CFD solution (details to be confirmed, see the plan).
7. Shocks: handled with LinStab2D's spatial filter (no base-flow smoothing).
8. GSA only (no resolvent).
9. Estimate memory and walltime before every run.

### Revised 2026-11-09 (user): follow the supervisor's GSA method

The supervisor's package (Micah Wong, `context/Supervisor_Package_GSA_Resolvent_FastSolver/README.md`) is the method.
This replaces items 3, 6 and 7 above:

- **Viscosity:** effective, MU = MU_Sutherland + sigma_visc * MU_t, with MU_t = mu_t(SU2 `Eddy_Viscosity`) / mu_inf.
  mu_t is frozen (only the Sutherland part is perturbed), and the total mu is used in heat conduction (Pr_eff = 0.72).
  **A: sigma_visc = 1 baseline**, with sensitivity runs at 0.2 and 0.5.
- **B: domain** truncated to x -9..+20 mm (-0.5 D..1.1 D), r 0..27 mm (1.5 D), uniform masked grid, tanh sponges at the
  edges. The nozzle exit plane is an inflow. The tunnel wall and far wake are left out.
- **BCs (method):** inflow Dirichlet row-only on all 5 variables (penalty 1e8); walls Dirichlet u, v, w and Neumann T
  (adiabatic); outlet/top Neumann on all variables (1e4) plus sponges. Axis: LinStab2D's axisymmetric treatment.
- **Shocks:** local shock smoothing of the base flow (`smooth_shock` method), not LinStab2D's filter.
- **C: operator** = LinStab2D's axisymmetric operator, with the grad(mu) terms (`dmudT.*dTdr`, `dmudT.*dTdz`,
  `d2mudT2` terms) replaced by gradients of the total mu field. Regression: identical to upstream when mu_t = 0.
- **D: alpha_rho = 0** to start; 3e-5 only if spurious density modes appear.
- **Eigen-solve:** 50 modes, tol 1e-8, p = 300, shifts with imaginary part +0.1; adjoint from the same factorisation;
  MUMPS shift-invert.

## Memory / time calibration (2026-11-09, m3f102, 16 threads, r2023b)

Cavity example GSA (direct + adjoint `eigs`, 50 modes, p = 300), refined 1-2.5x (`phase0/scale_cavity.sh`):

| DOF | nnz(L) | eigs (dir + adj) | peak RSS |
|---|---|---|---|
| 81,705 | 2.1M | 90 s | 4.2 GB |
| 182,555 | 4.9M | 322 s | 8.5 GB |
| 323,405 | 8.7M | 671 s | 23.7 GB |
| 504,255 | 13.7M | 1109 s | 38.7 GB |

Fit over all four: memory ~ DOF^1.26, time ~ DOF^1.39. Between the last two points memory grows ~DOF^1.5, which
is what the estimates use (conservative).

## Fast solver: `fast_solver/` (own MUMPS wrapper, validated 2026-11-09)

The supervisor package has the MUMPS MEX build for M3 but not the MATLAB wrapper (`LinStab2D_addon/`), so this is our
own implementation of the design described in its README (user, 2026-11-09).

- `MumpsFactor.m`: factorise once (MUMPS 5.7.3, complex unsymmetric, METIS `ICNTL(7)=5`, `CNTL(1)=0.1`, `ICNTL(14)=30`,
  rows equilibrated by their max |entry|); `solve` = A\b, `solve_ctranspose` = A'\b from the same factors
  (`ICNTL(9)=2` + conjugation). Solves pass an empty sparse matrix, so A is not copied per solve. Random-RHS check at
  factorisation.
- `make_shift_invert_solver.m`: types `'mumps'` (no silent fallback unless `allow_fallback`) and `'umfpack'`.
- `eigs_shift_invert.m`: drop-in for the two `eigs` calls; adjoint via invW (L0 - sigma I)^-H W; true residuals of
  every mode against the original L0.
- `fls_mumps_available.m`: adds `$FLS_MUMPS_DIR` and checks for `zmumpsmex`.
- MEX: `/scratch/qk75/LeonPhillips/PhD/GSA/mumps_m3/` (copied from `mumps_builds/m3`, SHA-256 checked; only system
  libgfortran/libgomp/libstdc++ needed). Set `OMP_NUM_THREADS` and `maxNumCompThreads` to the core count.
- Validation (`validate_fast_solver_cavity.m`, `run_validate.sh`; m3f102, 16 threads; output
  `GSA/fast_solver_validation/`):

  | Case | Solver | Time (direct + adjoint) | Peak RSS | vs LinStab2D's own `eigs` |
  |---|---|---|---|---|
  | cavity, 81,705 DOF | builtin `eigs` x2 | 89.7 s | | reference |
  | | umfpack wrapper | 92.2 s | | max\|dw\| 2.9e-14 / 1.6e-13 (adj), alignment 1-1e-15 |
  | | **mumps** | **41.5 s** (factorise 1.7 s, 1.46 GB) | 5.4 GB (all three) | max\|dw\| 1.0e-13 / 1.5e-13 (adj), alignment 1-1.2e-15 |
  | cavity x2.5, 504,255 DOF | builtin (scaling test) | 1109 s | 38.7 GB | |
  | | **mumps** | **335 s** (factorise 14.5 s, factors 12.4 GB) | **17.0 GB** | residuals 1.5e-12 / 2.7e-16 |

  Max eigen-residuals ||L v - w v|| / (||L||_1 ||v||): 2.1e-13 direct, 1.7e-13 adjoint (same as UMFPACK).

## Operator: `operator/` (axisymmetric LNS with variable mean viscosity, 2026-11-09)

Route C: LinStab2D's axisymmetric operator, regenerated symbolically so that grad(mu_bar) comes from the total
viscosity field.

- `derive_axy_lns.m` (Symbolic Math Toolbox, ~25 s): linearises the conservative axisymmetric compressible NS
  equations in LinStab2D's form (q = [rho, U_r, V_theta, W_z, T], p = rho T/c2, E = cv T + |u|^2/2, stress/Re,
  heat flux -(mu/c1) grad T, exp(1i m theta)) and writes two coefficient scripts (generated, don't edit):
  - `getCoeffsAxyMuT.m`: mu = mu(T), grad(mu_bar) = dmudT grad(T_bar)
  - `getCoeffsAxyVarMu.m`: mu_bar = MU(x, r) any field; grad(mu_bar) from inputs dMUdr, dMUdz; mu' = dmudT T'
    (Sutherland only, frozen eddy viscosity); heat conduction with the total MU (Pr_eff = Pr)
- `GetLinProblemVarMu.m`: GetLinProblem('axy') with a selectable coefficient script (default VarMu). dMUdx/dMUdr
  default to mesh.Dx*MU, mesh.Dy*MU.
- `BC_DirichletRowOnly.m`: the method's row-only Dirichlet (penalty 1e8).

**Checks (m3f102, `GSA/operator_dev/`):**

- Coefficients at 200 random points (`test_coeffs_vs_upstream.m`):
  - 171 of 175 coefficients equal LinStab2D's `getCoeffsLaminar.m` to round-off.
  - The 4 that differ are A0_25, A0_35, A0_45, A0_55 (relative 2e-1, 4e-2, 1e-2, 2e-3 on random inputs). With
    dmudT = 0 all 175 agree (4e-16), and the difference does not depend on m or d2mudT2.
  - So **upstream LinStab2D omits mu' x (base viscous stress) terms**: mu' div(tau_bar) in momentum and the matching
    viscous-work terms in energy. It keeps the grad(mu') and 1i m mu' terms. The supervisor's `GetOperators.m`
    includes these terms (e.g. `dmu_dT * d2u_dyy` in L25), and so does our derivation.
  - VarMu with grad(MU) = chain rule equals MuT to 2e-16. With MU + constant it equals MuT with that MU.
- Operator on example_5_Paraboloid (`test_operator_paraboloid.m`, m = 0 and 5, 10,000 points):
  - Our assembly with upstream's coefficient script equals `GetLinProblem('axy')` exactly (wiring check).
  - With dmudT = 0, MuT equals upstream to 1e-16. VarMu with chain-rule gradients equals MuT to 2e-20.
  - The full linearisation vs upstream: 6e-7 relative (Ma 0.01, so T hardly varies there).

**Decision (user, 2026-11-09): use the full linearisation** (`getCoeffsAxyVarMu`, the default in
`GetLinProblemVarMu`), including the mu' x base-stress terms that upstream LinStab2D omits. The regression rule
"identical to upstream when mu_t = 0" is replaced by the checks above: identical with dmudT = 0, and identical in the
171 coefficients that don't carry those terms.

## Base flow and grid: `baseflow/` (Phase 3, 2026-11-09)

- `export_baseflow.py` (pvbatch): SU2 nodes in x -12..23 mm, r 0..30 mm (662,880 nodes for meshFine t3), plus the
  wall_noslip segments. Output: `GSA/srp_p0j6.6_meshFine_t3/baseflow/su2_{flow,wall}.{bin,json}`.
- Geometry read from the wall (mm): nozzle exit / lip tip x = 1.876 (exit radius 2.27, lip fillet to 2.37), lip
  outer radius 2.88 from the body face x = 0.496 to the tip, body radius 9.0 (0.1 mm shoulder fillet), throat r 2.0
  at x ~ 0.
- `build_srp_case.m` + `srp_case_defaults.m` + `run_build.m` + `check_srp_case.m`:
  - Uniform grid, h = dr. The x grid is aligned with the exit plane; the face is staircased to within h. The first
    radial row is at dr/2 (LinStab2D symmetry).
  - Solid mask = body (x < 0.496, r < 9) + lip/nozzle block (x < 1.876, r < 2.88). The nozzle interior is outside
    the domain.
  - Boundary sets: `ji` (exit plane, r < 2.27: inflow), `wi` (all other masked edges: walls, with axis-aligned
    normals into the fluid, diagonal at corners), `ri` (x max, freestream inflow), `li` (x min), `ti` (r max).
  - Base flow is interpolated from the SU2 nodes (scatteredInterpolant, natural; nearest for the few points outside
    the hull) and nondimensionalised by D and freestream.
  - mu_t = SU2 `Eddy_Viscosity`/mu_inf, clipped >= 0, masked Gaussian 5 cells. MU = Sutherland (LinStab2D
    `sutherland_air`, T_0 = 128.067 K) + sigma_visc mu_t.
  - Local shock smoothing as the method, but with masked Gaussians (no leakage from the solid) and sigmas in mm.
  - Sponges as the method's GetSponge (tanh, max of three): x max (amp 5), x min (amp 1.5), r max (amp 5); width
    3 mm, transition 1 mm.
- Sizes: h = 0.1 mm gives 69,373 fluid points (347k DOF); h = 0.05 mm gives 277,434 points (1.39M DOF), built in ~30 s.
- `operator/BC_NeumannRows.m`: zero normal gradient rows (penalty 1e4) on l/r/t edges and on the masked walls
  (per-node normals).
- Check plots: `figures/newFigs/2026-11-09_srp_GSA_p0j6.6_meshFine_case-*.png`.
  - Variant A (method-like sigmas: filter 0.2, blend 0.3 mm, sensor everywhere) smooths the lip expansion and the
    start of the shear layer (Mach on the shear-layer line drops from 2.8 to ~1.4 at x = 2 mm).
  - Variant B (filter 0.1, blend 0.15 mm, sensor only for x > 3 mm) keeps the lip and shear layer as the CFD, and
    spreads the bow shock and Mach disk over ~3-4 cells without moving them.
- **Decision (user, 2026-11-09): variant B** (filter 0.1 mm, blend 0.15 mm, sensor threshold 0.15 of max |grad p|,
  sensor sigma 0.15 mm, sensor only for x > 3 mm). Now the default in `srp_case_defaults.m`; the h = 0.05 mm case
  is `case/case_h0p05_B.mat`.

## Phase 4 smoke test (2026-11-09, h = 0.1 mm, shift 2 + 0.1i, m = 0 and 1)

Driver `run/run_gsa_srp.m`; output `GSA/srp_p0j6.6_meshFine_t3/smoke_h0p1/`; figures
`figures/newFigs/2026-11-09_srp_GSA_p0j6.6_smoke-h0p1_m{0,1}_{spectrum,modes}.png`.

- **Mechanics work.**
  - 346,865 DOF; operator build ~5 s.
  - MUMPS factorisation 6-9 s (factors 4.6 / 6.7 GB for m = 0 / 1). 50 modes direct + adjoint in 150-160 s per shift.
  - Peak RSS 16.6 GB (both m in one process); wall 6 min.
  - Eigen-residuals <= 7e-16, and every direct eigenvalue pairs with a conj(adjoint) one to <= 4e-12.
  - BCs as the method: Dirichlet on 1,465 + 528 dofs, Neumann on 176 + 2,340.
- **The leading modes are numerical, not physical.**
  - Many modes are strongly "unstable" (omega_i up to 0.83) at St_D 0.15-0.47 for both m.
  - Their shapes are grid-scale (odd-even / checkerboard) oscillations: on the bow shock, in the jet core, and in the
    supersonic freestream ahead of the bow shock, where perturbations should not live.
  - 25-35 % of their energy sits in the sponges.
  - At h = 0.1 mm the variant-B shock smoothing (0.1 mm) leaves the shocks ~2 cells thick. Upstream, mu_t/mu_inf ~ 0.3,
    so Re_eff ~ 6e5 there: nothing damps grid-scale content in the inviscid supersonic regions.

## Grid check + mode diagnostics (2026-11-09, shift 2 + 0.1i)

New per-mode metrics in `run_gsa_srp.m`:
- `grid_frac`: energy at grid scale, |q - Sq|^2/|q|^2 with S a masked 3x3 [1 2 1]^2/16 filter
- `shock_frac`: energy within 3 cells of the shock sensor
- `sponge_frac`: energy where sponge > 0.5

`locate_modes.m` gives the energy per sponge and boundary set, and the location of max |p'|.
Launcher: `run/run_gsa.sh`. Output: `GSA/srp_p0j6.6_meshFine_t3/gridcheck2_h{0p1,0p05}/`.

Fixes found on the way (cfd-tools 1837142):
- Neumann rows on r max now use the inward normal (it was outward, i.e. anti-relaxing).
- MUMPS instances are now actually freed (`zmumps` needs an output argument). h = 0.05 peak RSS went from 64 to 45 GB
  with both m in one process.

**Results:**

| h | m | top omega | grid_frac | shock_frac | sponge_frac | where |
|---|---|---|---|---|---|---|
| 0.1 | 0, 1 | many up to +0.83i | 0.8-0.95 | 0.15-0.2 | 0.25-0.35 | checkerboard on the bow shock, jet core, freestream |
| 0.05 | 0 | 1.394+1.194i, ... | ~0.3 | ~0 | 0.8-1.0 | outflow sponge, max at the outflow/body-wall corner (x -9, r 9.08 mm) |
| 0.05 | 1 | 2.575+1.218i, ... | ~0.3 | ~0 | 0.95-1.0 | same |

- The h = 0.1 checkerboard shock modes are gone at h = 0.05: they were resolution artefacts.
- At h = 0.05 the window is filled by modes at the corner where the body boundary layer leaves through the outflow.
- Outflow sponge test (amp 1.5 -> 5, width 3 -> 5 mm, h 0.05, m 0; `spongetest_out5_h0p05/`): the leading
  eigenvalues barely move (1.3936+1.1941i -> 1.3941+1.1948i), and their sponge fraction drops to 0.5-0.8 as the modes
  spread. So they are a boundary-closure (corner) artefact, not something the sponge controls.
- No physical candidate yet in this window (none with low sponge and grid fractions).

## Nozzle interior + outflow at x = -27 mm (user, 2026-11-10)

- **Base flow:** `baseflow/su2L_*` (x -30..23 mm, 681,738 nodes). The export box now comes from the `GSA_BOX` env var.
- **Mask:** the bore is opened from the plenum inlet (x = -8.124 mm, r 2.879) through the converging-diverging nozzle
  to the exit. The inner wall is staircased from the mesh wall contour (min r of the wall segments per column), and
  fluid runs shorter than 3 points are closed.
- **Boundaries and sponge:**
  - New inflow set `ni` at the inlet face: Dirichlet on all variables (the method's inlet; **assumption, to be
    confirmed**). SU2 uses total conditions there.
  - Outflow at x = -27 mm with a long, gentle sponge (amp 1.5, width 12 mm, transition 4 mm).
- **Sizes:** h 0.1 gives 104,401 points (522k DOF); h 0.05 gives 417,542 points (2.09M DOF; MUMPS 31 / 45 GB for
  m 0 / 1; peak RSS 65 GB; 15-16 min per m).
- **Check plots:** `figures/newFigs/2026-11-10_srp_GSA_p0j6.6_case-nozzle-h0p05_*.png`. Nozzle Mach goes from ~0.3 in
  the plenum to sonic at the throat and supersonic at the exit, as in the CFD.
- **GSA (shift 2+0.1i, `noz_h0p1`, `noz_h0p05`):** the leading modes are still spurious.
  - They moved with the corner: max |p'| at (-27, 9.1) mm, the first rows above the body side wall at the outflow.
    93-100 % of their energy is in the x-min sponge, little of it on the wall or outflow rows themselves.
  - Nothing comes from the nozzle or its inlet (energy at `ni` ~0).
  - Test with a Dirichlet outflow (`outflow_bc`, h 0.1, m 0, `noz_outD_h0p1`): the leading eigenvalues barely change
    (1.3275+0.7080i -> 1.3334+0.7114i, 2.2513+0.6753i -> 2.2513+0.6750i). Only one mode (1.8147+0.8035i) went away.
  - So these are tied to the body-side wall boundary layer in the outflow region, not to the outflow BC type, the
    sponge strength, or the distance from the jet.

## Body-side slip wall (user, 2026-11-10)

- The r = 9 mm cylinder behind the face (x < 0.496 mm) is a slip wall (`si`): v' = 0, zero normal gradient for
  rho, u, w, T. The base-flow boundary layer within 0.5 mm of it is replaced by the values at 9.5 mm.
- Check: `figures/newFigs/2026-11-10_srp_GSA_p0j6.6_case-slip-h0p1_bodyside.png`. The sharp viscous layer is ~0.5 mm;
  the outer Mach rise up to ~10.5 mm is kept.
- **GSA (shift 2+0.1i; `slip_h0p1`, `slip05_h0p05`; h 0.05: 65 GB peak, 35 min for both m):**
  - The outflow-corner localisation is gone: sponge fractions are now 0.02-0.3, down from ~1.
  - But the leading eigenvalues are the same as in the no-slip runs (h 0.05, m 0: 1.4326+1.199i, 1.1735+0.746i,
    2.8138+0.678i, 1.9958+0.628i, ...).
  - Their eigenvectors are now almost pure grid-scale content (grid_frac 0.95-0.99), at:
    - the bow shock (11.7, 8.8)
    - the triple point (7.3, 4.5)
    - the axis barrel-shock crossing (4.8, 0)
    - the Mach disk (8.2-8.5, 0.7-1.8)
    - the nozzle converging wall (-1.8, 2.6)
  - Earlier their maxima sat at the outflow corner. They are grid-dependent (different at h 0.1), so they are
    numerical: odd-even modes of the central differences at shocks and steep gradients. The eigenvector maximum
    follows whichever boundary closure lets it grow most.

## Density damping alpha_rho (user, 2026-11-10)

- Implemented in `GetLinProblemVarMu` (`opts.alpha_rho`, run option `'alpha_rho'`, default 0): the method's term,
  applied in both x and r, on rho' only: `D_rho = -alpha_rho (hx^3 H^-1 Dxx' H Dxx + hr^3 H^-1 Drr' H Drr)`, H = W r,
  Drr with the scalar parity for m. `operator/test_alpha_rho.m` (h 0.1): only the rho-rho block changes, H-symmetric
  to 1e-14, negative definite, +4 % nnz. Damping rate (3e-5, h 0.1): checkerboard 0.25, 2 mm wave 6e-5.
- **GSA, h 0.05, m 0, shift 2+0.1i** (`arho3em5_h0p05`, `arho1em4_h0p05`; 59 GB peak, 18 min each):

  | slip05 (alpha 0) | alpha 3e-5 | alpha 1e-4 |
  |---|---|---|
  | 1.433+1.199i | 1.335+1.109i | 1.293+1.106i |
  | 1.173+0.746i | 0.991+0.669i | 1.076+0.394i |
  | 2.814+0.678i | 2.848+0.616i | 2.857+0.298i |
  | 1.996+0.628i | 1.967+0.579i | 1.963+0.631i |
  | - | 3.387+0.483i | 3.212+0.741i |

  (rows paired by frequency and location, not exactly).
- The leading modes survive with grid_frac 0.94-0.99. `locate_modes` (now with energy by variable) shows they are
  **density-dominated** (rho 45-92 % of the energy, almost all of it grid-scale), at the same places (triple point,
  axis crossing, Mach disk; with 1e-4 the leading one moves to the bow shock on the axis, (13.8, 0.3) mm).
- Reading: the damping is the right kind but too weak against the production at the shocks. Both scale as 1/h
  (production ~ jump / shock thickness, thickness ~ 2h), so neither refinement nor alpha_rho of this size wins.
  Next candidate: wider shock smoothing (so the base-flow shocks span many cells), with alpha_rho kept.

## Wider shock smoothing (user, 2026-11-10)

- `baseflow/smoothing_check.m` (no solve): variant B already covers every hotspot (blend weight 1), but leaves the
  shocks sharp: max drho/rho per cell 0.16-0.28 at h 0.05. Threshold 0.05 with filter 0.2 / blend 0.25 mm (W2):
  0.07-0.14; filter 0.3 / blend 0.35 mm (W3): 0.04-0.09. The sharpest points left (0.3/cell) are in the jet shear
  layer at x 5.6-6.3, r 3.3 mm (a density jump, not a shock); threshold 0.02 would smooth it too, so it was not used.
  Plot: `figures/newFigs/2026-11-10_srp_GSA_p0j6.6_smoothing-check_h0p05.png`.
- **GSA, h 0.05, m 0, shift 2+0.1i, alpha_rho 3e-5** (`smW2_h0p05`, `smW3_h0p05`; 59 GB, 19 min each):
  - W2: 1.997+1.281i, 1.080+0.940i, 2.344+0.386i, 1.120+0.318i, 1.848+0.235i, then decaying.
  - W3: 3.048+1.138i, 2.354+0.460i, 2.242+0.459i, 1.667+0.422i, 0.906+0.264i, 3.315+0.167i, then decaying.
  - The off-axis shock modes (triple point, off-axis bow shock, Mach disk flank) are gone, and grid_frac of modes
    2-8 falls to 0.75-0.82 with W3. The eigenvalues still change with the filter, so none is physical yet.
  - What is left (W3; plots `figures/newFigs/2026-11-10_srp_GSA_p0j6.6_smW3-h0p05-m0_{modes,where,spectrum}.png`):
    - The leading mode (grid_frac 0.99, mostly u, v) sits at the bow shock near the axis, (13.5, 1) mm, with
      grid-scale lobes filling the supersonic freestream between the bow shock and the inflow sponge. Physically
      nothing can travel upstream there; central differences let 2h waves do it (negative group velocity), so the
      bow shock and the inflow boundary trap them.
    - Modes 2-5, 7, 8 peak on the axis at the barrel-shock crossing (4.73, 0.03), where the base flow is already
      smooth (0.04/cell): likely the axis treatment (1/r terms at r = h/2) rather than sharpness.
    - With W2, three modes peak on the axis inside the nozzle at (1.38, 0.03), which is unsmoothed (x <= 3 mm).

## LinStab2D bug: (a)symmetric radial derivatives behind the body (found 2026-11-12)

- `AddDW2mesh` (LinStab2D eefce30) decides whether a radial fluid segment touches the axis from the first point
  of the whole column (`yl(1)`), not of the segment. On our masked mesh:
  - columns without the nozzle bore (x < -8.1 mm) start at r = r_body, so `Dy_symm`, `Dy_asymm`, `D2y_symm`,
    `D2y_asymm` stay **zero** for every point above the body;
  - bore columns (-8.1 < x < 0.5 mm) start on the axis, so the segment above the body gets an axis-reflection
    closure at r = r_body.
- The axisymmetric operator uses only these matrices for r derivatives (`CreateDiffMatrices_Axy`), so **every GSA
  run before this date had no (or wrong) radial derivatives in the region x < 0.5 mm, r > 9 mm** (outflow and
  body-side region). The outflow-corner and body-side boundary-layer modes found earlier were probably this.
- Fix: `fix_offaxis_symm` in `build_srp_case` (commit 7c13d1c) gives off-axis segments the plain `Dy`/`D2y`
  rows and rebuilds the mixed derivatives. `baseflow/test_stretched_grid.m` (Gadi, h 0.2/0.1): all radial
  derivatives now match analytic ones (uniform 2e-10 to 5e-4, stretched 2e-6 to 2e-5); before, relative error ~1.

## Gadi set-up and grid study (2026-11-12)

- Gadi: `/scratch/os85/lp3540/su2/GSA/` with `LinStab2D` (eefce30), `mumps_gadi` (supervisor's Gadi build of
  `zmumpsmex`), and the `su2L` base flow (SHA-256 checked). Driver `run/run_gsa_gadi.pbs` (TASK validate | test |
  gsa, PRESET from `run/presets/`). Cavity check on Gadi: MUMPS matches builtin eigs to 1.2e-13.
- Stretched grid (`p.grid.h_box_mm`): meshMedium / meshFine spacing in the SU2 box (x 0.5-20, r 0-9.5 mm), 0.05 mm
  outside, <= 5 % growth, exit plane on a node, first row at h_box/2; built with `DeformMesh` (x then r map).
- All grid-study runs share one base flow built on a uniform 0.0177 mm helper grid (W3 smoothing, mu_t filter
  0.25 mm) and interpolated linearly to each grid; alpha_rho 3e-5 with the local spacing.
