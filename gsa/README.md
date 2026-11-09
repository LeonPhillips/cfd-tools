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
