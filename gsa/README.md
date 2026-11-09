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
