<!-- Copied 2026-10-28 from the old run folder ../SA_p0j8.015_p0inf3.247/ (mesh v2), which was deleted.
     These assumptions apply to the v4 runs, with the changes listed in README.md (geometry, mesh, outlet). -->

# Case SA_p0j8.015_p0inf3.247 (srp_M2.5)

SA RANS, axisymmetric, steady. Freestream M2.5 first, then the jet is switched on.
Mesh: `../../mesh/srp_M2.5_coarse_v2.su2` (first cell 1.8e-7 m on all walls).

## Confirmed assumptions (user, 2026-10-27)
On top of the case-level assumptions in `../../README.md`:

- **Conditions:** p0inf = 3.247 bar, p0j = 8.015 bar, T0inf = T0j = 288.15 K.
- **Freestream** (isentropic, gamma 1.4, R 287.058):
  - T = 288.15/2.25 = 128.067 K, p = 324700/2.25^3.5 = 19003.9 Pa.
  - U = 567.16 m/s in -x (AOA = 180 deg), rho = 0.51694 kg/m3.
  - mu = 8.861e-6 Pa s (Sutherland), Re/m = 3.31e7.
  - Pitot pressure about 1.62 bar.
- **Turbulence:** standard SA-neg with ft2 (`SA_OPTIONS= NEGATIVE, WITHFT2`; SU2 rejects NEGATIVE alone), no compressibility correction, freestream nu_tilde = 3 nu.
- **Mesh:** remeshed so the nozzle-throat y+ is about 1 at p0j = 8.015 bar, with a first cell of 1.8e-7 m on all walls.
- **Boundaries:**
  - Freestream inlet: `MARKER_SUPERSONIC_INLET` with static T, p and U.
  - Outlet: `MARKER_OUTLET` at p = 19003.9 Pa.
  - Walls: no-slip walls adiabatic; slip walls `MARKER_EULER`; axis `MARKER_SYM`.
- **Numerics:**
  - AUSM+-up2 with MUSCL and Venkatakrishnan limiting; turbulence first order.
  - Implicit, adaptive CFL starting at 1, max 50 (min 0.5).
  - Steady. If the solution looks like it is oscillating because the flow is unsteady, stop and report to the user (no switch to dual time stepping without asking).
- **Force reference:** area = pi (9 mm)^2 = 2.54469e-4 m2, length = 18 mm.
- **Resources:** 24 MPI ranks on a Strudel desktop (m3f104, job 61510172, L4 GPU, 4 d 4 h limit), run in tmux with `run_su2.sh`; restart every 1000 iterations. The plan originally used 32 ranks on `comp`; the user changed this on 2026-10-27.

## Phases
1. **`1_freestream/`**
   - Setup: the whole domain starts at the uniform freestream (impulsive start), with `inlet_jet` as an adiabatic no-slip wall (nozzle plugged).
   - Order: first order for 500 iterations, then ramped to second order by iteration 1500.
   - Stop: the density residual has dropped **2 orders below its peak** (`REL_RMS_DENSITY <= -2`), checked once fully second order, with a cap of 20000 iterations.
2. **`2_jet/`**
   - Setup: restart from Phase 1, with `inlet_jet` set to `MARKER_INLET` total conditions (8.015 bar, 288.15 K, +x), stepped on directly. The CFL restarts at 1.
   - Order: first order for 500 iterations, then second order by iteration 1000.
   - Stop: rms[Rho] has dropped **3 orders below the jet-on peak**, with a cap of 30000 iterations.
   - Then check that the axial force and jet mass flow each stay within 1 % over the last 2000 iterations.

`test/` holds a 60 + 60 iteration check of both configs. It is not a result.

## Running
On the desktop node, from this folder: `tmux new -d -s srp_phase1 "bash run_su2.sh phase1"`
(stages: `test`, `phase1`, `phase2`). Phase 2 is started by hand after Phase 1 has been checked.
