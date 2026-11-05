# SST, SLAU2 + Venkatakrishnan-Wang, meshFine, p0j 7.0 and 7.4945 bar (Gadi)

Set up on 2026-11-05. These runs test the effect of the convective scheme: the same cases as the AUSM+UP2 meshFine
t3 solutions, rerun for 10,000 iterations with SLAU2 fluxes and the Venkatakrishnan-Wang limiter.

## Confirmed assumptions (user, 2026-11-05)

- **Cases:** SST (V2003m) on meshFine (772,018 cells). Two jet inlets, both total conditions with T0j 288.15 K:
  - **A:** p0j 7.0 bar
  - **B:** p0j 7.4945 bar

  Everything else is as in the AUSM+UP2 runs: freestream M2.5 from p0inf 3.247 bar, axisymmetric, first-order k-ω
  convection.
- **Scheme:** one combined setup:
  - `CONV_NUM_METHOD_FLOW= SLAU2`, with MUSCL
  - `SLOPE_LIMITER_FLOW= VENKATAKRISHNAN_WANG`
  - `VENKAT_LIMITER_COEFF= 0.05`, the same value. With VW it is scaled by each variable's global range.
- **Start:** that case's AUSM+UP2 meshFine t3 solution (iteration 19999), copied into `0_start/`. The 7.0 bar files
  were copied from M3 on 2026-11-05 and SHA-256 checked. The meshFine mesh hash is the same on M3 and Gadi.
- **Short CFL ramp, 10,000 iterations in total:**

  | Stage | Iterations | CFL |
  |---|---|---|
  | t1 | 1000 | 0.1 → cap 0.5 |
  | t2 | 1000 | 0.5 → cap 1 |
  | t3 | 8000 | 1 → cap 3 |

  The CFL adapt factors are as before: ( 0.5, 1.01, 0.1, cap, 0.001, 0 ).
- **Convergence:** nothing stops. `CONV_RESIDUAL_MINVAL= -20` means every stage runs its full ITER.
- **Live monitor:** every 1000 iterations, on the default srp_M2.5 lines:
  - r = 0.01, 1.3, 2.3, 3 and 5 mm, x 1–20 mm
  - the shear layer from (2.0, 2.34) to (7.84, 4.70) mm

  Probes sit 1 mm either side of the Mach disk at r = 0, 1 and 2 mm, re-placed per stage. The watcher is
  `monitor/watch_conv_gadi.sh`, which uses ParaView 6.1.1. The report is the chain figures (`plot_chain.py`), as in
  `figures/srp_M2.5/SST/p0j7.0/meshFine/`.

## Gadi resources

- One normal-queue node per case: 48 cores, 190 GB. SU2 runs on 46 ranks; the other 2 cores run the watcher.
- Walltime 0:50, so each case can charge at most 80 SU. Expected about 40–50 SU.
- A and B run at the same time.
- The test job (`test.pbs`, 0:20, at most 32 SU) is run first.

## Files

| File | What it does |
|---|---|
| `setup_slau2vw.sh` | Builds the run folders under `/scratch/os85/lp3540/su2/cases/srp_M2.5/runs/` (`TEST=1`: the test folder) |
| `run_chain_gadi.sh` | The chain driver: probes, SU2, watcher per stage, then `plot_chain.py` |
| `chain.pbs` | One case: `qsub -N <name> -v RUN=<run folder> chain.pbs` (title from `<run>/title.txt`) |
| `test.pbs` | Checks ParaView 6.1 against M3's 5.13 (`compare_pv.py`), then runs 2000 live iterations |
