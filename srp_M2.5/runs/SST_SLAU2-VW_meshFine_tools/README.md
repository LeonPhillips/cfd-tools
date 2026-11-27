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

## p0j 6.6 bar (user, 2026-11-06)

- **Case C:** p0j 6.6 bar (660000 Pa), T0j 288.15 K. Its AUSM+UP2 meshFine t3 config differs from the 7.0 bar one only
  in the jet inlet.
- **Start:** the 6.6 bar AUSM+UP2 meshFine t3 (iteration 19999). It was copied from M3 on 2026-11-06 and SHA-256
  checked.
- **Same process as A and B:** t1 1000 (cap 0.5), t2 1000 (cap 1), t3 8000, t4 40,000 and t5 50,000 (cap 3), so
  100,000 iterations. All five stages run in one job (`CASES=6.6 FULL=1 bash setup_slau2vw.sh`).
- **Faster:** one normalsr node (Sapphire Rapids, 104 cores, 500 GB): 102 SU2 ranks, `chain_sr.pbs`.
  - First a 2000-iteration test at 104 cores (`TEST=1 TEST_P=6.6`, 0:12, at most 42 SU).
  - Full run: 2:00 walltime, at most 416 SU. Expected 75–95 min and 260–330 SU, against 3.3 h and 311 SU on 46 ranks.

## CT sweep pilot (user, 2026-11-25)

- **Conditions:** from `cases/srp_M2.5/CT_sweep_conditions.md` on Market (Schlieren1 file names). The file-name
  CT, CTa and CTd are labels only: the RANS results are not matched to them (the line-loss estimates behind CTa/CTd
  are likely too low).
- **Boundary conditions:** p0inf = tunnel p0 absolute (gauge + 1.01325 bar). p0j = nozzle p01 after losses,
  0.81895 × target + 0.03094 bar. T0inf = T0j = 288.15 K. Freestream p = p0inf / 2.25^3.5 at M 2.5.

  | Label | p0inf | p0j | Freestream p | rho | Re/m |
  |---|---|---|---|---|---|
  | CT 0.1626 | 6.913 bar (5.90 g) | 4.535 bar (target 5.50) | 40,461.6 Pa | 1.10062 kg/m³ | 7.05e7 |
  | CT 0.5202 | 3.013 bar (2.00 g) | 7.402 bar (target 9.00) | 17,635.8 Pa | 0.47972 kg/m³ | 3.07e7 |

- **Start:** the converged 6.6 bar / 3.247 bar SLAU2-VW t5 solution (cumulative iteration 99,999), copied into
  `0_start/`. Everything else is the 6.6 bar setup: t1–t5 (100k), `chain_sr.pbs`.
- **Tests first:** `TEST=1 bash setup_ct_pilot.sh` builds `..._test` folders with t1 (1000, cap 0.5) and t2 (1000,
  cap 1); run each with `qsub -l walltime=00:15:00` (at most 52 SU). Then `bash setup_ct_pilot.sh` for the full runs.

- **Superseded (user, 2026-11-26):** the sweep now uses the purple fit (below). These runs are kept for reference.

## CT sweep (user, 2026-11-26)

- **p0j = 0.887732 × target − 0.044750 bar abs**: the purple line ("p0 (direct) fit, slope 0.888") in
  `Leon_pitotTube_data/pitot_pdrop_p0_vs_target.fig`. p0inf as above (gauge + 1.01325 bar).
- **Runs:** `/scratch/os85/lp3540/su2/cases/srp_M2.5/runs/CT-sweep/CT<label>_p0j<P>_p0inf<Q>/`, with a README there.
- **Build:** `bash setup_ct_sweep.sh CT0.1626=<start stage folder> ...`. The conditions come from
  `ct_sweep_matrix.txt` (label, tunnel gauge, nozzle target). Each case restarts from the converged solution given
  (the nearest case), then t1 1k, t2 1k and t3 8k (CFL 0.1 → 0.5, 0.5 → 1, 1 → 3).
- **Auto-stop (user, 2026-11-27: "Auto-stop, cap 250k"; limits 0.02 mm / 0.01):** run each case with
  `qsub -N <name> -l walltime=04:30:00 -v RUN=<case folder>,AUTOSTOP=1 chain_sr.pbs`. After t3, `run_chain_gadi.sh`
  adds 10k chunks (`t4_c01`, `t4_c02`, ...: CFL 1 → 3 from the previous restart) and, from 50k cumulative iterations
  on, runs `monitor/autostop_check.py` after each chunk. It stops when, over the last 10k iterations:
  - the Mach disk x and the bow shock x (axis, 1 D, 2 D) each move by less than 0.02 mm, and
  - the line mean |dMach| per 1000 iterations is below 0.01 on every line in every snapshot.

  The stagnation point is not used (it wanders in weak-jet cases). A PASS writes `<run>/AUTOSTOP_PASS.txt`. At 250k
  without a PASS it writes `AUTOSTOP_CAP.txt`, flags the case and exits 3. The first proposal (0.01 mm / 0.002) was
  dropped: replayed on finished runs, even the converged 6.6 bar run never passed it.
  `python autostop_check.py --replay <run>` replays the test on a finished chain (`TOL_MM`, `TOL_DMACH`, `WINDOW`
  override the limits).
- **Plan (2026-11-27):** the 29 cases run in three waves of about 10, each case restarting from the nearest finished
  case (combined relative change in p0inf and p0j).
- **Deferred start:** `bash setup_ct_sweep.sh CT0.1931=@CT0.1808` writes `0_start/PARENT.txt` instead of copying a
  restart. `run_chain_gadi.sh` copies the parent's final restart.dat and flow.vtu when the job starts, and stops if
  the parent has no `AUTOSTOP_PASS.txt` / `AUTOSTOP_CAP.txt`. Submit with `-W depend=afterany:<parent job>`.

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
| `chain_sr.pbs` | One case on one normalsr node (104 cores, 102 SU2 ranks), 2:00: `qsub -N <name> -v RUN=<run folder> chain_sr.pbs` |
| `setup_ct_pilot.sh` | CT sweep pilot folders, restarted from the 6.6 bar t5 (`TEST=1`: t1 + t2 test folders) |
| `setup_ct_sweep.sh`, `ct_sweep_matrix.txt` | CT sweep folders in `runs/CT-sweep/` (purple-fit p0j), each from a given start solution |
| `add_stage.sh` | Adds a continuation stage from the previous stage's restart |
| `test.pbs` | Checks ParaView 6.1 against M3's 5.13 (`compare_pv.py`), then runs 2000 live iterations |
