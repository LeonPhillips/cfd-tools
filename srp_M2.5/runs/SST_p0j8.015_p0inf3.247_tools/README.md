# k-omega SST mesh study, case p0j 8.015 bar / p0inf 3.247 bar (srp_M2.5)

Launched 2026-10-30 11:30 in tmux `srp_sst` on m3f104 (desktop job 61510172, ends 31 Oct 16:40). Log: `chain.log`.

## Confirmed assumptions (user, 2026-10-30)
- **The same as the SA case** (`../SA_p0j8.015_p0inf3.247_meshCoarse/README_case_assumptions.md`) except the turbulence model:
  conditions, gas, boundaries, AUSM+-up2 with MUSCL and Venkatakrishnan, first-order turbulence equations, and the meshes.
- **Model:** `KIND_TURB_MODEL= SST`, `SST_OPTIONS= V2003m`, no compressibility correction.
- **Turbulence inflow:**
  - freestream: TI 1 %, mu_t/mu 1
  - jet inlet: `MARKER_INLET_TURBULENT= ( inlet_jet, 0.02, 10.0 )`
- **Start (option B):** coarse runs from scratch, with no SA solution used anywhere.
- **Nothing stops on convergence:** `CONV_RESIDUAL_MINVAL= -20`, and the watchers only log.

## Schedule
- **meshCoarse:**
  - `1_freestream`: SA Phase 1 settings (nozzle plugged, first order for 500 iterations, then second order by 1500,
    CFL 1 -> 50), a fixed 20000 iterations.
  - `2_jet_staged`, jet on:
    - `t1_cap0.3`: first order for 500 iterations, then second order by 1000.
    - then `t2_cap0.5`, `t3_cap1` and `t4_cap3`.
    - Each stage is 50000 iterations, with floor 0.1 and a x1.01 ramp. Each starts at the previous cap (t1 at 0.1).
- **meshMedium, then meshFine, meshFineJet, meshFineJet2 and meshVeryFine:**
  - `0_interp`: SU2 restart interpolation of the previous mesh's last stage.
    - Medium comes from coarse t4.
    - FineJet comes from fine, and fineJet2 from fineJet.
    - VeryFine comes from fine.
  - Then `t1_cap0.5` (5000 iterations), `t2_cap1` (25000) and `t3_cap3` (20000).
- **Test on 2026-10-30:** 300 iterations of Phase 1 and 100 of t1, both exit 0, in `../SST_p0j8.015_p0inf3.247_meshCoarse/test/`.
  Speed: 0.055 s/iteration in Phase 1 and 0.039 s/iteration in t1 on 24 ranks, about the same as SA. Expected total about 13 h.

## Tooling (this folder)
- `setup_sst.sh` built all the configs from the SA ones.
- `run_sst_chain.sh` runs the chain and skips meshes that are already done, so it can be restarted.
- `watch_fj.sh`, `box_snap.py`, `box_change.py` and `plot_conv_roi.py` are copies of the meshFineJet tools.
  - The plot also handles rms[k] and rms[w] and puts the run name in the title.
  - They provide the probes p1-p7, the jet-box change (mean, p99 and top-1 % mean) and `convergence_latest.png`.
- Analysis written by the chain:
  - each mesh's `analysis/`: the last-two-stage line check, and SA vs SST lines
  - `../grid_compare/`: `sst_grid_lines.*`, and `box_SST_SA_0.025mm.npz` for the contour pages

## Phase 1 restart (2026-10-30 11:41)
- **First launch (11:30):** Phase 1 with the SA Phase 1 CFL (1 -> 50, x1.2 up / x0.1 down) diverged.
  - The CFL kept cycling between about 0.5 and 10-36, and rms[k] rose from -1 to +1.45 after second order switched in.
  - SU2 called MPI_Abort at iteration 2227 without printing an error.
  - The files are kept in `../SST_p0j8.015_p0inf3.247_meshCoarse/1_freestream/failed_cfl50/`.
- **Relaunched (user):** Phase 1 CFL starts at 0.1, ramps x1.01 to a **cap of 0.5**, is halved on trouble, **floor 0.1**.
  Still a fixed 20000 iterations. The rest of the chain is unchanged.
- **At iteration 4676:** past the old failure point. rms[Rho] about -3.7. The CFL is still dropping now and then
  (to 0.1-0.15), and rms[k] has crept up from -1.3 to -0.2.
