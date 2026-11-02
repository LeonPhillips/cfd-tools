# SA, p0j 8.015 bar, p0inf 3.247 bar: meshFineNozzle (2026-10-28)

Case assumptions are the same as `../SA_p0j8.015_p0inf3.247_meshCoarse/README_case_assumptions.md`.

## Mesh
`../../mesh/srp_M2.5_meshFineNozzle.su2` (built by `mesh_job_fineNozzle.sh`).

Confirmed with the user on 2026-10-28: option B, the whole nozzle + plenum box, the name meshFineNozzle, and extra lines inside the nozzle.
- It is meshFine (0.0177 mm box at x 0.5-20 mm, r 0-9.5 mm) plus the nozzle + plenum interior at x -8.2..0.5 mm, r 0..3 mm, refined from 0.05 to 0.0177 mm.
- **835,123 cells** (meshFine has 772,018). In the nozzle/plenum census region there are 97,833 cells (34,531 before).
- Quality: min angle 6.49 deg, volume ratio <= 38.6.

## Run
`auto_fineNozzle.sh` runs unattended in tmux `srp_auto` on m3f104 (log `auto.log`). The user asked for this so the laptop could be closed. It does, in order:
1. Waits for the mesh job.
2. `0_interp`: interpolates fine f3 onto this mesh (1 iteration, CFL 1e-8).
3. `run_fineNozzle.sh` runs the stages:
   - n1: cap 0.5, 5000 iterations
   - n2: cap 1, 25000 iterations
   - n3: cap 3, 20000 iterations

   All stages use floor 0.1 and are second order.
4. Line-profile analysis into `analysis/`.

The runs never stop on a convergence test: `CONV_RESIDUAL_MINVAL= -20` and the watcher has `STOP_FLAGS=NONE` (user rule).

Analysis, in `analysis/`:
- Standard lines: fine f3 vs n3 (`grid_lines`) and n2 vs n3 (`conv_lines_n2_n3`).
- Nozzle lines (axis x -8.2..2 mm; radial line at x = 1.9 mm, r 0..3): `nozzle_lines` and `conv_nozzle_n2_n3`.

## Results (finished 2026-10-28 20:09; all exit 0, about 0.14-0.15 s/iteration)
| stage | CD end | rms[Rho] | CFL | jet kg/s |
|---|---|---|---|---|
| n1_cap0.5 | -0.12269 | -5.03 | 0.5, no drops; nozzle re-mesh transient in CD settled | 0.023718 |
| n2_cap1 | -0.12244 | -6.51 | 1, no drops | 0.023718 |
| n3_cap3 | -0.12183 | -6.88 | 3, **no drops at all** | 0.023717 |

- Fine f3 had jet mass flow 0.023721 kg/s. Imbalance is 0.03 %.
- **Effect of the nozzle refinement (fine f3 vs n3):** mean |dMach| 0.0006 on the axis in the nozzle and 0.0012 on the exit line at x = 1.9 mm.
  On the standard lines it is 0.001-0.005, about the same as the n2 -> n3 iterative change (0.0006-0.0043).
  So the change is negligible, and the 0.05 mm nozzle cells were already adequate.
