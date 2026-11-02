# SA, p0j 8.015 bar, p0inf 3.247 bar: meshVeryFine (2026-10-28/29)

Case assumptions are the same as `../SA_p0j8.015_p0inf3.247_meshCoarse/README_case_assumptions.md`.

## Mesh
`../../mesh/srp_M2.5_meshVeryFine.su2` (built by `mesh_job_veryFine.sh`).

Confirmed with the user on 2026-10-28: option A, the same box, start from fineNozzle n3, stages v1-v3, and coarse included in the grid plots.
- Box at x 0.5-20 mm, r 0-9.5 mm at **0.0125 mm**. Medium / fine / very fine are 0.025 / 0.0177 / 0.0125 mm, a constant sqrt 2 ratio.
- The nozzle stays at 0.05 mm.
- 1,381,423 cells. Min angle 6.49 deg, volume ratio <= 19.6.

## Run
`auto_veryFine.sh` ran unattended in tmux `srp_auto2`: interpolation of fineNozzle n3, then `run_veryFine.sh` (v1 cap 0.5 for 5k, v2 cap 1 for 25k, v3 cap 3 for 20k iterations), then analysis into `analysis/`.
Nothing stops on convergence.

| stage | CD end | rms[Rho] | CFL drops | jet kg/s |
|---|---|---|---|---|
| v1_cap0.5 | -0.12191 | -5.39 | 0 | 0.023727 |
| v2_cap1 | -0.12179 | -6.07 | 0 | 0.023722 |
| v3_cap3 | -0.12151 | -6.60 | 0 (held 3) | 0.023722 |

## Grid study, against very fine v3
- Mean |dMach| per line: coarse 0.003-0.095, medium 0.002-0.025, fine 0.001-0.010.
- Iterative change v2 -> v3: 0.0002-0.003, so the remaining differences are from the mesh.
- Bow shock and Mach disk positions are converged (within about 2 cells).
- Jet-core peak Mach values still change by about 1 % per level.
- The shear-layer peak at r = 3 mm is not converging (the change grows from fine to very fine).
